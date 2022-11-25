# -*- coding: utf-8 -*-
"""
Created on Tue Apr  7 14:49:02 2020

@author: Mark
"""
import numpy as np
import matplotlib.pyplot as plt

try:
    from ..model.shapes import Sphere
    from ..model.generate_angle import AngleDist
    from ..model.rotation import rotate
    from ..model.scatterer import Scatterer
    from ..model.electron import Electron
    from ..model.source import Source

except ImportError as e:
    if str(e) == "attempted relative import with no known parent package":
        import os
        import pathlib
        os.chdir(pathlib.Path(__file__).parent.parent)

        from model.shapes import Sphere
        from model.generate_angle import AngleDist
        from model.rotation import rotate
        from model.scatterer import Scatterer
        from model.electron import Electron
        from model.source import Source
    else:
        raise
#%%

class Simulation():
    """The simulation class is for running a Monte Carlo simulation of
    many electrons travelling through a scattering medium. The simulation
    is run inside a 3D volume, called the Boundary (it is a sphere).
    The scattering medium, and the source of the electrons is a volume
    called the Source (it is a cylinder at the center of the Boundary).
    Electrons are generated at random positions within the Source. They
    have an initial kinetic energy, and are given a random initial
    direction. A distance is randomly chosen from an exponential
    distribution, which is defined by the Total Cross Section (total_xsect).
    Electrons can experience two kinds of scattering event: 1. elastic
    and 2. inelastic scattering. Both elastic and inelastic processes have
    their own corss sections. The total cross section is the sum of elastic
    and inelastic. After choosing a random distance from the exponential
    distribtuion, it is randomly chosen whether the collision is elastic
    or inelastic (this is based on the relative cross sections of the two
    processes). Then, depending on the scatteing type, a random angle is
    chosen from an angular spread distribution. Based on the chosen angle,
    a new velocity is calculated (note: if inelastic scttering, the kinetic
    energy decreases each event).
    The electron's new position is determined from the random distance
    divided by the previous velocity. If the new position is inside the
    volume of the source, then the move is accepted. The position and
    velocity are updated. A counter for number of times elastically
    or inelastically scattered is incremented.
    If the electron's new position is outside of the Source, then the move
    is rejected, and the electron is projected onto the Boundary by
    extending its prevous velocity to the Boundary.
    When an electron hts the boundary, its simulation is over, and the
    simultation moves on to calculate the next electron. If the electron
    does not intersect the Boundary ofter a user defined number of
    inelastic scattering events, then the simulation of the electon is over
    and the simulation moves to a new electron.
    This process is repeated for as many electrons as desired.
    """

    def __init__(self, **kwargs):

        if 'boundary_shape' in kwargs.keys():
            self.boundary = kwargs['boundary_shape']
        elif 'r' in kwargs.keys():
            self.boundary = Sphere(kwargs['r'])
        self.intersections = []
        self.intersected = False
        self.collected = []
        self.depth = 0 # these parameters are used if one wants to generate
        # electrons piece-wise inside the Disc. Then the depth inside the
        # Disc where the electrons are generated is given by self.depth
        # and the thickness of the piece is given by self.height
        self.height = 0.01
        self.n_event_cutoff = 50 # max number of inel events per electron
        self.parameters = {}

    def addSource(self,r: float,h: float, **kwargs):
        """Generate a source of electrons.

        A Source is a shape that emits electrons. It represents a volume in
        which the initial coordinates of the electrons are generated.
        Parameters
        ----------
        r : FLOAT
            The radius of a cylinder in nm.
        h : TYPE
            Thickness of a cylinder in nm.
        **kwargs:
            angle_distribution: STR
                The initial angular distribution of the electrons
                emitted from the source.
                May be one of ['Lambert', 'Sphere', 'Constant']
            initial_KE: float
                The initial kinetic ennergy of the electrons in eV.

        Returns
        -------
        None.
        """
        self.source = Source(r,h, **kwargs)

    def addScatteringMedium(self, shape):
        self.scatteringMedium = shape

    def addScatterer(self, d:float, inel_factor:float, inel_exp:float,
                     el_factor:float, el_exp:float, Z:int, **kwargs):
        """ Add a scatterer to the scattering medium.

        The scatterer represents a medium containing atoms or molecules
        that can scatter electrons. The Scatterer has a density in atoms/nm^3,
        it has an inelastic and elastic scattering cross section factor.
        The factor is used to approximate the cross section as a function of
        kinetic energy. The factors need to be determined by fitting data from
        NIST to a 1/E^x function. The scatterer also has two angular spread
        functions: one for elastic and one for inelastic scattering. The
        parameters el_angle and inel_angle represend the FWHM of the angular
        spread function (in radians).

        Parameters
        ----------
        d : FLOAT
            Atomic density of the scattering medium in atoms/nm^3.
        inel_factor : FLOAT
            Pre-factor used to calculate the inelastic scattering cross
            section.
        inel_exp : FLOAT
            Exponential factor used to calculate the inelastic scattering cross
            section.
        el_factor : FLOAT
            Pre-factor used to calculate the elastic scattering cross section.
        el_exp : FLOAT
            Exponential factor used to calculate the inelastic scattering cross
            section.
        Z : INT
            Atomic number of the scattering medium's atoms.
        **kwargs:
            Can be 'loss_function'. The value can be a string of a filename
            for a csv with two columns (energy loss, probability).
            It can also be a list of lists, where the first element is a list
            of energy loss values and the second element is a list of
            probabilities.

        Returns
        -------
        None.
        """
        self.scatterer = Scatterer(d, inel_factor, inel_exp,
                                   el_factor,el_exp, Z, **kwargs)
        self.scatterer.setXSect(self.source.initial_KE)

    def createElectron(self):
        """Create an electron.

        This step is performed at the beginning of each electron simulation.
        It generates an Electron object, where initial position is inside the
        Source. The electron carries all of its information
        (i.e. position, velocity, time, n_times elastic and inelastically
        scattered) in an arrary called e.vector.
        """
        self.e = Electron(self.source) # instantiate electron
        self.results = np.array([np.copy(self.e.vector)]) # create place to
        # keep the results

        """Update the cross sections of the scatterer using the electron's
        new kinetic energy."""
        self.scatterer.setXSect(self.e.kinetic_energy)

    def _step(self):
        """Simulate one step.

        This function simulates one 'step' in the simulation of an
        electron's path. It calls a method from the Scatterer, called Scatter()
        which determines the kind of scattering event (elastic or inelastic),
        the distance the electron travels before the scattering event, the new
        direction of the electron relative to its old direction, in polar coord
        (radians).
        """
        """First get the distance electron travels before being scattered, as
        well as the kind of scattering event and the angle of scattering."""
        kind, d, theta, phi = self.scatterer.Scatter()
        """The new time is the distance travelled, divided by the length of
        the velocity vector."""
        velocity = self.e.vector[3:6]
        position = self.e.vector[0:3].copy()
        time = d / np.linalg.norm(velocity)
        old_position = position.copy()

        """Determine the new position."""
        new_position = (position.copy() + (velocity * time))

        if self.scatteringMedium.inside(new_position):
            """Check if next potential scattering event is inside or outside
            scattering medium. If it is still inside the scattering medium,
            then the electron is scattered once again."""
            self.e.vector[0:3] = new_position
            """Determine the new velocity after scattering.
            The velocity vector needs to be transformed into a different
            coordinate system to perform the rotation."""
            new_velocity = rotate(velocity, theta, phi)
            self.e.vector[3:6] = new_velocity
            """Check the kind of scattering (elastic or inelastic) and increment
            the count accordingly."""
            if kind == 'elastic':
                """Icrement the elastic scattering count."""
                self.e.vector[8]+=1
            else:
                """Increment inelastic scattering count."""
                self.e.vector[7]+=1
                delta_E = self.scatterer.getDeltaKE()
                self.e.updateKineticEnergy(delta_E)
                if self.e.kinetic_energy > 0:
                    self.scatterer.setXSect(self.e.kinetic_energy)
                else:
                    self.e.kinetic_energy = 0
                    self.e.vector[3:6] = np.array([0,0,0])
            """Determine the new time."""
            self.e.vector[6] += time
            self.e.pathlength += d

        else:
            """This condition is run, if the electron's new position
            is not inside the scattering medium. Then the scattering event
            is rejected and the electron is projected onto the Sphere."""

            """To determine path length within scattering medium
            first find the interesection with the surface of the scatterer."""
            scat_inter = self.scatteringMedium.getIntersection(old_position.copy(),
                                                             velocity)

            """Get the distance from intersection to previous position."""
            d = np.linalg.norm(old_position - scat_inter)

            """Append this distance to the pathlengths."""
            self.e.pathlength += d

            """Then find coordinates where electron will intersect the boundary."""
            intersect = self.boundary.getIntersection(position, velocity)
            """Set the electron's current position coordinates to the
            intersection coordinates."""
            self.e.vector[0:3] = intersect
            """Get new time."""
            distance = np.linalg.norm(intersect - position)
            new_time = distance / np.linalg.norm(velocity)
            self.e.vector[6] = new_time
            """Collect all intersected electrons into a list."""
            self.intersections += [list(self.e.vector)]
            self.intersected = True


    def Simulate(self):
        """Simulate the complete path of an electron from source to boundary.

        This function runs a complete set of simulation steps on an
        electron. It starts by generating the electron, then repeatedly
        performs simulation steps until the electron has either intersected
        the Boundary, or it has been scattered the set number of times
        (n_event_cutoff). The calcuation results are saved in the attribute
        self.results.
        """
        self.createElectron()
        self.intersected = False # this consdition is changed inside of the
        # function self._step()
        while ((not self.intersected) # while elec. has not reached boundary
            & (self.e.vector[7] < self.n_event_cutoff)
            & (self.e.kinetic_energy > 0)): # and not scattered too many times
            self._step()
            self.results = np.append(self.results,
                                     np.array([np.copy(self.e.vector[:])]),
                                     axis=0)

    def simulateMany(self, n: int, *args: str):
        """Run simulations for multiple electrons.

        Parameters:
        -----------
        n: int
            The number of simulations to run.
        *args:
            'keep all': If this argument is provided, all the calculated
                paths are kept in memory. results are stored in the
                attirbute self.all_paths.
            'start finish': If this argument is provided, then the
                starting properties, final properties and pathlengths of
                each electron are stored in the attribute self.start_finish
                which is
                Otherwise, only the final step in the calculation is saved,
                in the attribute self.intersections.
        """
        self.intersections  = []
        if 'keep all' in args:
            """ This keeps all the intermediate positions and velocities of
            every electron."""
            self.all_paths = []
            for step in range(n):
                self.Simulate()
                self.all_paths += [np.copy(self.results)]
        elif 'start finish' in args:
            start = []
            finish = []
            pathlengths = []
            for step in range(n):
                #print(step)
                self.Simulate()
                start += [np.copy(self.results[0])]
                finish += [np.copy(self.results[-1])]
                pathlengths += [self.e.pathlength]
            self.start_finish = [np.array(start), np.array(finish), np.array(pathlengths)]
        else:
            for step in range(n):
                self.Simulate() # this keeps the final positions, velocities
                # and scatter count in the attribute called self.results
        self.intersections = np.array(self.intersections)
        self.parameters['n'] = n

    def densityFromP(self, P: float, T: float = 300) -> float:
        """Calculate density, given pressure and temperature.

        Parameters
        ----------
        P : FLOAT
            Pressure in mbar.
        T : FLOAT
            Templerature in Kelvin.
        Returns
        -------
        density : FLOAT
            The particle density in particles per nm^3
        """
        R = 138 # in units of nm^3/mbar/K/atom
        density = P / (T * R)
        return density #returned in atoms / nm^3

#%% Set-up the simulation
if __name__ == '__main__':
    sim_radius = 5000
    source_radius = 6000
    sample_nozzle_distance = 5000
    boundary_shape = Sphere(sim_radius)
    initial_KE = 1400
    sim = Simulation(boundary_shape = boundary_shape)
    sim.addScatteringMedium(boundary_shape)
    # argument is the radius of the spherical simulation evnrionment
    sim.addSource(source_radius,1,initial_KE = initial_KE) # source is a shape that emits electrons
    # arguments are radius and thickness (in nm) of a disc
    sim.addScatterer(0.01,1.551,0.831,0.95,1.12,4)
    sim.scatterer.angle_dist = {'elastic':AngleDist(kind = 'Rutherford',
                                                    energy = initial_KE,
                                                    Z = sim.scatterer.Z,
                                                    param = 0.1),
                           'inelastic':AngleDist(kind = 'Constant')}

#%% Run the simulation
#Potentially long calculation.

    sim.simulateMany(5000, 'start finish')
    results = sim.start_finish


#%%
# plot all intersections with the boundary
    xyz = sim.intersections[:,0:3]
    fig = plt.figure(figsize = (6,6))
    ax = fig.add_subplot(111, projection='3d')
    xs = xyz[:,0]
    ys = xyz[:,1]
    zs = xyz[:,2]
    ax.scatter(xs, ys, zs)

    ax.set_xlabel('X Label')
    ax.set_ylabel('Y Label')
    ax.set_zlabel('Z Label')

    plt.show()


#%%
#This plots the angular distribution of the collected electrons
    data = sim.start_finish[1]
    counts = {}
    old_theta = 0
    for theta in np.arange(0, np.pi/2,np.pi/40):
        print(theta)
        counts[round(theta / np.pi * 180,2)]=np.array(
                [i for i in data if ((np.arcsin(i[5]
                / np.sqrt(i[3]**2 + i[4]**2 + i[5]**2))
                < theta)
                &
                (np.arcsin(i[5]
                / np.sqrt(i[3]**2 + i[4]**2 +i[5]**2))
                > old_theta))])
        old_theta = theta

    all = {}
    for i in counts:
        u = [j for j in counts[i]]
        all[i] = u

    all = np.array([[k,len(v)] for k,v in all.items()])



    unscattered = {}
    for i in counts:
        u = [j for j in counts[i] if ((j[7] == 0) & (j[8] == 0))]
        unscattered[i] = u

    un = np.array([[k,len(v)] for k,v in unscattered.items()])

    scattered_1 = {}
    for i in counts:
        u = [j for j in counts[i] if (j[7] == 1)]
        scattered_1[i] = u

    once = np.array([[k,len(v)] for k,v in scattered_1.items()])

    scattered_2 = {}
    for i in counts:
        u = [j for j in counts[i] if (j[7] == 2)]
        scattered_2[i] = u

    twice = np.array([[k,len(v)] for k,v in scattered_2.items()])

    scattered_3 = {}
    for i in counts:
        u = [j for j in counts[i] if (j[7] == 3)]
        scattered_3[i] = u

    thrice = np.array([[k,len(v)] for k,v in scattered_3.items()])

    #plt.plot(all[:,0],all[:,1])
    plt.plot(un[:,0],un[:,1])
    plt.plot(once[:,0],once[:,1])
    plt.plot(twice[:,0],twice[:,1])
    plt.plot(thrice[:,0],thrice[:,1])

    plt.show

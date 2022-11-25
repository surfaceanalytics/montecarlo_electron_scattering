# -*- coding: utf-8 -*-
"""
Created on Thu Apr 23 16:49:49 2020

@author: Mark
"""
from .simulation import Simulation
from ..model.generate_angle import AngleDist
from ..model.shapes import Disc

class GasPhaseSimulation(Simulation):
    """A class for simulating electron scattering through a gas phase."""

    def __init__(self, **kwargs):
        """Construct GasPhaseSimulation object.

        Parameters:
        ----------
            **kwargs:
                sim_radius: float
                    The radius of the Boundary environment in nm. In this case,
                    the Boundary shape is a cylinder.
                sample_nozzle_distance: float
                    The distance between the top of the Source object and the
                    top of the Boundary, in nm.
                pressure: float
                    The pressure of the gas phase scattering medium in mbar.
                initial_KE: float
                    The initial kinetic energy of the electrons in eV.
                source_diameter: float
                    The diameter of the Source Disc in nm.
                source_angle_dist: str
                    The angular distribution of the electrons generated from
                    the Source. One of ['Sphere', 'Lambert', 'Constant'].
        """
        super().__init__(**kwargs)

        '''This will be used as the radius of the simulation boundary'''
        if 'sim_radius' in kwargs.keys():
            sim_radius = kwargs['sim_radius']
        else:
            sim_radius = 1000000

        if 'sample_nozzle_distance' in kwargs.keys():
            sample_nozzle_distance = kwargs['sample_nozzle_distance']
        else:
            sample_nozzle_distance = 300000

        ''' This is the pressure of the gas phase in mbar'''
        if 'pressure' in kwargs.keys():
            pressure = kwargs['pressure']
        else:
            pressure = 10

        ''' This is the initial kinetic energy of the simulated electrons,
        in eV.
        '''
        if 'initial_KE' in kwargs.keys():
            initial_KE = kwargs['initial_KE']
        else:
            initial_KE = 1100

        if 'source_diameter' in kwargs.keys():
            source_diameter = kwargs['source_diameter']
        else:
            source_diameter = 300000

        if 'source_angle_dist' in kwargs.keys():
            source_angle_dist = kwargs['source_angle_dist']
        else:
            source_angle_dist = 'Lambert'

        """Here we define the boundary_shape."""
        boundary_shape = Disc(sim_radius,2*sample_nozzle_distance, center=[0,0,0])

        """Then we instantiate a simulation object and pass it the
        boundary shape.
        """
        super().__init__(boundary_shape = boundary_shape)
        """Set the pressure for the simulation. It returns the density
        in atoms/nm^3, given the pressure in mbar.
        """
        d = self.densityFromP(pressure)

        source_radius = source_diameter / 2

        """This adds a source, i.e. a shape that emits electrons. Here it is a
        Disc. Arguments are disc radius and disc thickness.
        """
        self.addSource(source_radius,1,initial_KE = initial_KE,
                       angle_distribution = source_angle_dist)

        """This defines the scattering medium. In the case of gas phase
        scattering it should be the same as the boundary.
        """
        self.addScatteringMedium(boundary_shape)

        """This sets the parameters of the scatterer.
        Arguments are:
            denstiy in atoms/nm^3.
            inel_factor, and inel_exp which determines the inelastic cross .
            section from a fit to TPP2M and NIST data using the equation:
            sigma = factor / (KE ^ exp).
            el_factor, which determines the elastic cross section
            and atomic number.
            One can also pass 'loss_function' as a keyword, where the value
            can either be the filename of a csv file, or a list of lists,
            containing the x-values (energy loss) and y-values (probability)
            of a loss function.
        """
        self.addScatterer(d,1.551,0.831,0.95,1.12,4, **kwargs)
        self.scatterer.angle_dist = {'elastic':AngleDist(kind = 'Rutherford',
                                                        energy = initial_KE,
                                                        Z = self.scatterer.Z,
                                                        param = 0.9),
                               'inelastic':AngleDist(kind = 'Constant')}
        '''
        self.addScatterer(d,0.0046,0,0.000001,0,4, **kwargs)
        self.scatterer.angle_dist = {'elastic':AngleDist(kind = 'Constant'),
                               'inelastic':AngleDist(kind = 'Constant')}'''

        self.height = 1

        self.parameters = {'sim_radius':sim_radius,
                           'sample_nozzle_distance': sample_nozzle_distance,
                           'pressure':pressure,
                           'initial_KE':initial_KE,
                           'source_diameter':source_diameter
                           }

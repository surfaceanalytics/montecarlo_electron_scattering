# -*- coding: utf-8 -*-
"""
Created on Thu Apr 23 16:49:49 2020

@author: Mark
"""
import numpy as np
import matplotlib.pyplot as plt
import xlsxwriter

#%%

class Analysis:
    def __init__(self, results, **kwargs):
        """The analysis class takes a complete results set.
        It only works if the 'start_finish' keyword was used in the 
        simulation.
        The 'results' atgument should be an n x 3 array, where n is the number
        of electrons simulated. The first element in the 'results' array is the 
        information of the starting position of the ith electron. It is a 
        1 x 7 array, with the first three elements [0:3] being position 
        corrdinates, the next three positions [3:6] being the velocity vector, 
        and the next two elements [6:8] are the number of times elastically
        and inelastically scattered, respectively.
        """
        self.results = results
        if 'parameters' in kwargs.keys():
            self.parameters = kwargs['parameters']
        else:
            self.parameters = {}
        
    def _unpack(self,results):
        return np.append(np.concatenate((results[0],results[1]), 
                                    axis=1).T, [results[2]], axis=0).T
    
    def _repack(self,unpacked_results):
    
        return [unpacked_results[:,:9],
                unpacked_results[:,9:18],
                unpacked_results[:,18]]
        
    def selectElectrons(self,*args, **kwargs):
        
        if 'select_all' in args:
            self.selection = self.results
        else:
            self.selection = self.selectHittingBoundary()
        
        if 'accept_angle' in kwargs.keys():
            """ Angle in degrees"""
            angle = kwargs['accept_angle']
            self.parameters['accept_angle'] = angle
            self.selection = self.selectAngle(angle)
                     
        if 'accept_radius' in kwargs.keys():
            radius = kwargs['accept_radius']
            self.selection = self.selectRadius(radius)
            
        if 'nozzle_diameter' in kwargs.keys():
            radius = kwargs['nozzle_diameter'] / 2
            self.parameters['nozzle_diameter'] = radius * 2
            self.selection = self.selectRadius(radius)          
            
    def selectHittingBoundary(self):
        """This determines all the electrons that have intersected the 
        boundary.
        We select only the electrons that hit the boundary by choosing
        all electrons that have been scatterred less than the cut-off number
        of times. If the electrons have been scattered the cut-off number
        of times, it means they did not make it to the boundary before the 
        calculation was terminated.
        """
        max_events = max(self.results[1][:,7])
        if max_events == 0:
            max_events = 1
        
        unpacked = self._unpack(self.results)
    
        unpacked = np.take(unpacked, np.where(self.results[1][:,7]
                                              <max_events)[0].tolist(), axis=0)
        
        repacked = self._repack(unpacked)
        
        return repacked
    
    def selectAngle(self, angle):
        """This method selects only the electrons whose angle between the
        velocity vector and the z-axis is less than the acceptance angle.
        """
        def cos(line):
                theta = np.arccos(line[2] / (np.sqrt(line[0]**2 +
                                                     line[1]**2 + 
                                                     line[2]**2))) * 180 / np.pi
                return theta
        
        unpacked = self._unpack(self.selection)
        
        filtered = np.array([i for i in unpacked if (np.abs(cos(i[12:15])) < angle)])
        
        repacked = self._repack(filtered)
        return repacked
    
    def selectRadius(self, radius):
        """ Select only the selectrons whose intersection with the boundary
        is within the radius of the central z-axis."""
        
        unpacked = self._unpack(self.selection)
        filtered = np.array([i for i in unpacked
                         if (((i[9]**2 + i[10]**2) < (radius)**2)
                             & (i[11] > 0))])
        repacked = self._repack(filtered)
        return repacked


    def pathHistogram(self,*args, bins=100, **kwargs):
        """Gets the number of times an electron was inelastically scattered
        along its path. 
        For electrons that do not intersect the boundary, this should 
        represent the cutoff for number of times scattering in the simulation 
        before the path calculation was terminated.
        
        Use the argument 'show' to show the plot after the calculation is 
        complete, and the keyword argument 'x_limits' with a tuple 
        (x_min, x_max) to adjust the axis limits.
        
        To limit the electrons to only the ones that have interesected the
        boundary within a desired acceptance angle, provide the keyword arg
        'accept_angle' = angle, where angle should be in degrees.
        """

        self.selectElectrons(**kwargs)

        """Get a list of the number times an electron was inelastically 
        scattered for all electrons hitting the sphere.
        """
   
        inel_count = self.selection[1][:,7]
        
        pathlengths = self.selection[2][:]
            
        profiles = {}
        for i in range(int(max(inel_count))+1):
            hist = np.histogram(pathlengths[np.where(inel_count == i)], 
                                            bins=bins)
            profiles[i] = {'counts':hist[0], 'pathlength':hist[1][:-1]}
        
        if 'show' in args:
            fig = plt.figure(figsize=(5,5))
            ax = fig.gca()
            for k,v in profiles.items():
                ax.plot(v['pathlength'], v['counts'])
            if 'x_limits' in kwargs.keys():
                limits = kwargs['x_limits']
                ax.set_xlim(limits)
            ax.set_xlabel('Pathlength [nm]', linespacing=3)
            ax.set_ylabel('Electron count', linespacing=3)
            plt.xticks(rotation=90)
            plt.show()
        
        self.profiles = profiles
    
    def areasUnderProfiles(self, *args):
        """ This plots the total areas under the path histogram profiles.
        """
        
        areas = {}
        for k,v in self.profiles.items():
            areas[k] = np.sum(v['counts'])
        self.areas = areas
        if 'show' in args:
            fig = plt.figure(figsize=(5,5))
            ax = fig.gca()
            ax.scatter(list(areas.keys()),list(areas.values()))
            ax.set_xlabel('Number of times inelastically scattered', linespacing=3)
            ax.set_ylabel('Electron count', linespacing=3)
            plt.show()
            
    def plotStartFinish(self, x_lim=800000, n_electrons=500, **kwargs):
        hitting_boundary = self.selection[1]

        acceptance = hitting_boundary
        
        if n_electrons > len(acceptance):
            n_electrons = len(acceptance)
        
        fig = plt.figure(figsize=(10,10))
        xlim = x_lim
        ylim = xlim
        zlim = xlim
        ax = fig.gca(projection='3d')
        ax.set_xlabel('x distance', linespacing=3)
        ax.set_ylabel('y distance', linespacing=3)
        ax.set_zlabel('z distance', linespacing=3)
        for l, k in enumerate(acceptance[:n_electrons]):
            V = k
            x, y, z = V[0], V[1], V[2]
    
            ax.scatter3D(x, y, z, c='blue')
        
        for l,k in enumerate(self.results[0][:n_electrons]):
            V = k
            x, y, z = V[0], V[1], V[2]
    
            ax.scatter3D(x, y, z, c='orange')
            
        ax.set_xlim3d(-xlim, xlim)
        ax.set_ylim3d(-ylim,ylim)
        ax.set_zlim3d(-100000,zlim)
        
        ax.set_xticks([-xlim, 0, xlim])
        ax.set_yticks([-ylim, 0, ylim])
        ax.set_zticks([-zlim, 0, zlim])
        
        ax.view_init(10, 90)    
        plt.show()
        
    def writeExcel(self, file):
        file = file + '.xlsx'
         
        with xlsxwriter.Workbook(file) as workbook:
            
            if len(self.profiles) != 0:
                self.addProfilesToExcel(workbook)
           
            if len(self.areas) != 0:
                self.addAreasToExcel(workbook)
                 
            if hasattr(self, 'spectrum'):
                self.addSpectrumToExcel(workbook)
                
            if hasattr(self, 'partial_spectra'):
                self.addPartialSpectraToExcel(workbook)
                
            #workbook.close()
            

    def getColLabel(self, idx):
        column_mapping = ['A','B','C','D','E','F','G','H','I','J','K','L','M',
              'N','O','P','Q','R','S','T','U','V','W','X','Y','Z']
        if (idx+1) > len(column_mapping):
            label = column_mapping[int((idx+1) / len(column_mapping))-1]
            label += column_mapping[(idx+1) % len(column_mapping) -1]
            
        else:
            label = column_mapping[idx]
            
        return label
            
    def addProfilesToExcel(self, workbook):
        sheet_name = 'profiles'
        worksheet = workbook.add_worksheet(sheet_name)
        start_col = 0
        start_row = 3
        cols_per_profile = 2
        x_axis = 'path length [nm]'
        y_axis = 'counts'
        txt = ''
        for k,v in self.parameters.items():
            txt += str(k) + " : " + str(v) + ",\n"
            options = {'width':300,'height':200} 
        worksheet.insert_textbox(0,0,txt,options)
        
        chart = workbook.add_chart({'type':'scatter', 'subtype':'straight'})
        
        for k, v in self.profiles.items():
            worksheet.write(0,start_col, 'n times scattered :')
            worksheet.write(0,start_col +1, k)
            worksheet.write(1,start_col, x_axis)
            worksheet.write(1,start_col+1, y_axis)
            x = v['pathlength']
            y = v['counts']
            for i, val in enumerate(x):
                worksheet.write(start_row+i,start_col, val)
            for i, val in enumerate(y):
                worksheet.write(start_row+i, start_col+1, val)
            
            chart.add_series({
                'name':'n = '+str(k),
                'categories':'='+sheet_name+'!$'+self.getColLabel(start_col)+'$'+str(start_row)+':$'+self.getColLabel(start_col)+'$'+str(i+start_row),
                'values':'='+sheet_name+'!$'+self.getColLabel(start_col+1)+'$'+str(start_row)+':$'+self.getColLabel(start_col+1)+'$'+str(i+start_row)
                })
 
            start_col += cols_per_profile
            
        chart.set_title({'name':'Profiles'})
        chart.set_x_axis({'name':x_axis})
        chart.set_y_axis({'name':y_axis})
        chart.set_style(2)
        worksheet.insert_chart('D2',chart)
                  
    def addAreasToExcel(self, workbook): 
        worksheet = workbook.add_worksheet('areas')
        start_col = 0
        cols_per_profile = 2
        worksheet.write(0,start_col, 'n times scattered')
        worksheet.write(0,start_col +1, ' counts')
        i=1
        for k, v in self.areas.items():
            worksheet.write(i,start_col, k)
            worksheet.write(i,start_col+1, v)
            i+=1
           
    def addSpectrumToExcel(self, workbook):
        sheet_name = 'spectrum'
        worksheet = workbook.add_worksheet(sheet_name)
        x_axis = 'Kinetic energy [eV]'
        y_axis = 'Electron count'
        start_col = 0
        start_row = 2
        worksheet.write(0,start_col, x_axis)
        worksheet.write(0, start_col+1, y_axis)
        for i in range(len(self.spectrum[0])):
            kinetic_energies = self.spectrum[1]
            electron_counts = self.spectrum[0]
            worksheet.write(i+start_row,start_col, kinetic_energies[i])
            worksheet.write(i+start_row,start_col+1, electron_counts[i])
            
        chart = workbook.add_chart({'type':'scatter',
                                    'subtype':'straight'})
        chart.add_series({
            'name':'spectrum',
            'categories':'='+sheet_name+'!$'+self.getColLabel(start_col)+'$'+str(start_row)+':$'+self.getColLabel(start_col)+'$'+str(i+start_row),
            'values':'='+sheet_name+'!$'+self.getColLabel(start_col+1)+'$'+str(start_row)+':$'+self.getColLabel(start_col+1)+'$'+str(i+start_row)
            })
        chart.set_title({'name':'Simulated spectrum'})
        chart.set_x_axis({'name':x_axis})
        chart.set_y_axis({'name':y_axis})
        chart.set_style(2)
        worksheet.insert_chart('D2',chart)
        
        
    def addPartialSpectraToExcel(self, workbook):
        sheet_name = 'partial_spectra'
        worksheet = workbook.add_worksheet(sheet_name)
        chart = workbook.add_chart({'type':'scatter',
                                    'subtype':'straight'})
        
        x_axis = 'Kinetic energy [eV]'
        y_axis = 'Electron count'
        start_col = 0
        start_row = 2

        col1 = start_col
        worksheet.write(0,col1, x_axis)
        for idx, i in enumerate(self.partial_spectra[0][1]):
                worksheet.write(idx+start_row,col1, i)
        
        
        for idx, item in enumerate(self.partial_spectra.items()):
            spectrum = item[1]
            n = item[0]
            if idx < 30:
                col2 = col1+1+idx
                worksheet.write(0,col2, y_axis)
                worksheet.write(1,col2, 'n = '+str(n))
                
                for i in range(len(spectrum[0])):
                    
                    electron_counts = spectrum[0]
                    worksheet.write(i+start_row,col2, electron_counts[i])
                    chart.add_series({
                        'name':'n = '+ str(n),
                        'categories':'='+sheet_name+'!$'+self.getColLabel(col1)+'$'+str(start_row)+':$'+self.getColLabel(col1)+'$'+str(i+start_row),
                        'values':'='+sheet_name+'!$'+self.getColLabel(col2)+'$'+str(start_row)+':$'+self.getColLabel(col2)+'$'+str(i+start_row)
                        })
                    
        chart.set_title({'name':'Partial spectra'})
        chart.set_x_axis({'name':x_axis})
        chart.set_y_axis({'name':y_axis})
        chart.set_style(2)
        worksheet.insert_chart('D2',chart)
  
    def getAveragePathLength(self):
        total_length = 0
        total_counts = 0
        if self.profiles != 0:

            for k,v in self.profiles.items():
                counts = v['counts']
                lengths = v['pathlength']
                total_length += np.sum(np.multiply(counts, lengths))
                total_counts += np.sum(counts)
        if total_counts != 0:
            avg_path = total_length / total_counts
            return avg_path
        
    def showSpectrum(self, **kwargs):
        def kineticEnergy(vector):
            '''This function returns kinetic energy in eV
            '''
            mass = 9.109383E-31 # the electron mass in kg
            J_eV = 1.602176E-19
            speed = np.linalg.norm(vector[3:6]) / 1E+9 # here the speed needs
            # to be converted from nm/s to m/s
            KE = (1/2 * mass * speed**2) / J_eV
            return KE
        
        if 'collected_only' in kwargs.keys():
            if kwargs['collected_only'] == True:
                results = self.selection[1]
            else: 
                results = self.results[1]
        else: results = self.results[1]
        

        KEs = [kineticEnergy(v) for v in results]
        
        
        if 'energy_range' in kwargs.keys():
            _range = kwargs['energy_range']
        else: 
            _range = (0, max(KEs))
            
        if 'bins' in kwargs.keys():
            bins = kwargs['bins']
        elif 'step_size' in kwargs.keys():
            bins = int((_range[1] - _range[0]) / kwargs['step_size'])
        else:
            bins = 100   
            
        temp_hist = np.histogram(KEs, bins=bins, range=(_range))
        hist = (temp_hist[0], temp_hist[1][:-1])

        self.spectrum = hist
        #plt.plot(hist[1],hist[0])
        #plt.show()
        
        fig = plt.figure(figsize=(5,5))
        ax = fig.gca()
        ax.plot(hist[1],hist[0])
        ax.set_xlabel('Kinetic energy [eV]', linespacing=3)
        ax.set_ylabel('Electron count', linespacing=3)
        plt.xticks(rotation=90)
        plt.show()
        
    def showPartialSpectra(self, **kwargs):
        def kineticEnergy(vector):
            '''This function returns kinetic energy in eV
            '''
            mass = 9.109383E-31 # the electron mass in kg
            J_eV = 1.602176E-19
            speed = np.linalg.norm(vector[3:6]) / 1E+9 # here the speed needs
            # to be converted from nm/s to m/s
            KE = (1/2 * mass * speed**2) / J_eV
            return KE
        
        if 'collected_only' in kwargs.keys():
            if kwargs['collected_only'] == True:
                results = self.selection[1]
            else: 
                results = self.results[1]
        else: results = self.results[1]

        KEs = [kineticEnergy(v) for v in results]
        if 'energy_range' in kwargs.keys():
            _range = kwargs['energy_range']
        else: 
            _range = (0, max(KEs))
            
        if 'bins' in kwargs.keys():
            bins = kwargs['bins']
        elif 'step_size' in kwargs.keys():
            bins = int((_range[1] - _range[0]) / kwargs['step_size'])
        else:
            bins = 100   
            
        fig = plt.figure(figsize=(5,5))
        ax = fig.gca()    
        all_n = sorted(set([v[7] for v in results]))
        partial_spectra = {}
        for n in all_n:
            KEs = [kineticEnergy(v) for v in results if v[7] == n]
            temp_hist = np.histogram(KEs, bins=bins, range=(_range))
            hist = (temp_hist[0].copy(), temp_hist[1][:-1].copy())
            partial_spectra[n] = hist  
            ax.plot(partial_spectra[n][1],partial_spectra[n][0])

        self.partial_spectra = partial_spectra
        #plt.plot(hist[1],hist[0])
        #plt.show()

        ax.set_xlabel('Kinetic energy [eV]', linespacing=3)
        ax.set_ylabel('Electron count', linespacing=3)
        plt.xticks(rotation=90)
        plt.show()
      

            



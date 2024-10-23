'''
 Desc  : Reader and Writer for RENSHAPE
 Author: Matteo Borghesi <matteo.borghesi@mib.infn.it>
'''


import numpy as np
import h5py
from base import base_utilities
from scipy import integrate


class LazyWriter:
    
    '''
    Simple writer for the .lazy file (hdf5).
    
    '''
    
    def __init__(self, output_path=None,name=None,fname=None):
        

        if output_path is not None:
            output_path = base_utilities.fix_path(output_path)
            self._file = output_path+name  +'.lazy'
        if fname is not None:
            self._file = fname

        return
        
    def get_file_name(self):
        return self._file
    
    def __check_exist(self,path=None,group_name = None):

        with h5py.File(self._file,'r') as f:
            if path is None:
                names = list(f.keys())   
            else:
                names = list(f[path].keys() ) 
            if group_name in names:
                return True
            else:
                return False
    
    def set_general_info(self, dictionary=None):
        
        '''
        Write "dictionary" in the .lazy header
        
        '''
                
        with h5py.File(self._file,'a') as f:

            if self.__check_exist(group_name = "info") is False:
                general_info=f.create_group('info')
            else:
                general_info=f.get("info")
        
            for key, value in dictionary.items():
                self.__save_parameter(general_info,key,value)
        return
    
    
    def write_nuclide_data(self, nuclide_name=None, dtype="info",dictionary=None, vname = None, vvalue = None):        


        big_group = "nuclides"       
        
        
        with h5py.File(self._file,'a') as f:
        
            if self.__check_exist(path=None,group_name = big_group) is False:
                f.create_group(big_group)
            
            if self.__check_exist(path=big_group,group_name = nuclide_name) is False:
                group=f.create_group(big_group + "/"+ nuclide_name)
            else:
                group=f.get(big_group + "/"+ nuclide_name)        
                    
            nuclide_name = big_group + "/"+ nuclide_name
            
            if self.__check_exist(path=nuclide_name,group_name = dtype) is False:
                group=f.create_group(nuclide_name+"/"+dtype)
            else:
                group=f.get(nuclide_name+"/"+dtype)      
                
            if dictionary is not None:    
                for key, value in dictionary.items():
                    self.__save_parameter(group,key,value)            
            else:
                self.__save_parameter(group,vname,vvalue)
        return
        

    def __save_parameter(self,group,variable_name,value):
        try:
            group.create_dataset(variable_name, data=value)        
        except (ValueError):     #if exist, overwrite it
            del group[variable_name]
            group.create_dataset(variable_name, data=value) 
            #group[variable_name][...] = value 
        return
 
 
        with h5py.File(self.__file,'a') as f:       
            f.create_dataset(path+'/'+name,data=value)
        return

    def __delete_parameter(self,path,name):
        with h5py.File(self.__file,  "a") as f:
            del f[path+'/'+name]
        return



class LazyReader:

    '''
    Simple reader for the .lazy file (hdf5).
    
    '''
        
    def __init__(self, path_file):
        self._file = path_file
        
    def __check_exist(self,path=None,group_name = None):

        with h5py.File(self._file,'r') as f:
            if path is None:
                names = list(f.keys())   
            else:
                names = list(f[path].keys() ) 
            if group_name in names:
                return True
            else:
                return False

    def __convert_to_dict(self, group):
        dic = {}

        for key in group.keys():
            el = group[key]
            if el.shape != ():
                el = np.array(el)
            elif el.dtype == 'object':
                el = str(np.array(el).astype(str))
            else:
                el = float(np.array(el).astype(float))
            dic[key] = el
        return dic    

    
    def get_info(self):
        """
        Get the header of the .lazy file.
        """
        with h5py.File(self._file,'r') as f:
            return self.__convert_to_dict(f['info'])
    
    def get_nuclides_list(self, group_path = "nuclides"):
        with h5py.File(self._file,'r') as f:
            return list(f[group_path].keys())
            
    def get_parameters_labels(self,attempts=20):
        with h5py.File(self._file,'r') as f:
            names = self.get_nuclides_list()[:attempts]    
            temp = []
            for name in names:
                temp.append(list(f["nuclides/"+name+"/info"].keys()))
            lenght = np.array([len(a) for a in temp]) 
        return temp[np.argmax(lenght)]     

    def get_parameters(self,name=None, group_path="nuclides",sub_group = "info",variable_not_found=-2):
        """
        Return the values for the selected parameter for all the nuclides in the dataset.
        
        Parameters
        ----------
        name : string
            Name of the selected parameter. Default is None.
        group_path : string
            Group name where the parameter is saved. Default is "nuclides".
        sub_group : string
            Sub-group where the parameter is saved. Default is "info".
        variable_not_found : float
            If the parameter is not present for a given nuclide, the value for that nuclide will be set
            to this number. Default is -2.
            
        Returns
        -------
        temp : ndarray
            Array with the selected parameters for all nuclides. The i-th entry of the array
            corresponds to the i-th nuclide.
            
        Examples
        --------
        >>> reader.get_parameters(name="cumulative_thermal_fy_235u", variable_not_found = 0)
        
        """
        with h5py.File(self._file,'r') as f:
            temp = []
            for el in f[group_path].keys():
                try:
                    value = f[group_path][el][sub_group][name]
                except:
                    value = variable_not_found
                temp.append(np.array(value))
            temp = np.array(temp)
            if temp.dtype == "O":
                temp = temp.astype(str)
        return temp
        
    def get_data(self,name="dN_dE_tot", group_path="nuclides",sub_group = "data",E_min = 0, E_max = 12000,E_step = None):
        """
        Return the arrays labeled as *name* present in the dataset. The arrays are zero-padded and re-arranged 
        as matrix.
        
        Parameters
        ----------
        name : string
            name of the variable to collect for each nuclide in the dataset. Default is 
            "dN_dE_tot".
        group_path : string
            Group name where the parameter is saved. Default is "nuclides".
        sub_group : string
            Sub-group where the parameter is saved. Default is "data".
        E_min : float
            The lowest energy for each spectrum in keV. Default is 0.
        E_max : float
            The highest energy for each spectrum in keV. Default is 12000.   
        E_step : float
            The energy step for each spectrum. If None, the value is taken from the header
            of the dataset. Default is None.
            
        Returns
        -------
        data : ndarray
            Data matrix. Each row corresponds to a nuclide. 
        energies : ndarray
            Array with the energy in keV. len(energies) = data.shape[1].
        pos_ok : ndarray
            The indices to match the rows of data with the output of get_parameters and get_nuclide.
            e.g. :  data[i] = get_nuclide(loc=pos_ok[i])["dN_dE_tot"]
                    data[i] -> get_parameters("Q")[pos_ok[i]]
        pos_notok : ndarray
            The indices of the nuclides which are not present in data. 
        
        """
        if E_step is None:
            E_step = self.get_info()["E_step"]
            
        energies = np.arange(0,E_max,step=E_step) #energy in the dataset
        pos_min, pos_max = np.where((energies>=E_min)&(energies<=E_max))[0][[0,-1]]
        data_lenght = int(pos_max-pos_min)
        with h5py.File(self._file,'r') as f:
            temp = []
            pos_ok = []
            pos_notok = []
            i = 0
            for el in f[group_path].keys():
                try:
                    value = np.array(f[group_path][el][sub_group][name])
                    pos_ok.append(i)
                except:
                    pos_notok.append(i)
                    i+=1
                    continue
                value = value[pos_min:pos_max+1]
                value = np.append(value,np.zeros(data_lenght-len(value)+1))  
                temp.append(value)
                i+=1
        energies = energies[pos_min:pos_max+1]
        data = np.array(temp)
        pos_ok = np.array(pos_ok)
        pos_notok = np.array(pos_notok)
        return data,energies,pos_ok,pos_notok
    
    def get_nuclide(self,name=None,loc=None,group_path="nuclides"):
        subgroup1_name = "info"
        subgroup2_name = "data"
        with h5py.File(self._file,'r') as f:
            if name is not None:
                nuclide_name = name
            if loc is not None:
                nuclide_name = self.get_nuclides_list()[loc]
            
            dic1 = self.__convert_to_dict(f[group_path][nuclide_name][subgroup1_name])
            if self.__check_exist(path=group_path+"/"+nuclide_name,group_name = subgroup2_name) is True:
                dic2 = self.__convert_to_dict(f[group_path][nuclide_name][subgroup2_name])
                dic1.update(dic2)
            dic1['lazy_name'] = nuclide_name
        return dic1
        
    def _evaluate_total_spectrum(self,loc,thr=0.2):
        """
        Return the neutrino spectrum from the neuclide in loc position.
        The uncertainty is given considering also the uncertainty on the BRs.
        """
        nuc = self.get_nuclide(loc=loc)
   
        # Normalize the transition spectrum to 1 insted of it's BR.
        y =(nuc["transition_dN_dE"].T/nuc["transition_intensity"]).T
        y_er = (nuc["transition_unc_dN_dE"].T/nuc["transition_intensity"]).T
    
        unc_BR = []
        for i in range(len(nuc["transition_intensity"])):
            unc_BR.append(base_utilities.short_to_std(nuc["transition_intensity"][i],nuc["transition_unc_intensity"][i]))
        unc_BR = np.array(unc_BR)

        #substitute the BR with unknown uncertainty with an uncertainty equal to thr*BR
        ii = np.where(unc_BR==None)[0]
        unc_BR[ii] = nuc["transition_intensity"][ii]*thr
    
        y_er2 = np.sum( ((y.T*unc_BR).T)**2 + ((y_er.T*nuc["transition_intensity"]).T)**2, axis=0)**0.5
    
        return np.sum((nuc["transition_dN_dE"]),axis=0),y_er2
        
        
    def get_nu_spectra(self,E_min = 0, E_max=12e3, unc_BR = True, default_unc = 0.2, 
                       force_normalization=True,thr_norm=0.99):
        '''
        Return the neutrino spectra of the nuclides in the dataset with their uncertainties.
        
        Parameters
        ----------
        E_min : float
            The lowest energy for each spectrum in keV. Default is 0.
        E_max : float
            The highest energy for each spectrum in keV. Default is 12000.   
        unc_BR : bool
            Evaluate the uncertanties for each spectrum considering also the branching 
            ratio uncertanties. Default is True.
        default_unc : float
            If a nuclide does not have uncertanties or one of its levels do not have a known error, default_unc
            is used as percentage error instead. Default is 0.2.
        force_normalization : bool
            If True, each spectrum with an integral value greater than *thr_norm* will be normalized to one.
            Default is True.
        thr_norm : float
            Only if force_normalization is True. Default is 0.99.
            
        Returns
        -------
            energy : ndarray
                Array with the energy in keV. len(energies) = spectra.shape[1].
            spectra : ndarray
                Data matrix. Each row corresponds to spectrum for a nuclide. 
            spectra_er : ndarray
                Data matrix. Each row corresponds to the uncertanties related to the spectrum for a nuclide. 
            posOK : ndarray
                The indices to match the rows of *spectra* with the output of get_parameters and get_nuclide.
                e.g. :  spectra[i] = get_nuclide(loc=pos_ok[i])["dN_dE_tot"]
                        spectra[i] -> get_parameters("Q")[pos_ok[i]]
            posnotOK : ndarray    
                The indices of the nuclides which are not present in *spectra*. 
        
        '''
    
        spectra, energy, posOK, posnotOK = self.get_data("dN_dE_tot",E_min = E_min, E_max=E_max)
    
        if unc_BR is False:
            spectra_er, _, _, _ = self.get_data("unc_dN_dE",E_min = E_min, E_max=E_max)
        
            #if the uncertainty is unkwnown, use the defalt_unc
            pos_to_fix = np.where(np.mean(spectra_er,axis=1)==0)[0]
            for p in pos_to_fix:
                spectra_er[p] = spectra[p]*default_unc
        else:
            spectra_er = np.zeros((len(posOK),int(E_max)))
            i = 0
            for p in posOK:
                try:
                    _, temp = self._evaluate_total_spectrum(self,p,thr=default_unc)
                    temp = temp[:spectra.shape[1]]  
                except:
                    #if the uncertainty is unkwnown, use the defalt_unc
                    temp = spectra[i]*default_unc
                spectra_er[i,:len(temp)] = temp
                i += 1
          
        if force_normalization is True:
            for i in range(spectra.shape[0]):
                I = integrate.simpson(spectra[i],x=energy)
                if I >= thr_norm:
                    spectra[i] = spectra[i]/I
                    spectra_er[i] = spectra_er[i]/I        
    
        return energy, spectra, spectra_er, posOK, posnotOK
     

    def get_total_spectrum(self,labels=["cumulative_thermal_fy_235u","cumulative_thermal_fy_239Pu","cumulative_thermal_fy_241Pu","cumulative_fast_fy_238u"],
                           labels_unc=["unc_ct_235u","unc_ct_239Pu","unc_ct_241Pu","unc_cf_238u"],
                           ffs = [0.564,0.076,0.304,0.056], ffs_unc = None, do_sum = True,
                           spectra=None,spectra_er=None,posOK=None):

        '''
        Return the neutrino spectrum with its uncertainty given the fission fractions.
        
        Parameters
        ----------
        labels : list of string
            Labels for the cumulative fission fractions to use. Default is ["cumulative_thermal_fy_235u","cumulative_thermal_fy_239Pu","cumulative_thermal_fy_241Pu","cumulative_fast_fy_238u"].
        labels_unc : list of string
            Labels for the uncertainties associated with the fission fractions. There must be consistency between labels, labels_unc and ffs.
            Default is ["unc_ct_235u","unc_ct_239Pu","unc_ct_241Pu","unc_cf_238u"].
        ffs : list of float
            Fission fractions to use. There must be consistency between labels, labels_unc and ffs. Default is  [0.564,0.076,0.304,0.056].
        ffs_unc : list of float or None
            Uncertainties for the fission fraction. If None, no uncertainty is used. Default is None.
        do_sum : bool
            If True, sum all the spectra to return the total antineutrino spectrum. Otherwise return a matrix. Default is True.
        spectra : ndarray
            The matrix containing all the spectra to be combined. It should be the output of "get_nu_spectra". Default is None.
        spectra_er : ndarray
            The matrix containing all the uncertainties for the spectra. It should be the output of "get_nu_spectra". Default is None.
        posOK : ndarray
            The indices to match the rows of *spectra* with the output of get_parameters and get_nuclide. It should be the output of "get_nu_spectra". Default is None.
            
        Returns
        -------
            spectrum : ndarray
                The reactor antineutrino spectrum
            spectrum_err : ndarray
                The uncertainties in the reactor antineutrino spectrum
        '''
        
    

        sum_cfy = np.zeros(len(posOK))
        sum_cfy_unc = np.zeros(len(posOK))
    
        for i in range(len(labels)):
            cfy_selected = self.get_parameters(labels[i] ,variable_not_found=0)
            sum_cfy += cfy_selected[posOK]*ffs[i]
        
            sum_cfy_unc += (self.get_parameters(labels_unc[i],variable_not_found=0)[posOK]*ffs[i])**2
            if ffs_unc is not None:
                sum_cfy_unc += (cfy_selected*ffs_unc[i])**2
        
        sum_cfy_unc = sum_cfy_unc**0.5

        spectra_cfy = (spectra.T * sum_cfy).T
    
        spectra_cfy_err = ((spectra_er.T * sum_cfy).T)**2
        spectra_cfy_err += ((spectra.T*sum_cfy_unc).T)**2

        if do_sum is True:
            spectrum = np.sum(spectra_cfy,axis=0)
            spectrum_err = np.sum(spectra_cfy_err,axis=0)**0.5   
            return spectrum,  spectrum_err
        else:
            return spectra_cfy, spectra_cfy_err

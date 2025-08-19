'''
 Author: Matteo Borghesi <matteo.borghesi@mib.infn.it>
 Desc  : Collection of utilities used (?) in other codes 
'''

import os, errno
from operator import methodcaller
import numpy as np
import string
import uncertainties
import re
import periodictable
from base import base_utilities as bu

SYMBOL_LIST = []
for el in periodictable.elements: 
    SYMBOL_LIST.append(el.symbol)

def lazy_to_ensdf(list_of_names=None,dtype=".ensdf"):
    
    #remove the unwanted _m from the nuclides names
    pos_cor = locate_word(list_of_names,'_')  
    for i in pos_cor:
        list_of_names[i] = list_of_names[i][:-3]
        
    #make all uppercase and add the correct dtype
    for i in range(len(list_of_names)):
        list_of_names[i] = list_of_names[i].upper() + dtype

    list_of_names = np.unique(np.array(list_of_names))
    return list_of_names    
    
def get_ensdf_file(lname=None,ensdf_path=None,dtype = ".ensdf"):
    '''
    Return all files with the .ensdf format in ensdf_path. 
    If lname is not None, return the files which are also present in the lname lazy file.
    '''
    if lname is None:
        lista = np.array(listdir(ensdf_path))
        temp = []
        for i in range(len(lista)):
            if lista[i].find(dtype) > 0:
                temp.append(lista[i])
        return np.array(temp)
    else:
        reader = lazy_handler.LazyReader(lname)
        nuclides_for_betashape = bu.lazy_to_ensdf(reader.get_nuclides_list())
        path_files = np.array(listdir(ensdf_path))
        
        common_files = np.intersect1d(nuclides_for_betashape,path_files)
        not_present =  np.array(list(set(nuclides_for_betashape) - set(common_files)))
        return common_files, not_present
        
        
def create_dict(z = None,a=None,m=0,nfile=None, 
                data = None,symbol_list=SYMBOL_LIST,example_dic=None):
    '''
    Convert the data from ENDF/B sub library into ENSDF-like data for BetaShape.
    
    Parameters
    ----------
    z : int
        Atomic number. Default is None.
    a : int
        Atomi mass. Default is None.
    m : int
        Metastable index. m=0 no metastable, m=1 first metastable... Default is None.
    nfile : string
        Optional. Path for the endf/b sub-library. Using (z,a,m), it will extract the corresponding data matrix.
        In that case, data variable must be None. Default None.
    data : np.array
        Optional. Required if nfile is None. Data matrix (n,5), with n is the number of levels 
        in the decay scheme. Each row must have the following syntax: [end point energy (MeV), uncertainty (MeV),
        intensity, uncertainty , type of transition].
        intensity <=1.
        uncertainty is expressed as 1 standard deviation.
        transition: 'a' for allowed, '1u' for first-forbidden unique, '2u' for second-forbidden unique...
    example_dic : dictionary
        Is the example dictionary for BetaShape. See CmdBetaShape class. It's the output of the method
        get_example_dictionary()
            
    Returns
    -------
    dic : dictionary
        Dictionary containing all the previous information written in an ENSDF-like data format.
        When used as input for the evaluate_decay method of CmdBetaShape, it returns the
        beta and neutrino spectra with the characteristics specified by the input data matrix.
    '''
    
    if data is None:
        endf_sub_reader = endf_reader.EndfBSubLibraryReader(nfile)
        _,_, data = endf_sub_reader.get_element(z=z,a=a,m=m)
    
    if len(data.shape) == 1:
        data = np.array([data])
    E0 = (data[:,0].astype(float)*1e3).flatten()  #keV
    d_E0 = (data[:,1].astype(float)*1e3).flatten()  #keV

    I = (data[:,2].astype(float)*100).flatten()
    d_I = (data[:,3].astype(float)*100).flatten()
    trans = (data[:,4]).flatten()
    
    dic = example_dic
    #=========
    dic["P"] = str(a)+symbol_list[z].upper()
    dic["D"] = str(a)+symbol_list[z+1].upper()
    dic["JP_P"] = "0+"

    dic["Q"] = np.max(E0)
    dic["dQ"] = d_E0[np.argmax(E0)]
        
    res = bu.std_to_shorthand(dic["Q"],dic["dQ"])

    dic["dQ"] = res[1] #sistemato per ensdf
    dic["Q"] = res[0]        

    dic["IB"] = np.round(I,4)
    dic["d_IB"] = np.round(d_I,4)
  
    temp1 = []
    for i in range(len(dic["IB"])):
        if dic["d_IB"][i] == 0:
            temp1.append((dic["IB"][i],0))
        else:
            temp1.append(bu.std_to_shorthand(dic["IB"][i],dic["d_IB"][i]))
    temp1 = np.array(temp1)

    dic["IB"] = temp1[:,0]
    dic["d_IB"] = temp1[:,1].astype(int)

    dic["E_D"] = float(dic["Q"])- E0
    dic["dE_D"] = d_E0
    temp1 = []
    for i in range(len(dic["E_D"])):
        temp1.append(bu.std_to_shorthand(dic["E_D"][i],dic["dE_D"][i]))
    temp1 = np.array(temp1)
    dic["E_D"] = temp1[:,0]
    dic["dE_D"] = temp1[:,1]


    temp = []
    for el in trans:
        if el == "a":
            temp.append("0+")
        if el == "1u":
            temp.append("2-")
        if el == "2u":
            temp.append("3+")
        if el == "3u":
            temp.append("4-")
    dic["JP_D"] = temp
    return dic



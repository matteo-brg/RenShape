#!/usr/bin/env python
import numpy as np
import argparse
from os import listdir
import periodictable
import os


from base import base_utilities as bu
from rw import endf_reader
from rw import lazy_handler
from process import wrappers

def lazy_to_ensdf(list_of_names=None,dtype=".ensdf"):

    #remove the unwanted _m from the nuclides names
    pos_cor = bu.locate_word(list_of_names,'_')  
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
                data = None,reader = None,symbol_list=None,example_dic=None):
    '''
    Convert the data from ENDF/B sub library into ENSDF-like data for BetaShape.
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

    if reader is not None:
        dic["Q"] = reader.get_nuclide(name=str(a)+symbol_list[z])["Q"]
        dic["dQ"] = reader.get_nuclide(name=str(a)+symbol_list[z])["unc_Q"]
    else:
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






        

def main():
    
    usage='createLazyFile.py -lp /path/to/.lazy/file -jp /path/to/JEFF/file -esp /path/to/ensdf/folder -ep /path/to/ENDFB/file -bp /path/to/betashape/folder'
    parser = argparse.ArgumentParser(description='Collect data on the B- decaying nuclides from various databases', usage=usage)

    parser.add_argument("-lp", "--lazy_path"   , dest="lazy"   , type=str , help="path to the lazy file to write", default = None, required = True)
    parser.add_argument("-jp", "--jeff_path"   , dest="jeff"   , type=str , help="path to the jeff file", default = None, required = False)
    parser.add_argument("-jr", "--jeff_release"   , dest="jeff_r"   , type=str , help="jeff release", default = "", required = False)
    parser.add_argument("-br", "--betashape_release"   , dest="beta_r"   , type=str , help="betashape release", default = "", required = False)
    parser.add_argument("-esp", "--ensdf_path"   , dest="ensdf"   , type=str , help="path to the ensdf folder", default = None, required = False)
    parser.add_argument("-ep", "--endf_path"   , dest="endf"   , type=str , help="path to the endf-B sub-library file", default = None, required = False)
    parser.add_argument("-bp", "--betashape_path"   , dest="betashape"   , type=str , help="path to betashape directory", default = None, required = False)
    parser.add_argument("-bo", "--betashape_output"   , dest="betashape_output"   , type=str , help="path to betashape output", default = None, required = False)
    parser.add_argument("-bc", "--betashape_config"   , dest="betashape_config"   , type=str , help="input parameters for betashape", default = "myEstep=1 nu=1", required = False)
    parser.add_argument("-fix", "--ensdf_fix"   , dest="fix"   , type=int , help="try to fix the missing/theoretical spectra", default = 1, required = False)
    parser.add_argument("-ovr", "--overwrite"   , dest="overwrite"   , type=int , help="overwrite existing data", default = 1, required = False)
    parser.add_argument("-t", "--type"   , dest="type"   , type=str , help="select data_beta for evaluating the electron spectra or data_nu for evaluating the neutrino spectra", default = "data_nu", required = False)

        
    args = parser.parse_args()    
    overwrite = bool(args.overwrite)
    fix = bool(args.fix)

    
    
    LW =lazy_handler.LazyWriter(fname=args.lazy)   
    
    #Process the Jeff file
    if args.jeff is not None:
        bu.log("Processing the jeff data in "+args.jeff,level=0)
        JR =endf_reader.JeffReader(args.jeff)
        diz_JEFF =JR.get_cfy_from_list()
        
        if "0n" in diz_JEFF.keys(): #remove neutron data
            del diz_JEFF["0n"]
        bu.log(str(len(diz_JEFF.keys())) + " nuclides found with cfy greater than 0", level = 1)
        
        bu.log("Writing jeff data on " + args.lazy, level = 1)
        for key, value in diz_JEFF.items():
            LW.write_nuclide_data(nuclide_name=key,dtype="info",dictionary=value) 
        if args.jeff_r is not None:
            LW.set_general_info({"JEFF_release": args.jeff_r})
         
        bu.log("Done!", level = 1)
        
    else:
        bu.log("jeff path not given, I will assume cfy data are already present in "+args.lazy, level=0)
     
    LR =lazy_handler.LazyReader(args.lazy)  
    
    #Process the ENSDF file(s)
    if args.ensdf is not None:  
        bu.log("Processing the ensdf data in "+bu.fix_path(args.ensdf),level=0)
        lista, _ = get_ensdf_file(lname=args.lazy,ensdf_path=bu.fix_path(args.ensdf))
        diz_ensdf = {}
        R =endf_reader.EnsdfReader()
        
        for el in lista:
            nfile = bu.fix_path(args.ensdf)+el
            R.open_file(nfile)
            diz_ensdf.update(R.get_dict())   
        bu.log(str(len(diz_ensdf.keys())) + " nuclides with cfy > 0 have an associated ensdf file", level = 1)
        bu.log("I will save their Q and half-life data", level = 2)
        
        bu.log("Writing ensdf data on " + args.lazy, level = 1)
        for key, value in diz_ensdf.items():
            LW.write_nuclide_data(nuclide_name=key,dtype="info",dictionary=value) 
    else:
       bu.log("ensdf path not given, I will assume Q data are already present in "+args.lazy, level=0)  
       
       
    #Configure betashape 
    if (args.endf is not None) | (args.ensdf is not None):     
        CmdBEtashape = wrappers.CmdBetaShape()
        CmdBEtashape.set_betashape(args.betashape)

        if args.betashape_output is None:
            CmdBEtashape.set_save_path(os.getcwd())
        else:
            CmdBEtashape.set_save_path(args.betashape_output)

        bu.log("Betashape: ",level=0)
        bu.log("Betashape folder: " + args.betashape, level=1)
        bu.log("Options: " +args.betashape_config,level=1)
        bu.log("Dummy output directory: " +CmdBEtashape._save_path,level=1)
        
        Estep = float(args.betashape_config[args.betashape_config.find("myEstep=")+8:].split(" ")[0]) # in keV
        LW.set_general_info({"E_step": Estep})
        
        example_dic = CmdBEtashape.get_example_dictionary()
        
        LW.set_general_info({"Betashape_options": args.betashape_config,"Betashape_version":args.beta_r})

    
    #Process the data in the ENDF-B sub-library with betashape. 
    if args.endf is not None:  
        bu.log("Evaluating the neutrino spectrum given the data in "+args.endf,level=0)    
    
    
        ER = endf_reader.EndfBSubLibraryReader(args.endf) 
        
        symbol_list = []
        for el in periodictable.elements: 
            symbol_list.append(el.symbol)
                
        c = 1
        cmax = len(LR.get_nuclides_list())   
        for el in LR.get_nuclides_list():
            bu.log("["+str(c)+"/"+str(cmax)+"] Processing " +str(el),level=1)
            nuc =LR.get_nuclide(el)
            

            if ("dN_dE_tot" in nuc.keys()) & (overwrite is False): #FIXME
                bu.log("Aldready processed, I will skip it", level=2)
                c +=1
                continue
            
            
            if "z" not in nuc.keys():
                bu.log("No cumulative fission yield found", level=2)
                c +=1
                continue
                
            _,tag, data = ER.get_element(z=nuc["z"],a=nuc["n"]+nuc["z"],m=nuc["m"])
            
            if tag is None:
                bu.log("No data was found in the database",level=2)
                
            elif tag == "discreet":
                bu.log("Running betashape...",level=2)
                dic = create_dict(z=int(nuc["z"]),a=int(nuc["n"]+nuc["z"]),m=int(nuc["m"]),data=data,symbol_list=symbol_list,example_dic=example_dic)
                res_diz = CmdBEtashape.evaluate_decay(dic,boptions=args.betashape_config)
                
                res_list = res_diz[0]  #FIXME
                dizio = CmdBEtashape.convert_output_into_dic(res_list,tipo=args.type)
                
                bu.log("Writing data...",level=2)
                LW.write_nuclide_data(nuclide_name = nuc["lazy_name"],dtype="data",dictionary=dizio)
                LW.write_nuclide_data(nuclide_name = nuc["lazy_name"],dtype="info",vname="tag", vvalue = "endf_b")
                LW.write_nuclide_data(nuclide_name = nuc["lazy_name"],dtype="info",vname="Emax", vvalue = np.max(dizio["transition_Emax"]))
                                
                #salva i dati in continuo

            elif tag == "continuum":
                bu.log("No experimental data, just theoretical calculation",level=2)
                energy = data[:,0].astype(float)*1e3 #in keV
                nu_spectrum = np.interp(np.arange(0,energy[-1]+Estep,Estep),energy,data[:,2].astype(float))
                dizio = {"dN_dE_tot" : nu_spectrum/(np.sum(nu_spectrum)*Estep),
                         "unc_dN_dE" : np.zeros(len(nu_spectrum)),
                         "transition_Emax" : energy[-1],
                         "dN_dE_tot_c" : nu_spectrum}
                         
                bu.log("Writing data...",level=2)
                LW.write_nuclide_data(nuclide_name = nuc["lazy_name"],dtype="data",dictionary=dizio)
                LW.write_nuclide_data(nuclide_name = nuc["lazy_name"],dtype="info",vname="tag", vvalue = "endf_b_c")
                LW.write_nuclide_data(nuclide_name = nuc["lazy_name"],dtype="info",vname="Emax", vvalue = energy[-1])

            
            c+=1
    else:
        bu.log("endf path not given. Neutrino spectra will not be evaluated", level=0)
         

    #Process the data in the ENSDF database. 
    if fix is True:
        bu.log("I will try to add the missing data and to fix theoretical calc using ensdf data", level=0)  
        if args.ensdf is None:
            bu.log("Error: ensdf path not given!", level=0)
            return
         

        c = 1
        cmax = len(LR.get_nuclides_list())   
        
        ensdf_lazy_names = lazy_to_ensdf(list_of_names=LR.get_nuclides_list())
        
        for el in LR.get_nuclides_list():
            bu.log("["+str(c)+"/"+str(cmax)+"] Processing " +str(el),level=1)  
            nuc =LR.get_nuclide(el)
            
            if "z" not in nuc.keys():
                bu.log("No cumulative fission yield found", level=2)
                c +=1
                continue
                
                
            
            if "tag" not in nuc.keys():
                bu.log("Nu spectrum not found",level=2)
                                               
                if os.path.isfile(args.ensdf + lazy_to_ensdf(list_of_names=[nuc["lazy_name"]])[0]):
                    bu.log("An ensdf may file exist!", level=3)
                    res_diz =CmdBEtashape.evaluate_decay(ensdf_path=args.ensdf,ensdf_name=lazy_to_ensdf(list_of_names=[nuc["lazy_name"]])[0],boptions=args.betashape_config)
                    metastable = int(nuc["m"])
                    
                    if metastable in res_diz.keys():                      
                        res_list = res_diz[metastable]  #FIXME
                        dizio = CmdBEtashape.convert_output_into_dic(res_list,tipo=args.type)
                        bu.log("Writing data...",level=4)
                        LW.write_nuclide_data(nuclide_name = nuc["lazy_name"],dtype="data",dictionary=dizio)
                        LW.write_nuclide_data(nuclide_name = nuc["lazy_name"],dtype="info",vname="tag", vvalue = "ensdf") 
                        LW.write_nuclide_data(nuclide_name = nuc["lazy_name"],dtype="info",vname="Emax", vvalue = np.max(dizio["transition_Emax"]))
                    else:
                        bu.log("No data found for the nuclide",level=3)   
                    
                else:
                    bu.log("No ensdf file was found.", level=3)  
                    
            elif (nuc["tag"] == "ensdf") & (overwrite is False):
                bu.log("Aldready processed, I will skip it", level=2)
                c +=1
                continue
                
            elif nuc["tag"] == "endf_b_c":                    
                bu.log("Theoretical spectrum found, i will try to overwrite it with ensdf data",level=2) 
                
                if os.path.isfile(args.ensdf + lazy_to_ensdf(list_of_names=[nuc["lazy_name"]])[0]):
                    bu.log("An ensdf may file exist!", level=3)
                    res_diz =CmdBEtashape.evaluate_decay(ensdf_path=args.ensdf,ensdf_name=lazy_to_ensdf(list_of_names=[nuc["lazy_name"]])[0],boptions=args.betashape_config)
                    metastable = int(nuc["m"])
                    
                    if metastable in res_diz.keys():                      
                        res_list = res_diz[metastable]  #FIXME
                        dizio = CmdBEtashape.convert_output_into_dic(res_list,tipo=args.type)
                        bu.log("Writing data...",level=4)
                        LW.write_nuclide_data(nuclide_name = nuc["lazy_name"],dtype="data",dictionary=dizio)
                        LW.write_nuclide_data(nuclide_name = nuc["lazy_name"],dtype="info",vname="tag", vvalue = "ensdf") 
                        LW.write_nuclide_data(nuclide_name = nuc["lazy_name"],dtype="info",vname="Emax", vvalue = np.max(dizio["transition_Emax"]))
                    else:
                        bu.log("No data found for the nuclide",level=3)   
                    
                else:
                    bu.log("No ensdf file was found.",level=3)       
                    
            c += 1   
    else:
        bu.log("No fixing with ensdf data", level=0)  
    
if __name__ == "__main__":
    main()
    
   

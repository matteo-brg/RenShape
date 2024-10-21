'''
Old script, NOT used in RENSHAPE

'''

#!/usr/bin/env python
import numpy as np
from base import base_utilities as bu
from rw import lazy_handler
import h5py
from process import wrappers
from os import listdir
import os
from os import listdir
import re
import argparse
from scipy import integrate


def find_file_to_save(path=None,desin="tot_myEstep"):
    '''
    find all file inside the directory in path with desin in the name somewhere
    '''
    lista = [f.path for f in os.scandir(path)]
    temp = []
    for i in range(len(lista)):
        if lista[i].find(desin) > 0:
            temp.append(lista[i])
    return temp

def get_betashape_output(fname=None,dtype="normal"):
    '''
    get the betashape output and change it in a uniform format
    '''
    
    cmdBeta = wrappers.CmdBetaShape()
    if dtype == "normal":
        labels, data = cmdBeta.get_data(fname=fname,dtype=dtype)
        label_to_find = "dNtot/dEnu"
    if dtype == "nupartial":
        labels, _, data = cmdBeta.get_data(fname=fname,dtype=dtype)
        label_to_find = "dN/dEcalc."
    col_save = labels.index(label_to_find)
    col_save = [col_save, col_save+1]
    data_save = data[:-1,col_save]
    Emax = data[-1,0]
    Estep = data[2,0]-data[1,0]
    return data_save, Emax, np.array(labels)[col_save],Estep


def convert_nuname_betashape_lazy(nome):
    '''
    change the betashape filename to make it compatible with .lazy file
    '''
    atomic_number = re.findall(r'\d+', nome)[0]
    name = re.findall(r'\D+', nome)[0]

    metastable_index = ""
    pos = bu.findOccurrences(nome,"_")
    if len(pos) != 0:
        metastable_index = "_"+str(int(nome[pos[0]+1:])+1)+"m"
    else:
        pos = bu.findOccurrences(nome,"m")
        if len(pos) != 0:
            if pos[-1] == len(nome)-1:
                metastable_index = "_1m"
    return atomic_number+name+metastable_index
    
def calculate_full_spectrum(lname=None,nuclide_name=None,Estep=None):
    reader = lazy_handler.LazyReader(lname)
    nuclide =reader.get_nuclide(nuclide_name)
    
    spectrum_trans = nuclide["transition_dN_dE"]
    energy2 = np.arange(0,spectrum_trans.shape[1],Estep)

    for i in range(spectrum_trans.shape[0]):
        spectrum_trans[i,:] = spectrum_trans[i,:]/integrate.simps(spectrum_trans[i,:],energy2)

    spectrum_calc = (spectrum_trans.T * nuclide["transition_intensity"]).T
    spectrum_calc = np.sum(spectrum_calc,axis=0)
    spectrum_calc = spectrum_calc/integrate.simps(spectrum_calc,energy2)
    
    return spectrum_calc, np.zeros(len(spectrum_calc))
    
    
    
def get_intensity(fname = None):
    '''
    get the intensity (and its error) from the fname (betashape output)
    '''
    with open(fname,"rb") as f:    
        header = f.readline().decode(errors='replace')
        compact = header.replace(" ","")
        while compact[:6] != "E(keV)":
            text = f.readline().decode(errors='replace')
            compact = text.replace(" ","")
            header+=text
        header =header.split("\n")[:-1] 
        
        for line in header:
            lista = line.split()
            try:
                pos = lista.index("Intensity:")
                intensity = float(lista[pos+1])
                try:
                    error = float((lista[pos+2].replace("(","")).replace(")",""))
                except(IndexError):
                    error = np.nan
            except(ValueError):
                pass
        return intensity, error


def main():
    usage='runBetashape.py -f lazy_file -p ensdf_path -c config_file -o output_dir'
    parser = argparse.ArgumentParser(description='Collect all the info about the B- decaying nuclides from the live-chart with cumulative fission yield >0', usage=usage)

    parser.add_argument("-f", "--file"      , dest="fname"   , type=str , help=".lazy file"      , required = True)
    parser.add_argument("-p", "--path"      , dest="dir_beta"   , type=str , help="path to the betashape directory"      , required = True)
    parser.add_argument("-t", "--type"    , dest="type"    , type=str , help="data to save (full, partial)"  , required = True)  
    parser.add_argument("-l", "--label"    , dest="label"    , type=str , help="tag attached to each data"  , required = False,default = "ensdf")  

    args = parser.parse_args()    
    
    lname = args.fname
    betashape_out_dir = args.dir_beta

    
    #========= convert the filename inside the betashape_out_dir
    bu.log("Fixing the nuclides names in "+ betashape_out_dir+"...",level=0)
    betashape_nuclides_names =[ f.name for f in os.scandir(betashape_out_dir) if f.is_dir() ]
    lazy_nuclides_names = []
    for el in betashape_nuclides_names:
        lazy_nuclides_names.append(convert_nuname_betashape_lazy(el))
        
        
    #========= set the type of data you want to save. "tot_myEstep" = only total spectrum. "myEstep" also save the transitions.
    if args.type == "full":
        desin = "myEstep"
        bu.log("Saving both the total nu spectrum and the spectrum for each transitions (slow) ",level=0) 
    if args.type == "partial":
        desin = "tot_myEstep"    
        bu.log("Saving only the total nu spectrum (fast) ",level=0)     
        

    writer = lazy_handler.LazyWriter(fname = lname)
    reader = lazy_handler.LazyReader(lname)

    bu.log("Start merging datas ("+ betashape_out_dir+"->"+lname+"...",level=0) 
    dir_to_process = [ f.path for f in os.scandir(betashape_out_dir) if f.is_dir() ]
    bugged_dir_list = []  #list of folders without the tot_myEstep file

    imax = len(dir_to_process)
    for i in range(len(dir_to_process)):

        file_to_process = find_file_to_save(dir_to_process[i],desin=desin)
        if len(bu.locate_word(np.array(file_to_process),"tot")) == 0:  #check if the "tot_myEstep" file exist
            bugged_dir_list.append(dir_to_process[i])
            
        dic = {}
        trans_data = []
        er_trans_data = []
        trans_energy = []
        trans_intensity = []
        er_trans_intensity = []
        for file in file_to_process:
            if file.find("tot") != -1:
                dtype = "normal"
            else:
                dtype = "nupartial"
            data_mat, E_max, labels, Estep = get_betashape_output(fname=file,dtype=dtype)   
        
            if dtype == "normal":
                pos = np.where(labels == "dNtot/dEnu")[0]
                dic["dN_dE_tot"] = (data_mat[:,pos].T).flatten()
                dic["unc_dN_dE"] = (data_mat[:,pos+1].T).flatten()
            #writer.write_nuclide_data(nuclide_name = lazy_nuclides_names[i],dtype="info",vname="Emax",vvalue=E_max)
        

            if dtype == "nupartial":
                flag = False
                if "dN/dEexp." in labels:
                    pos = np.where(labels == "dN/dEexp.")[0]
                    flag = True
                else:
                    pos = np.where(labels == "dN/dEcalc.")[0]
                trans_data.append((data_mat[:,pos].T).flatten())
                er_trans_data.append((data_mat[:,pos+1].T).flatten())
                trans_energy.append(E_max)
                I, er_I = get_intensity(file.replace("_myEstep",""))
                trans_intensity.append(I)
                er_trans_intensity.append(er_I)
                temp = []
                for el in trans_data:
                    temp.append(el.shape[0])
                dim_max = max(temp)
                for j in range(len(trans_data)):
                    trans_data[j] = np.append(trans_data[j],np.zeros(dim_max-len(trans_data[j])))
                    er_trans_data[j] = np.append(er_trans_data[j],np.zeros(dim_max-len(er_trans_data[j])))
                dic["transition_dN_dE"] = np.array(trans_data)
                dic["transition_unc_dN_dE"] = np.array(er_trans_data)
                dic["transition_intensity"] = np.array(trans_intensity)
                dic["transition_unc_intensity"] = np.array(er_trans_intensity)
                dic["transition_Emax"] = np.array(trans_energy)

        bu.log("["+str(i+1)+"/"+str(imax)+"] Processing "+ dir_to_process[i]+"...",level=1)
        
        try: #FIXME #if the nuclide data is already present, skip it
            dummy = reader.get_nuclide(name=lazy_nuclides_names[i])["transition_dN_dE"]
            already_exist = True
        except:
            already_exist = False
        if already_exist == True:
            bu.log("Already has data, I will skip it!",level=2)
            continue
            
        if flag is True:
            writer.write_nuclide_data(nuclide_name = lazy_nuclides_names[i],dtype="info",vname="exp",vvalue=1) 
        writer.write_nuclide_data(nuclide_name = lazy_nuclides_names[i],dtype="info",vname="Emax",vvalue=E_max)
        writer.write_nuclide_data(nuclide_name = lazy_nuclides_names[i],dtype="info",vname="label", vvalue = args.label)
        writer.write_nuclide_data(nuclide_name = lazy_nuclides_names[i],dtype="data",dictionary=dic)


    bu.log("Done!",level = 0)
    
    if len(bugged_dir_list) != 0:
        bu.log("Warning: the following folders do NOT have the _tot_myEstep_ file", level=0)
        for f in bugged_dir_list:
            bu.log(f,level=1)
        if args.type == "partial":
            bu.log("If you want to evaluate them from the transition files, launch this script with -t full option", level=2)
        if args.type == "full":
            bu.log("I will try to fix it, evaluate the total spectrum from the transition files...", level=2)
            for i in range(len(bugged_dir_list)):
                bu.log("Fixing "+bugged_dir_list[i], level=3)
                nuname = convert_nuname_betashape_lazy( bugged_dir_list[i][bu.findOccurrences(bugged_dir_list[i],"/")[-1]+1:] )
                try:
                    spectrum_total, dummy_error = calculate_full_spectrum(lname,nuname,Estep)
                except Exception as error:
                    bu.log("Some error occurred! -> ", level=4)
                    print(error)
                    continue
                writer.write_nuclide_data(nuclide_name = nuname,dtype="data",vname="dN_dE_tot",vvalue=spectrum_total)
                writer.write_nuclide_data(nuclide_name = nuname,dtype="data",vname="unc_dN_dE",vvalue=dummy_error)
                writer.write_nuclide_data(nuclide_name = nuname,dtype="info",vname="Emax",vvalue=np.arange(0,len(spectrum_total),Estep)[-1])
                writer.write_nuclide_data(nuclide_name = nuname,dtype="info",vname="fixed",vvalue=1)
                
    bu.log("Fixing the lazy file header...",level=0)
    reader = lazy_handler.LazyReader(lname)
    labels = reader.get_info()["labels"]
    writer.set_general_info({"E_step":Estep})
    writer.set_general_info({"labels":labels+", Emax, fixed, label"})
    bu.log("All done!",level=0)
    return 
    

    


if __name__ == "__main__":
    main()

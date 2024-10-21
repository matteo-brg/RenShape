'''
Old script, NOT used in RENSHAPE

'''

#!/usr/bin/env python
import numpy as np
from base import base_utilities as bu
from rw import lazy_handler
from process import wrappers
from os import listdir
import os
import argparse



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
    if lname is None:
        lista = np.array(listdir(ensdf_path))
        temp = []
        for i in range(len(lista)):
            if lista[i].find(dtype) > 0:
                temp.append(lista[i])
        return np.array(temp)
    else:
        reader = lazy_handler.LazyReader(lname)
        nuclides_for_betashape = lazy_to_ensdf(reader.get_nuclides_list())
        path_files = np.array(listdir(ensdf_path))
        
        common_files = np.intersect1d(nuclides_for_betashape,path_files)
        not_present =  np.array(list(set(nuclides_for_betashape) - set(common_files)))
        return common_files, not_present



def main():
    usage='runBetashape.py -f lazy_file -p ensdf_path -c config_file -o output_dir'
    parser = argparse.ArgumentParser(description='Collect all the info about the B- decaying nuclides from the live-chart with cumulative fission yield >0', usage=usage)

    parser.add_argument("-f", "--file"      , dest="fname"   , type=str , help=".lazy file"      , required = False, default=None)
    parser.add_argument("-p", "--path"      , dest="dir_ensdf"   , type=str , help="path to the ensdf directory"      , required = True)
    parser.add_argument("-c", "--config"    , dest="config"   , type=str , help="config file"      , required = True)
    parser.add_argument("-o", "--output"    , dest="save_path"    , type=str , help="output directory"  , required = True)  
    parser.add_argument("-overwrite", "--overwrite"    , dest="overwrite"    , type=str , help="overwrite files"  , required = False, default=False)      
    args = parser.parse_args()    
    
    lname = args.fname
    ensdf_path = args.dir_ensdf
    dtype = ".ensdf"
    config = args.config
    betashape_out_dir = args.save_path
    config_dic = bu.Config2Dict(config).get_parameters()
    
    bu.log("Finding "+dtype+" files to process with betashape in "+ensdf_path+" ...",level=0)
    fname_to_process, file_not_found = get_ensdf_file(lname=lname,ensdf_path=ensdf_path)  #list of file present both in betashape dir and in the .lazy file.
    imax = len(fname_to_process)
    bu.log(str(imax)+" files to process.", level=1)
    
    bu.log("Setting up betashape...",level=0)
    cmdBeta = wrappers.CmdBetaShape(config_dic["path_to_betashape"])
    bu.log("Creating or setting the output dir in "+betashape_out_dir,level=1)
    os.system('mkdir -p' + betashape_out_dir)
    cmdBeta.set_save_path(betashape_out_dir)
    
    bu.log("Setting the required options...",level=1)
    cmd_for_betashape = str(config_dic["betashape_opt"])

    f_not_processed = []
    f_processed = []
    f_not_found = []    
        
    listina =[f.name for f in os.scandir(betashape_out_dir) if f.is_dir() ]
    if len(listina) != 0:
        for i in range(len(listina)):
            listina[i] = listina[i].upper() + ".ensdf"    
        
    bu.log("Running betashape...",level=1)
    imax = len(fname_to_process)
    for i in range(imax):
        bu.log("["+str(i+1)+"/"+str(imax)+"] Processing "+ fname_to_process[i]+"...",level=2)
        if ((fname_to_process[i] in listina) & (args.overwrite is False)):
            bu.log("File already exist! Skip!",level=3)
            continue
        cmdBeta.run_betashape(fpath=ensdf_path, fname=fname_to_process[i], options=cmd_for_betashape,verbose=False)
        bu.log(cmdBeta.get_result_state(),level=3)
        
        if cmdBeta.get_result_state() == "file processed":
            f_processed.append(fname_to_process[i])
        if cmdBeta.get_result_state() == "file not processed":
            f_not_processed.append(fname_to_process[i])
    
    text_log = "# \n" + "# ENSDF files NOT processed: \n"
    text_log += ", ".join(f_not_processed) + "\n"
    text_log += "# \n" + "# ENSDF files NOT found: \n"
    text_log += ", ".join(list(file_not_found)) + "\n"
    text_log += "# \n" + "# ENSDF files found: \n"
    text_log += ", ".join(f_processed) + "\n"    
    bu.log("Creating .log file...",level=0)
    with open(bu.fix_path(betashape_out_dir)+"README"+".log","w") as f:
        f.write(text_log)
    return 
    

    


if __name__ == "__main__":
    main()

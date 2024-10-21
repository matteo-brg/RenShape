'''
Old script, NOT used in RENSHAPE

'''

#!/usr/bin/env python
import numpy as np
import warnings
warnings.filterwarnings('ignore') #FIXME

from base import base_utilities as bu
from rw import lazy_handler
import h5py
from process import wrappers
import argparse
from datetime import date

def find_cumulative_fission_yield(parent_list = None):
    '''
    find all the nuclides with a cumulative fission yield > 0
    '''
    dictionary = {}
    metastable_nuclides_names = np.array([])  
    non_metastable_nuclides_names = np.array([])    
    variable_names = ["cumulative_thermal_fy","unc_ct","cumulative_fast_fy","unc_cf"]
    
    for parent in parent_list:
        nuchart = wrappers.CmdNuChart()
        nuchart.get_fission_yields(kind="cumulative",parent=parent)
        '''
        variable_names = ["cumulative_thermal_fy","unc_ct"]    
        if parent == "238u":
            variable_names = ["cumulative_fast_fy","unc_cf"]
        '''
        for kk in range(len(variable_names)//2):
        
            #========= find the cumulative fission yield >0
            ii = nuchart.get_index_non_empty(variable_names[kk*2])
            cfy= nuchart.get_array(variable_names[kk*2])[ii].astype(float)
            er_cfy = nuchart.get_array(variable_names[kk*2+1])[ii].astype(float)

            #========= write the name of the nuclides as mass number + symbol
            nuclides_names = (nuchart.get_array("a_daughter")[ii].astype(str)+nuchart.get_array("element_daughter")[ii]).astype(str)
            nuclides_names = nuclides_names.astype("<U8")
            nuclides_names_original = nuclides_names

            #========= fix the name of the metastable nuclides. e.g. 166Ho_1m
            uniq, index, counts =np.unique(nuclides_names,return_index=True,return_counts=True)
            index_non_unique =np.uint32(np.where(np.in1d(np.arange(len(nuclides_names)),index,invert = True))[0])

            index_metastable = index_non_unique   
            index_normal =np.uint32(np.where(np.in1d(np.arange(len(nuclides_names)),index_metastable,invert = True))[0])
    
            #update the array
            metastable_nuclides_names = np.append(metastable_nuclides_names,nuclides_names_original[index_metastable])
            non_metastable_nuclides_names = np.append(non_metastable_nuclides_names,nuclides_names_original[index_normal])
 
            i = 1
            while len(index_non_unique) != 0:
    
                for j in index_non_unique:
                    if i == 1:
                        nuclides_names[j] = nuclides_names[j]+"_"+str(i)+"m"
                    else:
                        nuclides_names[j] = nuclides_names[j][:-3]+"_"+str(i)+"m"
                i += 1
                uniq, index =np.unique(nuclides_names,return_index=True)
                index_non_unique =np.uint32(np.where(np.in1d(np.arange(len(nuclides_names)),index,invert = True))[0])
                
            #========= update the dictionary with the info found 
            j = 0
            for nuclide in nuclides_names:
                temp = {}
                temp[variable_names[kk*2]+"_"+parent] = cfy[j]
                temp[variable_names[kk*2+1]+"_"+parent] = er_cfy[j]
                j +=1
                try:
                    dictionary[nuclide].update(temp)
                except:
                    dictionary[nuclide] = temp

    return dictionary



def find_nuclide_beta(E_thr = None,nuclides_names=None,names_lcn=None,names=None):
    '''
    find all the beta-decaying nuclides with a Q-value >= E_thr
    '''
    dictionary = {}
    nuchart = wrappers.CmdNuChart()
    nuchart.get_nu_chart()  #return the ground states
    
    atomic_mass = (nuchart.get_array("z")+nuchart.get_array("n")).astype(str)
    symbol = nuchart.get_array("symbol").astype(str)
    nuclides_names_total = np.char.add(atomic_mass,symbol)

    #========= find the position of all the beta-decaying nuclides (iiB) with a Q value greater than E_thr (iiQbm)
    # the nuclides which satisfy all the requiremets have position (ii_thr)

    iiB = np.where((nuchart.get_array("decay_1")=="B-")|(nuchart.get_array("decay_2")=="B-")|(nuchart.get_array("decay_3")=="B-"))[0]

    Qbm = nuchart.get_array("qbm")[nuchart.get_index_non_empty(column_name="qbm")].astype(float)
    iiQbm = nuchart.get_index_non_empty(column_name="qbm")[np.where(Qbm>E_thr)[0]]

    ii_thr = np.intersect1d(iiB,iiQbm)
  
    # [nuclides_names_beta_thr] are the names of the nuclides: 1) B- decay 2) with a cumulative fission yield  >=0 3) with a Q >= E_Thr
    nuclides_names_beta_thr, _, index_ar2 = np.intersect1d(nuclides_names,nuclides_names_total[ii_thr],return_indices=True) 
    
    #========= update the dictionary with the info
    data_mat = nuchart.get_array(names_lcn)[ii_thr[index_ar2]]

    j = 0
    for nuclide in nuclides_names_beta_thr:
        temp = dict(zip(names, data_mat[j,:]))
        dictionary[nuclide] = temp
        j += 1
    return dictionary



def find_nuclide__meta_beta(E_thr = None,nuclides_names=None,
                            names_lcn1=None,names1=None,
                            names_lcn2=None,names2=None,
                            verbose=True):
    '''
    find all the beta-decaying nuclides metastable with a Q-value >= E_thr
    '''
    dictionary = {}
    nuchart = wrappers.CmdNuChart()

    i_max = len(nuclides_names)

    for i in range(i_max):
        nuclide_name = nuclides_names[i][:-3] #delete the _nm for livechart compatibility
        if verbose is True:
            print("[",i+1,"/",i_max,"] ", nuclide_name, " |-> ", nuclides_names[i])
        flag = True
        while flag is True:
            try:
                nuchart.get_nu_chart(nuclide_name)
                flag = False
            except:
                print("Waiting for server response...")
                
        #save the usual parameters for the metastable nuclide from the ground state
        temp = dict(zip(names1, nuchart.get_array(names_lcn1)[0]))
        temp["metastable"] = int(nuclides_names[i][-2])
        
        if temp[names1[0]] == '':  #skip the ones without a Q value for B- decay
            continue
    
        #add the parameters for the metastable nuclide from the levels table
        flag = True
        while flag is True:
            try:
                nuchart.get_nuclear_levels(nuclide_name)
                flag = False
            except:
                print("Waiting for server response...")
                

        # find the location of proper the energy shift for the corresponding metastable state 
        # |-> [pos_beta][temp["metastable"]])
            
        pos_beta = np.where((nuchart.get_array("decay_1")=="B-")|(nuchart.get_array("decay_2")=="B-")|(nuchart.get_array("decay_3")=="B-"))[0]
        try:
            array_to_save =nuchart.get_array(names_lcn2)[pos_beta][temp["metastable"]]
        except(IndexError):    #if the metastable nuclide is not beta-decaying, skip it
            continue
        
        dict2 = dict(zip(names2, array_to_save)) #add the new parameters for the metastable nuclide with the proper name
        temp.update(dict2)
    
        #handle ripl_shift FIXME
        try:
            if (float(temp[names1[0]]) + float(temp[names2[0]])) >= E_thr:
                dictionary[nuclides_names[i]] = temp
            else:
                pass
        except:
            temp["ripl_shift"] = 1
            dictionary[nuclides_names[i]] = temp
    return dictionary 




def main():
    
    usage='script.py -p save_path -n name -c config_path -not "write here some notes if you want"'
    parser = argparse.ArgumentParser(description='Collect all the info about the B- decaying nuclides from the live-chart with cumulative fission yield >0', usage=usage)

    parser.add_argument("-p", "--path"   , dest="save_path"   , type=str , help="Output path"      , required = True)
    parser.add_argument("-n", "--name"    , dest="fname"    , type=str , help="name of the .lazy file (hdf5)"  , required = True)
    parser.add_argument("-c", "--config"    , dest="config"    , type=str , help="path to the .conf file"  , required = True) 
    parser.add_argument("-not", "--notes"    , dest="note"    , type=str , help="some notes"  , required = None, default = '') 
    
    args = parser.parse_args()    
    
    #========= Set the initial parameters
    dic_config_file = bu.Config2Dict(args.config)
    dic_config_file = dic_config_file.get_parameters()
    
    save_path = bu.fix_path(args.save_path)
    
    E_thr = dic_config_file["E_thr"]
    name_to_save_lcn = ["qbm","unc_qb","z","n","symbol"]
    name_to_save = ["Q","unc_Q","z","n","symbol"]
    name_to_save_meta_lcn = ["energy","unc_e"]
    name_to_save_meta = ["E_excited","unc_E_excited"]
    
    if dic_config_file["variable_to_save_lcn"] != '':
        name_to_save_lcn.extend(dic_config_file["variable_to_save_lcn"])
        name_to_save.extend(dic_config_file["variable_to_save"])
    
    if dic_config_file["variable_to_save_meta_lcn"] != '':
        name_to_save_meta_lcn.extend(dic_config_file["variable_to_save_meta_lcn"])    
        name_to_save_meta.extend(dic_config_file["variable_to_save_meta"])  
           
    #========= Cumulative fission yield
    bu.log("Collecting the cumulative fission yield of 235u, 238u, 239Pu and 241 Pu...", level=0) 
    dictionary = find_cumulative_fission_yield(parent_list = ["235u","238u","239Pu","241Pu"])
    bu.log("Done!", level=1)

    dic_keys = np.array(list(dictionary.keys()))
    nuclides_names_non_meta = dic_keys[np.uint32(np.where(np.in1d(np.arange(len(dic_keys)),bu.locate_word(dic_keys,"_"),invert = True))[0])]
    nuclides_names_meta = dic_keys[bu.locate_word(dic_keys,"_")]
    bu.log(str(len(dic_keys)) + " nuclides found!",level=1)
    bu.log(str(len(nuclides_names_non_meta))+ " nuclides found in a non metastable state",level=2)
    bu.log(str(len(nuclides_names_meta)) + " nuclides found in a metastable state",level=2)

    
    #========= B- non metastable
    bu.log("Collecting B- nuclides (non metastable) with Q >= "+ str(E_thr) + " [keV]...",level=0)
    dictionary2 = find_nuclide_beta(E_thr = E_thr,nuclides_names=nuclides_names_non_meta,names_lcn=name_to_save_lcn,names=name_to_save)
    nuclides_deleted2 = nuclides_names_non_meta[np.uint32(np.where(np.in1d(nuclides_names_non_meta,np.array(list(dictionary2.keys())),invert = True))[0])]
    bu.log("Done!",level=1)
    bu.log(str(len(nuclides_deleted2) )+ " nuclides (non metastable) previously found do not satisfy the required conditions.",level=2)
    
    
    #========= B- metastable
    bu.log("Collecting B- nuclides (metastable) with Q >= "+ str(E_thr) + " [keV]...",level=0)
    dictionary3 = find_nuclide__meta_beta(E_thr = E_thr,nuclides_names=nuclides_names_meta,
                        names_lcn1=name_to_save_lcn,names1=name_to_save,
                        names_lcn2=name_to_save_meta_lcn,names2=name_to_save_meta)
    nuclides_deleted3 = nuclides_names_meta[np.uint32(np.where(np.in1d(nuclides_names_meta,np.array(list(dictionary3.keys())),invert = True))[0])]
    bu.log(str(len(nuclides_deleted3) )+ " nuclides (metastable) previously found do not satisfy the required conditions.",level=2)
        
    #========= Create one dictionary
    bu.log("Merging data...",level=0)
    dictionary_final = bu.merge_dictionaries(dictionary,dictionary2)
    dictionary_final = bu.merge_dictionaries(dictionary_final,dictionary3)
    
    bu.log("Cleaning the useless nuclides...",level=1)

    temp = []
    for key in dic_keys:
        if name_to_save[0] in dictionary_final[key]:
            pass
        else:
            temp.append(key)
            dictionary_final.pop(key)
    bu.log("In the end, only "+ str(len( np.array(list(dictionary_final.keys()))))+ "/"+ str(len(dic_keys))+ " nuclides remain.",level=2)    
    #========= Create .log file
    
    dic_keys = np.array(list(dictionary_final.keys()))
    text_log = "# These are nuclides name to use to search the appropriate ENSDF file: \n"
    text_log += ", ".join(dic_keys[np.uint32(np.where(np.in1d(np.arange(len(dic_keys)),bu.locate_word(dic_keys,"_"),invert = True))[0])]) + '\n'
    text_log += "# \n" + "# These are nuclides name in the .lazy file: \n"
    text_log += ", ".join(dic_keys) + "\n"
    text_log += "# \n" + "# These are nuclides non metastable which have cumulative fission yield > 0 but do not satisfy the required conditions.: \n"
    text_log += ", ".join(nuclides_deleted2) + "\n"
    text_log += "# \n" + "# These are nuclides metastable which have cumulative fission yield > 0 but do not satisfy the required conditions.: \n"
    text_log += ", ".join(nuclides_deleted3) + "\n"
    bu.log("Creating .log file...", level=0)
    with open(save_path+args.fname+".log","w") as f:
        f.write(text_log)
    bu.log("Done!",level=1)

    #========= Save the dictionary in a .lazy (hdf5) file
    bu.log("Saving the .lazy file in "+ save_path + "...",level=0)
    writer = lazy_handler.LazyWriter(output_path=save_path,name=args.fname)
    
    name_to_save.extend(name_to_save_meta)    
    name_to_save.extend(["metastable","ripl_shift"])
    dic_header = {}
    dic_header["labels"] =  ", ".join(name_to_save)
    dic_header["date"] = date.today().strftime("%m/%d/%y")
    if args.note != "":
        dic_header["note"] = args.note
    dic_header["threshold_energy"] = E_thr
    writer.set_general_info(dic_header)
    
    for key, value in dictionary_final.items():
        writer.write_nuclide_data(key,"info",value)
    bu.log("All Done!",level=0)

if __name__ == "__main__":
    main()

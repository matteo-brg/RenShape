'''
 Desc  : Wrapper for Betashape and the Livechart Data API
 Author: Matteo Borghesi <matteo.borghesi@mib.infn.it>
'''

import warnings
warnings.filterwarnings('ignore') #FIXME

import pandas as pd
import urllib.request
import numpy as np
from base import base_utilities as bu 
from os import listdir
import subprocess
import os
import re



class CmdNuChart(object):
    '''
    Wrapper for https://www-nds.iaea.org/relnsd/vcharthtml/api_v0_guide.html
    Currently not used in RENSHAPE
    
    '''
    
    def __init__(self,http="https://nds.iaea.org/relnsd/v0/data?"):
        self._url = http
        self._df = None
        
    def _return_csv(self,url):
        req = urllib.request.Request(url)
        req.add_header('User-Agent', 'Mozilla/5.0 (X11; Ubuntu; Linux x86_64; rv:77.0) Gecko/20100101 Firefox/77.0')
        self._df = pd.read_csv(urllib.request.urlopen(req))
        return
    
    def get_df(self):
        return self._df
    
    def _ask_df(self,url):
        url = self._url + url
        self._return_csv(url)
        return

    def get_nu_chart(self,nuclide="all"):
        self._ask_df("fields=ground_states&nuclides="+nuclide)
        return
    
    def get_nuclear_levels(self,nuclide=None):
        self._ask_df("fields=levels&nuclides="+nuclide)
        return    
    
    def get_fission_yields(self,kind='cumulative',parent=None):
        if kind == 'cumulative':
            name = "fields=cumulative_fy"
        if kind == "independent":
            name = "fields=independent_fy"
        self._ask_df(name+"&parents="+parent)
        return
        
    def get_df_labels(self):
        return list(self._df.columns)
    
    def get_number_of_elements(self):
        return len(self.get_df())
    
    def get_index_non_empty(self,column_name = None):
        return np.where(np.isnan(np.array(self.get_df()[column_name])) == False)[0]
    
    def get_array(self,column_name=None):
        data = np.array(self.get_df()[column_name])

        if len(data.shape) == 2:
            for col in range(data.shape[1]):
                #data[:,col][data[:,col]==" "] = np.nan
                try:
                    data[:,col] = data[:,col].astype(float)
                except:
                    data[:,col] = data[:,col].astype(str)
                    
        if len(data.shape) == 1:
            for col in range(data.shape[0]):
                try:
                    data[col] = float(data[col])
                except:
                    data[col] = str(data[col])    
                            
        return data
    
    def get_decays(self,nuclide=None,dtype="bm"): #135xe
        self._ask_df("fields=decay_rads&nuclides="+nuclide+"&rad_types="+dtype)
        return
        
        
class CmdBetaShape(object):
    '''
    Wrapper for Betashape. It allows to use Betashape directly from python with a simple
    I/O. It also enable Betashape to precess a generic python dictionary.
    
    '''

    def __init__(self,betashape_path=None):
        self._betashape_path = betashape_path
        self._save_path = "."
        self._message = None
        
    def set_betashape(self,path=None):
        self._betashape_path = path
        return
        
    def set_save_path(self,path=None):
        self._save_path = path
        return
        
    def get_example_dictionary(self):
        dic = {"P": "12B",
           "D": "12C",
           "E_P": 0.0,
           "dE_P": "",
           "JP_P" : "1+",
           "half_life": 20.20,
           "half_life_units": "MS",
           "dhalf_life": 2,
           "Q": 13369.4,
           "dQ": 13,
           "NR": 1.,
           "dNR": "",
           "BR": 1,
           "dBR": "",
           "E_D": [0],
           "dE_D": [""],
           "JP_D": ["0+"],
           "L": "",
           "IB": [98.216],
           "d_IB": [28],
           "trans_type":""}
        return dic
        
    def run_betashape(self,fpath=None,fname = None, options = "",verbose=True):
        bpath = bu.fix_path(self._betashape_path)
        original_files = [f for f in listdir(bpath)] # name of the files in the betashape folder BEFORE launching betashape
        
        file_to_process = bu.fix_path(fpath) + fname  #copy the file to process inside the betashape folder
        subprocess.call(["cp",file_to_process, bpath])
        
        cmd = "betashape "+fname + " " + options  #launch betashape
        message =subprocess.Popen(cmd,stdout=subprocess.PIPE,stderr=subprocess.PIPE,shell=True,cwd=self._betashape_path).communicate()
        if verbose is True:
            print(message[0].decode('ascii'))
            
        full_files = [f for f in listdir(bpath)] # name of the files in the betashape folder AFTER launching betashape
        ii =np.uint32(np.where(np.in1d(np.array(full_files),np.array(original_files),invert = True))[0])
        new_files = np.array(full_files)[ii]     # file created by betashape
        
        file_dir = None
        for filen in new_files:
            if os.path.isdir(bpath+filen) is True:
                file_dir = filen
            subprocess.call(["mv",  bpath+filen, self._save_path])
        self._message = message[0].decode('ascii')
        
        '''
        new_files = np.delete(new_files,new_files==file_dir)  
        if file_dir is not None:    #copy ALL the betashape output inside the nuclide folder (file_dir).
            for filen in new_files:
                subprocess.call(["mv",  self._save_path+"/"+filen, self._save_path+"/"+file_dir])
        '''
        return  
        
    def get_result_state(self,thr=800):
        if self._message is None:
            return "Not executed"
        if len(self._message) == 0:
            return "file not found"
        if len(self._message) >= thr:
            return "file processed"
        if len(self._message) < thr:
            return "file not processed"
        
    def __read_header(self,fname=None):
        with open(fname,"r",encoding="ISO-8859-1") as f:    
            nrow = 1
            header = ""
            try:
                header += f.readline()
            except:
                pass
            compact = header.replace(" ","")
            while compact[:6] != "E(keV)":
                text = f.readline()
                compact = text.replace(" ","")
                header+=text
                nrow +=1
            header =header.split("\n")[:-1]
            labels = header[-1]
            labels = labels.split("  ")
            labels_true = [x for x in labels if( (x != '') & (x != ' '))]
            for i in range(len(labels_true)):
                labels_true[i] = labels_true[i].replace(" ","")
        return header, labels_true, nrow
        
        
    def create_dummy_ensdf(self,dic,path,debug=False):
        D_number = re.findall(r'\d+', dic["D"])[0]
        D_nuclide_record = "{0:>3}".format(D_number)+dic["D"][len(D_number):]   
        
        P_number = re.findall(r'\d+', dic["P"])[0]
        P_nuclide_record = "{0:>3}".format(P_number)+dic["P"][len(P_number):] 
        
        identification_record = D_nuclide_record + bu.entry(" ",9-len(D_nuclide_record)) + bu.entry(dic["P"]+ " B- DECAY:"+str(dic["half_life"])+ " "+dic["half_life_units"],30) + "\n" #FIXME

        parent_record =P_nuclide_record + bu.entry(" ",7-len(P_nuclide_record))+"P "+ bu.entry(dic["E_P"],10)+bu.entry(dic["dE_P"],2)
        parent_record += bu.entry(dic["JP_P"],18) + bu.entry(str(dic["half_life"])+ " "+ dic["half_life_units"],10) + bu.entry(dic["dhalf_life"],6)
        parent_record += bu.entry("",9) +bu.entry(dic["Q"],10) + bu.entry(dic["dQ"],2) +bu.entry("",4)+"\n"

        normalization_record = D_nuclide_record+ bu.entry(" ",7-len(D_nuclide_record))+"N " + bu.entry(dic["NR"],9)+bu.entry(dic["dNR"],2)+bu.entry("",11)
        normalization_record += bu.entry(dic["BR"],8)+bu.entry(dic["dBR"],2)+bu.entry(1.0,10)+"\n"

        production_norm_record = D_nuclide_record+ bu.entry(" ",6-len(D_nuclide_record))+ "PN" +"\n"
    
        dummy_ensdf = identification_record+parent_record+normalization_record+production_norm_record
        
        for i in range(len(dic["E_D"])):
            level_record = D_nuclide_record+ bu.entry(" ",7-len(D_nuclide_record))+"L "+ bu.entry(dic["E_D"][i],10)+bu.entry(dic["dE_D"][i],2)+ bu.entry(dic["JP_D"][i],18)
            level_record += bu.entry("",16) + bu.entry(dic["L"],9) + "\n"
        
            beta_record = D_nuclide_record+ bu.entry(" ",7-len(D_nuclide_record))+ "B " + bu.entry("",12)+bu.entry(dic["IB"][i],8)+bu.entry(dic["d_IB"][i],2)+bu.entry("",46)
            beta_record += bu.entry(dic["trans_type"],2) + "\n"
            
            dummy_ensdf += level_record + beta_record
        
        if debug is True:
            return dummy_ensdf
        
        with open(path+'/dummy.ensdf', 'w') as f:
            f.write(dummy_ensdf)
        return      


    def evaluate_decay(self,dictionary=None,ensdf_path=None,ensdf_name=None,verbose=False,boptions="myEstep=1 nu=1",rmdir=True):
        
        folder_name = "dummyFolder101"
        path = bu.fix_path(self._save_path)
        original_files = [f for f in listdir(path)] 
        if folder_name in original_files:
            print("Error!" ,folder_name," is present! Please remove it!")
            return
        message =subprocess.Popen("mkdir "+folder_name,stdout=subprocess.PIPE,stderr=subprocess.PIPE,shell=True,cwd=path).communicate()
        
        old_save_path = self._save_path
        self.set_save_path(path+folder_name)
        if dictionary is not None:
            self.create_dummy_ensdf(dictionary,path+folder_name)
            self.run_betashape(fpath=path+folder_name,fname="dummy.ensdf",options=boptions,verbose=verbose)
        if ensdf_path is not None:
            self.run_betashape(fpath=ensdf_path,fname=ensdf_name,options=boptions,verbose=verbose)
        #self.run_betashape(fpath="./"+folder_name,fname="dummy.ensdf",options="myEstep="+str(Estep)+ " nu=1 "+boptions,verbose=verbose)
    
        full_files = [f for f in listdir(path+folder_name)] 
        file_dir = []
        for filen in full_files:
            if os.path.isdir(path+folder_name+"/"+filen) is True:
                file_dir.append(filen)

        output_dir = {}
        for el in file_dir:        
            output_list = self.get_data_from_folder(folder_name=path+folder_name+"/"+el)
            output_dir[self.__find_metastable(el)] = output_list
            

        
        if rmdir is True:
            message =subprocess.Popen("rm -r "+folder_name,stdout=subprocess.PIPE,stderr=subprocess.PIPE,shell=True,cwd=path).communicate()
    
        self.set_save_path(old_save_path)
        return output_dir
        
    def get_data_from_folder(self,folder_name=None):
        full_files = [f for f in os.listdir(folder_name)]    
        output_list = []
        pos =np.intersect1d(bu.locate_word(np.array(full_files),"trans"),bu.locate_word(np.array(full_files),"myEstep"))
        for i in pos:
            dic = {}
            file_to_open = folder_name+ "/" + full_files[i]
            label, data_beta, data_nu = self.get_data(file_to_open,dtype="nupartial")

            output_message = self.__read_header(file_to_open.replace("_myEstep","") )[0] 
            dic["label"] = label
            dic["data_beta"] = data_beta
            dic["data_nu"] = data_nu
            dic["output_message"] = output_message
            output_list.append(dic)    
        return output_list
        
    def __find_transition_type(self,list_of_lines):

        tipo = ""
        for line in list_of_lines:
            if "allowed" in line:
                tipo = "a"    
            if "Calculation of" in line:
                if "allowed" in line:
                    tipo = "a"
                elif "forbidden unique":
                    number = re.findall(r'\d+', line)[0]
                    tipo = str(number) + "u"
            if "forbidden non-unique" in line:
                number = re.findall(r'\d+', line)[0]
                tipo = str(number) + "nu"
        return tipo      
        
    def convert_output_into_dic(self,output_list=None,tipo="data_nu"):
    
        Emax = []
        transition_intensity = []
        transition_dN_dE = []
        transition_unc_dN_dE= []
        transition_intensity = []
        transition_unc_intensity = []
        transition_type = []
    
        for diz in output_list:
            key = "dN/dEcalc."
            pos = np.where(np.char.find(diz["label"],key)==0)[0][0]
            transition_dN_dE.append(diz[tipo][:,pos][:-1])
            Emax.append(diz[tipo][-1,0])
            transition_unc_dN_dE.append(diz[tipo][:,pos+1][:-1])
    
            for line in diz["output_message"]:
                word = "Intensity: "
                pos= line.find(word)
                if pos!= -1:
                    new_line = line[pos+len(word):]
                    pos2=new_line.find("(")
                    if pos2 == -1:
                        I = float(new_line)
                        dI = 0
                    else:
                        I = float(new_line[:pos2])
                        dI = int(new_line[pos2+1:new_line.find(")")])          
            transition_intensity.append(I)
            transition_unc_intensity.append(dI)
            transition_type.append(self.__find_transition_type(diz["output_message"]))
    
        Emax = np.array(Emax)
        transition_unc_intensity = np.array(transition_unc_intensity)
        transition_intensity = np.array(transition_intensity)
        
        temp = np.zeros((len(transition_dN_dE),max([el.shape[0] for el in transition_dN_dE])))
        temp2 = temp.copy()
        for i in range(temp.shape[0]):
            temp[i,:len(transition_dN_dE[i])] = transition_dN_dE[i] 
            temp2[i,:len(transition_unc_dN_dE[i])] = transition_unc_dN_dE[i] 
        transition_dN_dE = temp
        transition_unc_dN_dE = temp2
    
        transition_dN_dE_total = np.sum(transition_dN_dE,axis=0)
        unc_dN_dE = np.sum(transition_unc_dN_dE**2,axis=0)**0.5
    
        dizio ={"transition_Emax" : Emax,
                "transition_intensity" : transition_intensity,
                "transition_dN_dE" : transition_dN_dE,
                "transition_unc_dN_dE" : transition_unc_dN_dE,
                "transition_intensity" : transition_intensity,
                "transition_unc_intensity" : transition_unc_intensity,
                "dN_dE_tot" : transition_dN_dE_total,
                "unc_dN_dE": unc_dN_dE,
                "transition_type" : transition_type
                }
        return dizio
        
        
    def get_data(self,fname=None,dtype="normal"):
        header,labels,nrow = self.__read_header(fname)
        
        if dtype == "normal":
            data = np.loadtxt(fname,skiprows=nrow)
            return labels, data
        if dtype == "nupartial":
            with open(fname,"r") as f:    
                a = f.readlines()
                a = np.array(a)
                index = np.where(a=='\n')[0]

            nrows1 = index[index>nrow][0] - nrow
            data1 =  np.loadtxt(fname,skiprows=nrow,max_rows=nrows1)
            nrow2 = index[index>nrow][-1] +2   
            data2 =  np.loadtxt(fname,skiprows=nrow2)
            return labels, data1, data2
    
    def __find_metastable(self,name):
        meta = None
        if name.find("_") != -1:
            meta = int(name[name.find("_")+1:]) +1
            return meta
        
        if name[-1] == "m":
            return 1
        else:
            return 0

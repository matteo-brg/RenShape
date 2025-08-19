'''
 Desc  : Different readers for different nuclear database
 Author: Matteo Borghesi <matteo.borghesi@mib.infn.it>
'''

import numpy as np
import periodictable
import re
import uncertainties


def zaid(z,a):
    return 1000*z+a

def element(zaid):
    z = int(zaid/1000)
    a = int(zaid-1000*z)
    return z,a, str(a)+periodictable.elements[z].symbol

def convert(stringa):
    if (("E" not in stringa) & ("e" not in stringa)) :
        stringa = stringa.replace("+","e+")
        stringa = stringa.replace("-","e-")
    return float(stringa)


class EndfBReader:
    
    '''
    Simple ENDF reader. It should extract all the relevant information for evaluating the B- decay.
    '''
    
    def __init__(self, path=None):
        

        self._path = path
        self._endf = None
        self._endf_total = None
        self._dic = {}
        self._MAT = []
        self._MF = []
        self._MT = []
        self._pos_section = []
        self._nline = None
        if self._path is not None:
            self.open_file()
        return

   
    def open_file(self,path=None):
        if path is not None:
            self._path = path
        with open(self._path,"r") as f:
            self._endf_total = f.readlines()
        return  

    def get_element(self,name,verbose=False):
        '''
        name : string
            name must be in the following format: Element+mass number. e.g. Si35. 96YM (for metastable state)
        '''
        line = []
        i = 0
        flag = False
        for el in self._endf_total:
            endmarker = el[71:80]
            if el[3:11].replace(" ","").replace("-","") == name:
                line.append(i-4)
                i +=1
                flag = True
            if ((endmarker == "8  099999") & (flag == True)):
                line.append(i-1)
                break
            i +=1
        if len(line)==0:
            print("Nuclide not found")
            self._endf = None
            return
        else:
            if verbose == True:
                print("The nuclide can be found between lines", line[0],line[1])
            self._endf = self._endf_total[line[0]:line[1]]
            self._MAT = []
            self._MF = []
            self._MT = []
            self._pos_section = []
            for el in self._endf:
                self._MAT.append(int(el[66:70]))
                self._MF.append(int(el[70:72]))
                self._MT.append(int(el[72:75]))
                self._pos_section.append(el[60:66])  
            self._MAT = np.array(self._MAT)
            self._MF = np.array(self._MF)
            self._MT = np.array(self._MT)
            self._pos_section=np.array(self._pos_section)
            return

    
    
        return
    
    
    def __convert(self,stringa):
        stringa = stringa[0]+stringa[1:].replace("+","e")
        stringa = stringa[0]+stringa[1:].replace("-","e-")
        return float(stringa)
        
    def get_header(self):
        if self._endf is None:
            print("Select a nuclide first!")
            return
        header_start = np.where(self._MT==451)[0][0]
        header_end = np.where(self._MT==451)[0][-1]
        return self._endf[header_start:header_end-1]
        
    def get_decay(self):
        '''
        Returns a dictionary with the decay data. 
        Works only for B-.
        Decay data is a matrix in the endf/b sub library data format.
        '''
        if self._endf is None:
            print("Select a nuclide first!")
            return
        ii = np.where(self._MT==457)[0]

        decay_raw_data =self._endf[ii[0]:ii[-1]+1]
        chunks, chunk_size = len(decay_raw_data[0][:66]), len(decay_raw_data[0][:66])//6 
        temp = []  #temp is decay data divided in group of 6.
        for el in decay_raw_data:
            temp.append([el[i:i+chunk_size] for i in range(0, chunks, chunk_size) ])

        temp = np.ndarray.flatten(np.array(temp))  #now temp is just flatten with all the information

        self._dic["zaid"] = int(self.__convert(temp[0]))
        self._dic["LIS"] = int(self.__convert(temp[2]))
        self._dic["LISO"] = int(self.__convert(temp[3]))
        self._dic["half-file"] = self.__convert(temp[6])  #seconds
        self._dic["d_half-file"] = self.__convert(temp[7])
        self._dic["decay_data"] = np.zeros((1,5))

        skip = 12+int(self.__convert(temp[10]))
        self._dic["S_p"] = self.__convert(temp[skip])  #spin of the parent 
        self._dic["P_p"] = self.__convert(temp[skip+1]) #parity of the parent
        ndk = int(self.__convert(temp[skip+4]))
        skip = skip +6
        
        for i in range(ndk//6):
            rtype = self.__convert(temp[skip+i*6])
            rfs = temp[skip+1+i*6]
            Q = self.__convert(temp[skip+2+i*6])
            er_Q = self.__convert(temp[skip+3+i*6])
            BR = self.__convert(temp[skip+4+i*6])
            er_BR = self.__convert(temp[skip+5+i*6])
            if rtype == 1:  #save only the info for the beta decay
                self._dic["Q"] = Q*1e-6 #in Mev
                self._dic["d_Q"] = er_Q*1e-6 #in MeV
                self._dic["I_tot"] = BR   #total BR for beta decay
                self._dic["d_I_tot"] = er_BR  
        skip = skip + ndk
        if len(temp[skip:]) == 0:
            print("Warning: No decay information available")
            return self._dic
        while True:
            stype = self.__convert(temp[skip+1])
            lcon = self.__convert(temp[skip +2])
            NER = int(self.__convert(temp[skip+5]))
            skip = skip + 12
            if stype == 0: #skip gamma information
                skip = skip + 6*3*NER
                continue
            elif stype == 1:
                for i in range(NER):
                    jump = 12
                    Emax = self.__convert(temp[skip+i*jump])*1e-6 #in MeV
                    d_Emax = self.__convert(temp[skip+1+i*jump])*1e-6 #in MeV
                    tipo = self.__convert(temp[skip+7+i*jump])
                    I = self.__convert(temp[skip+8+i*jump])#*1e-2 
                    d_I = self.__convert(temp[skip+9+i*jump])#*1e-2
                    if tipo == 1:
                        tipo = "a"
                    else:
                        tipo = str(int(tipo-1))+"u"
                    self._dic["decay_data"] = np.vstack((self._dic["decay_data"],np.array([Emax,d_Emax,I,d_I,tipo])))
                self._dic["decay_data"] = self._dic["decay_data"][1:,:]
                if np.sum(self._dic["decay_data"][:,2].astype(float))>=2: #fix wrong normalization
                    self._dic["decay_data"][:,2] = (self._dic["decay_data"][:,2].astype(float)*1e-2).astype(str)
                    self._dic["decay_data"][:,3] = (self._dic["decay_data"][:,3].astype(float)*1e-2).astype(str)
                break

        return self._dic    
    
class JeffReader:
    
    '''
    Reader for extract cumulative fission yields from JEFF Fission Yields database (Neutron Induced Fission Product Yields).
    
    '''
    
    def __init__(self, path=None):
        
        '''
        some info here
        
        '''
        self._path = path
        self._endf = None
        self._dic = {}
        self._id_row = 8459                        #identifier for cumulative fission yield (cfy) data.
        self._thermal_energy = 2.53e-02            #eV, required to identify thermal and fast cfy
        self._fast_energy = 4e5                    #eV, required to identify thermal and fast cfy
        self._very_fast_energy = 1.4e7             #eV, required to identify thermal and fast cfy
        self._delimiter = [0,11,22,33,44,55,66]
        self._jeff_version = ""
        if self._path is not None:
            self.open_file()
        return
        
    def open_file(self,path=None):
        if path is not None:
            self._path = path
        with open(self._path,"r") as f:
            self._endf = f.readlines()

            
    def get_cfy_from_isotope(self,z,a,isotope_label=None,dic={}):

        zz = zaid(z,a)
        row = []
        for i in range(len(self._endf)):
            try:
                number = float(self._endf[i][:12])
            except:
                try:
                    number = convert(self._endf[i][:12])
                except:
                    continue
            if number == zz:
                row.append(i)
        for r in row:
            if int(self._endf[r][-10:-2]) == self._id_row:
                right_row = r #the row after which are writter the cfy of the chosen isotope
        #get the info inside a dictionary
        i = 1
        info = []
        if isotope_label is None:
            isotope_label = element(zaid(z,a))[-1]
        
        while True:
            line = self._endf[right_row+i]
               
                
            if (((line[:12].strip() == "")) | (convert(line[:12])==0 )) :
                break
            if convert(line[:12]) == self._thermal_energy:
                label = "cumulative_thermal_fy_"+isotope_label
                label_std = "unc_ct_"+isotope_label
                i += 1
                continue
            if convert(line[:12]) == self._fast_energy:
                label = "cumulative_fast_fy_"+isotope_label
                label_std = "unc_cf_"+isotope_label
                i += 1
                continue
            if convert(line[:12]) == self._very_fast_energy:
                break
                #label = "cumulative_ff_fy_"+isotope_label
                #label_std = "unc_cff_"+isotope_label
                #i += 1
                #continue

            for j in range(len(self._delimiter)-1):
                info.append(line[self._delimiter[j]:self._delimiter[j+1]])
                if len(info)==4:
                    for el in info: 
                        if len(el.replace(" ","")) == 0:
                            info = []             
                
                if len(info)==4:               
                    z,a, name = element(convert(info[0]))
                    isomeric = int(convert(info[1]))
                    cfy = convert(info[2])
                    cfy_std = convert(info[3])
                    if isomeric != 0:
                        name += "_"+str(isomeric)+"m" 
                    if name in dic.keys():
                        dic[name].update({label:cfy, label_std:cfy_std, "z":z,"n":int(a-z),"m":isomeric})
                    else:
                        dic[name] = {label: cfy, label_std:cfy_std,"z":z,"n":int(a-z),"m":isomeric}
                    info = []
            i += 1

        return dic
    
    def get_cfy_from_list(self,z_list=[92,92,94,94],a_list=[235,238,239,241],
                          isotope_labels=["235u","238u","239Pu","241Pu"]):
        
        diz = {}
        for i in range(len(z_list)):
            diz =self.get_cfy_from_isotope(z=z_list[i],a=a_list[i],isotope_label=isotope_labels[i],dic=diz)
        return diz    
    
    
    

class EnsdfReader:
    
    '''
    Reader to extract some data (e.g. Q value, half life, normalization for the B- decay) from a ENSDF dataset.
    
    '''
    
    def __init__(self, path=None):
        self._path = path
        self._endsf = None
        self._dic = {}
        
        self._map = {}
        for el in periodictable.elements:
            self._map[el.symbol.upper()] = el.symbol
        if self._path is not None:
            self.open_file()
        return
        
    def open_file(self,path=None):
        if path is not None:
            self._path = path
        with open(self._path,"r") as f:
            self._endsf = f.readlines()
            
    def get_time_conv(self,label=None):
        conv = None
        if label == "MS":
            conv = 1e-3
        if label == "S":
            conv = 1
        if label == "D":
            conv = 24*60*60
        if label == "M":
            conv = 60
        if label == "H":
            conv = 60*60
        if label == "Y":
            conv = 360*24*60*60
        return conv
    
    def _get_half_file(self,string):
        string = string[39:49]
                
        try:
            final_list = []
            for elem in string.split():
                try:
                    final_list.append(float(elem))
                except ValueError:
                    pass
            time = final_list[0]
            #time =float(re.findall(r"[-+]?(?:\d*\.*\d+)", stringa)[0])
            string =re.sub("\.", "", string)
            string ="".join(string.split())
            time_label = re.findall(r"\D+", string)[-1]
            time *=self.get_time_conv(time_label)
        except:
            time = -1
        return time
    
    def _get_Q(self,list_of_string):
        Q = float(list_of_string[-2])
        try:
            Q_std = int(list_of_string[-1])
            x = uncertainties.ufloat_fromstr(str(Q)+"("+str(Q_std)+")")  # Short-hand notation
            Q = x.nominal_value
            Q_std = x.std_dev
        except:
            Q_std = -1
            
        #fix metastable
        try:
            energy_level = float(list_of_string[2])
        except:
            energy_level = 0
        #Q += energy_level
        return Q,Q_std, energy_level
        
    def _get_norm(self,string):
        try:
            norm = float(string[31:39])
        except:
            norm = 0
        try:
            norm_std = float(string[39:41])
        except:
            norm_std = 0
        return norm,norm_std
        
    def get_dict(self):
        label = ""
        i = 1
        dic = {}
        flag = False
        for liner in self._endsf:    
            line = liner.split()
            if len(line)<2:
                continue
            if (line[1] =="P"):
                Q, Q_std, energy_level = self._get_Q(line)
                time = self._get_half_file(string=liner)
                
                name = str(re.findall(r"\d+", line[0])[0])+self._map[re.findall(r"\D+", line[0])[0]]
                name += label
                if i != 1:
                    label = "_"+str(i-1)+"m"  #this should work. The first chunk of data should be the  normal state, then
                                              #the others should follow...
                    name += label
                dic[name] = {"Q": Q, "unc_Q":Q_std,"half_life_sec":time,"energy_level":energy_level}
                if flag == True:
                    dic[name]["BR"] = norm
                    dic[name]["unc_BR"] = norm_std    
                    flag = False
                i+=1
            if (line[1] =="N"):
                norm, norm_std = self._get_norm(string=liner)
                try:
                    name
                    if "BR" in dic[name].keys():
                        flag = True
                    else:
                        dic[name]["BR"] = norm
                        dic[name]["unc_BR"] = norm_std
                except:
                    flag = True

        return dic    
 
class EndfBSubLibraryReader:
    
    '''
    Reader to extract the data from the ENDF/B decay data sub-library (from the NNDC GitLab server).
    The database is not in the ENDF-6 data format.
    
    '''
    
    def __init__(self, path=None):
        
        '''
        some info here
        
        '''
        self._path = path
        
    def set_path(self,path):
        self._path = path

    def get_element(self,z=None,a=None,m=0):
        '''
        Return the beta-decay data of a nuclide given z,a and m.
        The data can be tagged as discreet or continuum.
        If former, data = info for each level of the decay --> end point energy (MeV), uncertainty (MeV), intensity, uncertainty, and type of transition: 'a' for allowed, '1u' for first-forbidden unique, and so on.
        If latter, data = spectrum from a CGM calculation --> energy, electron spectrum, antineutrino spectrum.
        '''
        nfile = self._path
        n_line = 1
        label = str(int(z))+";"+str(int(a))+";"+str(int(m))
        with open(nfile,"r") as f:    
            while True:
                line =f.readline()
                if label in line:
                    line =f.readline()
                    if (("discreet" in line) | ("continuum" in line)):
                        n_data = 0
                        while True:
                            line2 = f.readline()
                            if "----" in line2:
                                break
                            n_data +=1
                        break
                n_line += 1
                if line =="":
                    return None, None, None
    
        data =np.loadtxt(nfile,delimiter=";",dtype="U",skiprows=n_line+1,max_rows=n_data)
        if len(data.shape) == 1:
           data = np.array([data])
        if data.shape[1] == 3:
            tag = "continuum"
        else:
            tag = "discreet"
        if data.shape[0] == 1:
           data = data.flatten()
        return n_line, tag, data
    
    
    
    

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
    stringa = stringa.replace("+","e")
    stringa = stringa.replace("-","e-")
    return float(stringa)


class EndfBReader:
    
    '''
    Simple ENDF/B reader. It should work with JENDL also (they are both in the ENDF-6 data format).    
    
    '''
    
    def __init__(self, path=None):
        

        self._path = path
        self._endf = None
        self._dic = {}
        self._MAT = np.array([])
        self._MF = np.array([])
        self._MT = np.array([])
        self._pos_section = np.array([])
        self._nline = None
        if self._path is not None:
            self.open_file()
        return
        
    def open_file(self,path=None):
        if path is not None:
            self._path = path
        with open(self._path,"r") as f:
            self._endf = f.readlines()
            
        for el in self._endf:
            self._MAT = np.append(self._MAT,int(el[66:70]))
            self._MF = np.append(self._MF, int(el[70:72]))
            self._MT = np.append(self._MT,int(el[72:75]))
            self._pos_section = np.append(self._pos_section,el[60:66])  
        self.__fix_pos()      
        return
        
    def __fix_pos(self):
        i=0
        pos = []
        for el in self._pos_section:
            if (el[-1] != "0") & (el[-4:].find("+") == -1) & (el[-4:].find("-") == -1):
                pos.append(i)
            i += 1
        pos = np.array(pos)
        ii = np.where(self._MT == 457)[0]
        ii =np.intersect1d(pos,ii)
             
        i = 0
        idel = []
        nline = np.array(self._pos_section)[ii]
        for el in nline:
            if el.strip() == "":
                idel.append(i)
            i +=1
        nline = np.delete(nline,idel).astype(int)
        ii = np.delete(ii,idel)
        self._pos_section = ii
        self._nline = nline
        return

    
    def __convert(self,stringa):
        stringa = stringa.replace("+","e")
        stringa = stringa.replace("-","e-")
        return float(stringa)
        
    def get_header(self):
        header_start = np.where(self._MT==451)[0][0]
        header_end = np.where(self._MT==451)[0][-1]
        return self._endf[header_start:header_end-1]
        
    def get_decay(self,loc=4):
        for el in self.get_header():
            pos = el.find("Jpi")
            if pos != -1:
                self._dic["JP_P"] = el[pos+4:pos+12].strip()
            pos = el.find("Parity:")
            if pos != -1:
                self._dic["JP_P"] = el[pos+8:pos+16].strip()        
          
            pos = el.find(" E:")
            if pos != -1:
                self._dic["E_P"] = self.__convert(el[pos+3:el.find("eV")])*1e-3
            pos = el.find("Excitation Energy:")
            if pos != -1:
                self._dic["E_P"] = self.__convert(el[pos+19:pos+45])*1e-3
        
            pos = el.find(" Mode:")
            if pos != -1:
                self._dic["mode"] = el[pos+6:el.find("    ")]  
                
        lineQ = self._pos_section[1]+1
        self._dic["Q"] = self.__convert(self._endf[lineQ][23:33])*1e-3 #keV
        self._dic["d_Q"] = self.__convert(self._endf[lineQ][34:44])*1e-3 #keV
               
        Emax = []
        d_Emax = []
        I = []
        d_I = []
        tipo = []

        lineBeta = self._pos_section[loc]+2
        for i in range(self._nline[loc]):
            line = self._endf[lineBeta+i*2]
            Emax.append(self.__convert(line[1:11]))
            d_Emax.append(self.__convert(line[11:22]))
        
            line = self._endf[lineBeta + i*2+1]
            I.append(self.__convert(line[22:33]))
            d_I.append(self.__convert(line[33:44]))
            tipo.append(self.__convert(line[12:23]))
    
        self._dic["Emax"] = np.array(Emax)*1e-3  #keV
        self._dic["d_Emax"] = np.array(d_Emax)*1e-3  #keV
        self._dic["I"] = np.array(I)
        self._dic["d_I"] = np.array(d_I)
        self._dic["tipo"] = np.array(tipo)
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
            if (line[:12].strip() == ""):
                break
            if float(line[:12]) == self._thermal_energy:
                label = "cumulative_thermal_fy_"+isotope_label
                label_std = "unc_ct_"+isotope_label
                i += 1
                continue
            if float(line[:12]) == self._fast_energy:
                label = "cumulative_fast_fy_"+isotope_label
                label_std = "unc_cf_"+isotope_label
                i += 1
                continue
            if float(line[:12]) == self._very_fast_energy:
                break

            for j in range(len(self._delimiter)-1):
                info.append(line[self._delimiter[j]:self._delimiter[j+1]])
        
                if len(info)==4:
                    z,a, name = element(float(info[0]))
                    isomeric = int(float(info[1]))
                    cfy = float(info[2])
                    cfy_std = float(info[3])
            
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
    Reader to extract some data (e.g. Q value and half life) from a ENSDF dataset.
    
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
        Q += energy_level
        return Q,Q_std
    
    
    def get_dict(self):
        label = ""
        i = 1
        dic = {}
        for liner in self._endsf:    
            line = liner.split()
            if len(line)<2:
                continue
            if (line[1] =="P"):
                
                Q, Q_std = self._get_Q(line)
                time = self._get_half_file(string=liner)
                
                name = str(re.findall(r"\d+", line[0])[0])+self._map[re.findall(r"\D+", line[0])[0]]
                name += label
                dic[name] = {"Q": Q, "unc_Q":Q_std,"half_life_sec":time}
                label = "_"+str(i)+"m"  #this should work. The first chunk of data should be the  normal state, then
                                        #the others should follow...
                i+=1
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
    
    
    
    

'''
 Author: Matteo Borghesi <matteo.borghesi@mib.infn.it>
 Desc  : Collection of utilities used (?) in other codes 
'''

import os, errno
from operator import methodcaller
import numpy as np
import string
import uncertainties

class Config2Dict(object):
    '''
    This class is used to convert in dictionary a config file written in the form:
    
    field1 = value1
    field2 = value2

    '''
    def __init__(self, filename):
        if not os.path.exists(filename):
            raise IOError(errno.ENOENT, os.strerror(errno.ENOENT), filename)
            return
         
        self.__filename=filename
        self.__pars=dict()

        self.__open()
        
        return

    def get_parameters(self):
        return self.__pars

    def __open(self):
        with open(self.__filename) as f:
            # Read config file as a unique string
            text=f.read();
            
            # Remove all the comments and split the string in rows
            nocomments = [s.split('#')[0] for s in text.split('\n')]            
            
            entries = {}
            for el in nocomments:
                pos = el.find("=")
                key = el[:pos].strip()
                value = el[pos+1:].strip()
                entries[key] = value

            self.__pars.clear()

            # Check for float in the entries 
            self.__pars={d: float(entries[d]) if isfloat(entries[d]) else entries[d] for d in entries}

            # Check for list in the entries 
            self.__pars={d: self.__pars[d] if isfloat(self.__pars[d]) or not (',' in self.__pars[d])
                         else self.__pars[d].split(',') for d in self.__pars}

            # Check for float in the lists 
            self.__pars={d: self.__pars[d] if not isinstance(self.__pars[d], list)
                         else [float(l) if isfloat(l) else l for l in  self.__pars[d]]  for d in self.__pars}

            
        return
        

def fix_path(path):
    '''
    fix a str path
    '''
    if path[-1] != "/":
        path += "/"
    return path
    
    
'''
 Check if a string represents float value
''' 
def isfloat(s):
    try: 
        float(s)
        return True
    except ValueError:
        return False    
    
    
def entry(element,lenght):
    out = "{0:<"+str(lenght)+"}"
    out = out.format(str(element))
    out = out[:lenght]
    return out    
    
    
    
def locate_word(array,word):
    '''
    given a numpy array of str, returns all the indexes of the entries which contain 'word' somewhere. Not very fast.
    
    e.g.: locate_word(["hi", "nice3", "is not nice", dude"], word='nice') = [1,2]
    '''
    i = 0
    temp = []
    for el in array:
        if el.find(word) != -1:
            temp.append(i)
        i+=1
    return np.array(temp)
    
def merge_dictionaries(dic1,dic2,keys=None):
    '''
    dic1 and dic2 are dictionaries of dictionary.
    If keys are None, keys are the common keys between dic1 and dic2.
    Else, keys are the key of dic1 which will be updated by dict2[key]. 
    '''
    dic3 = dic1.copy() #for safety
    if keys is None:
        keys = list(set(dic3).intersection(set(dic2)))
    for k in keys:
        dic3[k] = {**dic3[k] , **dic2[k]}
    return dic3
    
def log(message='',level=1):
    out = '|'+level*2*'-'
    if level==-1:
        out = '|'+20*'-'
    else:
        out +=' '+ message
    if level == 0:
        out = message
    print(out)
    return


def findOccurrences(s, ch):
    '''
    Find all the occurences of a character (ch) in a string (s)
    '''
    return [i for i, letter in enumerate(s) if letter == ch]
    
    
class ShorthandFormatter(string.Formatter):
    '''
    from https://pythonhosted.org/uncertainties/user_guide.html
    '''
    def format_field(self, value, format_spec):
        if isinstance(value, uncertainties.UFloat):
            return value.format(format_spec+'S')  # Shorthand option added
        # Special formatting for other types can be added here (floats, etc.)
        else:
            # Usual formatting:
            return super(ShorthandFormatter, self).format_field(
                value, format_spec)
                

def std_to_shorthand(number,std):
    '''
    Given a value as number +/- std, this function returns 
    the shorthand notation as required by the ensdf data format
    value(unc)
    (see ensdf manual)
    '''
    digit=2
    frmtr = ShorthandFormatter()
    x=uncertainties.ufloat(number,std)
    out = frmtr.format("{0:."+str(digit)+"u}", x)
    unc = out[out.find("(")+1:-1]

    loc = unc.find(".")
    if loc != -1:
        if unc[loc+1:] == "0":
            digit = 1
            out = frmtr.format("{0:."+str(digit)+"u}", x)
        else:
            value = out[:out.find("(")]
            value = value[:value.find(".")]
            if len(value)>1:
                pass
            else:
                out = out[:out.find("(")]+"("+unc.replace(".","")+")"
    value = out[:out.find("(")]
    unc = out[out.find("(")+1:-1]
    if len(unc) > 2:
        unc = unc[:unc.find(".")]
        
    return value, unc 


def short_to_std(number,short):
    '''
    Given a value as number(short), this function returns
    the standard notation for uncertanty-> number +/- std
    '''
    if np.isnan(number) | np.isnan(short):
        return None
    try:
        if isinstance(number, str) is False:
            value3 = str(number)
        short3 = "("+str(int(short))+")"

        x = uncertainties.ufloat_fromstr(value3+short3)
    except:
        value2 = str(number)
        digit = str(abs(int(value2[value2.find("e")+1:])))
        temp = "{:."+digit+"f}"
        value4 = temp.format(number)
        short4 = "("+str(int(short))+")"
        x = uncertainties.ufloat_fromstr(value4+short4)
    return float(x.std_dev)


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

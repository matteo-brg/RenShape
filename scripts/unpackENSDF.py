#!/usr/bin/env python
import numpy as np
from base import base_utilities
import argparse
from os.path import exists


def unpack(fname=None,out_dir = None,dtype=".ensdf"):
    
    out_dir=base_utilities.fix_path(out_dir)

    data_file =open(fname,'r') 

    a = data_file.readlines()
    a = np.array(a)

    index = np.where(a=='\n')[0]
    if index[-1] == (len(a)-1):
        index = index[:-1]

    c = 0
    for i in range(len(index)):
        text = a[c:index[i]]
        name = text[0].split()[1]
        name = out_dir + name
        if exists(name+dtype):
            temp = ["\n"]
            temp.extend(text)
        else:
            temp = text
        with open(name+dtype,"a") as file:
            file.writelines(temp)
        c = index[i]+1

    data_file.close()
    return



def main():
    
    usage='unpackENSDF.py -f file_to_unpack -o output_dir'
    parser = argparse.ArgumentParser(description='Unpack a txt file with all the ENSDF data into different files', usage=usage)

    parser.add_argument("-f", "--file"   , dest="fname"   , type=str , help="file to unpack"      , required = True)
    parser.add_argument("-o", "--output"    , dest="save_path"    , type=str , help="output directory"  , required = True)
    
    args = parser.parse_args()    
    
    unpack(args.fname,args.save_path)

if __name__ == "__main__":
    main()

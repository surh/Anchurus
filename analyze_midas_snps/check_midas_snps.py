#!/usr/bin/env python
# Copyright (C) 2017 Sur Herrera Paredes

# coding: utf-8

# In[2]:


import os
import argparse


# In[3]:


def check_sample_dir(path, snp_dir = 'snps/',
                     species_file = 'species.txt',
                     output_dir = 'output/'):
    """Check that there are snps called for species defined by MIDAS"""
    
    snp_dir = "".join([path,"/", snp_dir, "/"])
    species_file = "".join([snp_dir,species_file])
    output_dir = "".join([snp_dir,output_dir,"/"])
    #print(snp_dir)
    #print(species_file)
    #print(output_dir)
    #print("Checking {}".format(snp_dir))
    
    # Check if dir exists
    if not os.path.isdir(snp_dir):
        print("\tNo SNP dir")
        return("No SNP dir")
    
    # Check that species file exists and reads it
    species = []
    if os.path.isfile(species_file):
        with open(species_file) as fh:
            for l in fh:
                species.append(l.rstrip())
        fh.close()
        #print("\tFound {} species".format(str(len(species))))
    else:
        print("\tNo species file")
        return("No species file")
    
    if os.path.isdir(output_dir):
        present = []
        absent = []
        missing = False
        for s in species:
            snp_file = output_dir + s + ".snps.gz"
            if os.path.isfile(snp_file):
                present.append(s)
            else:
                missing = True
                absent.append(s)
        
        if missing:
            print("\tSome species missing")
            return("Some species missing")
        else:
            print("\tComplete")
            return("Complete")
            
    else:
        print("\tNo output dir")
        return("No output dir")
            
            
            
    


# In[4]:


if __name__ == "__main__":
    # Required arguments
    # Does not work in ipython
    # parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    # required = parser.add_argument_group("Required arguments")
    # required.add_argument("--indir", help ="Directory where sample files are located",
    #                       type = str, required = True, default = './')
    # args = parser.parse_args()

    # For iPython
    args = argparse.Namespace
    args.indir = "./"
    
    # Get list of files in dir
    samples = os.listdir(args.indir)
    #samples
    
    # Check ecery dir
    for s in samples:
        sample_path = args.indir + "/" + s
        if os.path.isdir(sample_path):
            print(s,end="")
            res = check_sample_dir(sample_path)
        else:
            print("\t".join([s,"Skipping"]))


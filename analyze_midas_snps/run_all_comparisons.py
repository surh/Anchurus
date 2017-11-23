#!/usr/bin/env python
# Copyright (C) 2017 Sur Herrera Paredes

import pandas as pd
import argparse
import os
import sutilspy


if __name__ == '__main__':
    # comparisons_file = '/home/sur/micropopgen/exp/2017/today4/comparisons.txt'
    # indir = '/godot/hmp/midas/merged.snps/'
    # map_file = '/home/sur/micropopgen/exp/2017/today4/full_map.txt'
    # mk_bin = '/home/sur/micropopgen/src/micropopgen/analyze_midas_snps/MKtest.py'
    # outdir = './out/'
    
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    required = parser.add_argument_group("Required arguments")
    required.add_argument("--comparisons_file", help = "File with comparisons to make",
                          type = str, required = True)
    required.add_argument("--indir", help = "Directory wher output from merge_midas is stored. One subdirectory per species",
                   type = str, required = True)
    required.add_argument('--map_file', help = "Mapping tab-delimited file, must contain and ID and Groups columns",
                      type = str, required = True)
    required.add_argument('--mk_bin', help = 'Executable for MKtest.py',
                    type = str, required = True)
    required.add_argument('--outdir', help = 'Output directory. Will create a subdirectory per species')
    args = parser.parse_args()
        
    # Read list of comparisons
    comparisons = pd.read_csv(args.comparisons_file, sep="\t")
    comparisons.head()
    
    # For every comparison
    for i, r in comparisons.iterrows():
        #print(i)
        #print(r)
        
        species_indir = ''.join([args.indir,'/',
                                 r['Species'],
                                 '/'])
        #print(species_dir)
        
        # Create output directory
        species_outdir = ''.join([args.outdir,
                                  '/',r['Species'],
                                  '/'])
        try:
            os.mkdir(species_outdir)
            print("Creating output directory")
        except FileExistsError:
            print("Directory already exists")
        
        suffix = r['A'] + '_' + r['B']
        suffix = suffix.replace(' ','.')
        #print(suffix)
            
        species_outfile = ''.join([species_outdir,
                                   '/','mk_results.',
                                   suffix,'.txt'])
        species_tables = ''.join([species_outdir,'/',
                                  'mk_tables.',suffix,
                                  '.txt'])
        species_log = ''.join([species_outdir,'/',
                               suffix,'.log'])
        species_err = ''.join([species_outdir,'/',
                               suffix,'.err'])
        
        group1 = ''.join(["'",r['A'],"'"])
        group2 = ''.join(["'",r['B'],"'"])
        
        cmd = [args.mk_bin, '--indir', species_indir,
              '--metadata_file', args.map_file,
               '--group1', group1,
               '--group2', group2,
               '--outfile', species_outfile,
               '--tables', species_tables,
               '1>', species_log,
               '2>', species_err]
        cmd = ' '.join(cmd)
        #print(cmd)
        sutilspy.io.run_command(cmd)
        
        #/micropopgen/src/micropopgen/analyze_midas_snps/MKtest.py --indir /godot/hmp/midas/merged.snps/Porphyromonas_sp_57899/ --metadata_file map.txt --group1 'Supragingival plaque' --group2 'Tongue dorsum' --outfile Porphyromonas_sp_57899_mk_results.txt --tables Porphyromonas_sp_57899_mk_tables.txt 1> log 2> err
    

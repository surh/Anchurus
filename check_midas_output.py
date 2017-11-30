#!/usr/bin/env python
# Copyright (C) 2017 Sur Herrera Paredes

# This script checks whether the output files from midas for
# species and/or snp are present

# Planned features:
# 1.Takes a directory as input, as well as an instruction indicating
# if it should take every sub-directory as a sample or as a single sample
# 2.Takes a list of options to check (species, snps or all [eventually genes])
# 3.Queries database and checks whether results are consistent with records in
# database
# 4.Updates database with results of checks.

import argparse
import os


def check_sample_dir(directory, args):
    """Takes a directory corresponding to a single sample and checks it"""

    # Define whether passed directory was a sample directory or the
    # actual output directory.
    species_dir = ''.join([directory, '/', 'species'])
    snps_dir = ''.join([directory, '/', 'snps'])
    # Not sure if the following makes sense

    if not os.path.isdir(snps_dir):
        snps_dir = './'

    sp_check = snp_check = None
    # Check species
    if 'species' in args.which:
        print("Check species")
        print("\t{}".format(species_dir))
        sp_check, sp_msg = check_species_output(species_dir, args.nspecies)
        print("\tSpecies:{}\t{}".format(sp_check, sp_msg))

    # Check snps
    if 'snps' in args.which:
        snp_check, snp_msg = check_snps_output(snps_dir)

    return sp_check, snp_check


def check_species_output(directory, nspecies):
    """Checks that the MIDAS species output is present and consistent"""
    # If the directory does not exist
    if not os.path.isdir(directory):
        return False, 'speciesdir'

    # If the species profile file does not existt
    sp_profile = ''.join([directory, '/species_profile.txt'])
    if not os.path.isfile(sp_profile):
        return False, 'speciesprof'

    # Check file starts correctly and count lines
    with open(sp_profile) as f:
        # Check header
        header = f.readline()
        header = header.rstrip()
        exp_header = "\t".join(['species_id', 'count_reads', 'coverage',
                                'relative_abundance'])
        if header != exp_header:
            return False, 'profhead'

        # Count lines
        nlines = 0
        for l in f:
            nlines = nlines + 1
        if nlines != args.nspecies:
            return False, 'nspecies'
    f.close()

    return True, 'all'


def check_snps_output(dir):
    """Check that the MIDAS SNPs output is present and consistent"""

    return True, 'all'


def get_sample_dirs(args):
    """Gets list of subdirectories withina directory, and checks that
    there are no non-directory entries"""

    files = os.listdir(args.indir)

    dirs = []
    not_dirs = []
    for f in files:
        path = "".join([args.indir, "/", f])
        if os.path.isdir(path):
            dirs.append(path)
        else:
            not_dirs.append(path)

    if args.notdirs == 'fail':
        try:
            if len(not_dirs) > 0:
                raise ValueError
        except:
            print(("ERROR: The passed directory ({}) contains non-directory "
                   "entries").format(args.indir))
            raise

    return dirs


def process_arguments():
    """Reads command line arguments, and performs checks and processes
    if any of the arguments require it"""

    # Read arguments
    parser_format = argparse.ArgumentDefaultsHelpFormatter
    parser = argparse.ArgumentParser(formatter_class=parser_format)
    required = parser.add_argument_group("Required arguments")
    required.add_argument("--indir", help=("Path to directory where MIDAS "
                                           "output is located"),
                          required=True)
    required.add_argument("--type", help=("Indicates whether the directory "
                                          "is from a single sample or if "
                                          "it contains a number of sample "
                                          "sub-directories"),
                          choices=['single', 'multi'],
                          required=True)
    parser.add_argument("--which", help=("Indicates which MIDAS output to "
                                         "check"),
                        choices=['species', 'snps', 'all'],
                        default='all')
    parser.add_argument("--notdirs", help=("Indicates what to do in case "
                                           "that a 'multi' directory is "
                                           "passed, and it contains some "
                                           "non directory entries"),
                        choices=['ignore', 'fail'],
                        default='ignore')
    parser.add_argument("--nspecies", help=("Number of species in database"),
                        type=int,
                        default=5952)

    # Still needs output options
    print("Reading arguments")
    args = parser.parse_args()

    # Check and process arguments
    # Check that a directory is passed
    try:
        if not os.path.isdir(args.indir):
            raise ValueError
    except:
        print("ERROR:Indir ({}) is not a directory".format(args.indir))
        raise

    # Prepare list of directories
    if args.type == 'multi':
        args.dirs = get_sample_dirs(args)
    elif args.type == 'single':
        args.dirs = [args.indir]
    else:
        print("ERROR: Incorrect type of directory passed")
        raise ValueError

    # Define which comparisons to make
    if args.which == 'all':
        args.which = ['species', 'snps']

    return args


if __name__ == "__main__":
    # Arguments
    args = process_arguments()
    # print(args.dirs)

    # Check directories
    for d in args.dirs:
        print("Processing directoy:\n\t{}".format(d))
        checks = check_sample_dir(d, args)

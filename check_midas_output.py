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
    # if not os.path.isdir(species_dir):
    #     species_dir = './'
    # if not os.path.isdir(snps_dir):
    #     snps_dir = './'

    # Define which comparisons to make
    if args.which == 'all':
        args.which = ['species', 'snps']

    sp_check = snp_check = None
    # Check species
    if 'species' in args.which:
        sp_check = check_species_output(species_dir)

    # Check snps
    if 'snps' in args.which:
        snp_check = check_snps_output(snps_dir)

    return sp_check, snp_check


def check_species_output(dir):
    """Checks that the MIDAS species output is present and consistent"""
    check = True
    return check


def check_snps_output(dir):
    """Check that the MIDAS SNPs output is present and consistent"""
    check = True
    return check


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

    return args


if __name__ == "__main__":
    # Arguments
    args = process_arguments()

    # Prepare list of directories
    if os.type == 'multi':
        dirs = get_sample_dirs(args)
    elif os.type == 'single':
        dirs = [args.indir]
    else:
        print("ERROR: Incorrect type of directory passed")
        raise ValueError

    # Check directories
    for d in dirs:
        checks = check_sample_dir(d, args)

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


def get_sample_dirs(args):
    """Gets list of subdirectories withina directory, and checks that
    there are no non-directory entries"""

    files = os.listdir(args.indir)

    dirs = []
    not_dirs = []
    for f in files:
        path = "".join([args.indir, "/", f])
        if os.path.isdir(path):
            dirs.append(f)
        else:
            not_dirs.append(f)
    if args.notdirs == 'fail':
        try:
            if len(not_dirs) > 0:
                raise ValueError
        except:
            print(("ERROR: The passed directory ({}) contains non-directory "
                   "entries").format(args.indir))
            raise

    return dirs


if __name__ == "__main__":
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

    # Check that a directory is passed
    try:
        if not os.path.isdir(args.indir):
            raise ValueError
    except:
        print("ERROR:Indir ({}) is not a directory".format(args.indir))

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
        check_sample_dir(d,args)

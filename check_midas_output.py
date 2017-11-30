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
import gzip


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
    sp_msg = snp_msg = None
    # Check species
    if 'species' in args.which:
        print("Check species")
        print("\t{}".format(species_dir))
        sp_check, sp_msg = check_species_output(species_dir, args.nspecies)
        print("\tSpecies:{}\t{}".format(sp_check, sp_msg))

    # Check snps
    if 'snps' in args.which:
        print("Check SNPs")
        print("\t{}".format(snps_dir))
        snp_check, snp_msg = check_snps_output(snps_dir)
        print("\tSNPs:{}\t{}".format(snp_check, snp_msg))

    return [sp_check, sp_msg, snp_check, snp_msg]


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
        nlines = count_remaining_lines(f)
        if nlines != args.nspecies:
            return False, 'nspecies'
    f.close()

    return True, 'all'


def check_snps_output(directory):
    """Check that the MIDAS SNPs output is present and consistent"""
    # If the directory does not exist
    if not os.path.isdir(directory):
        return False, 'snpdir'

    # Check and read list of species with SNPs
    sp_file = ''.join([directory, '/', 'species.txt'])
    if not os.path.isfile(sp_file):
        return False, 'specieslist'
    else:
        sp_list = []
        with open(sp_file) as f:
            for s in f:
                sp_list.append(s.rstrip())
        f.close()
        nspecies = len(sp_list)
        # print(nspecies)
        # print(sp_list)

    # Check summary file
    summary_file = ''.join([directory, '/', 'summary.txt'])
    if not os.path.isfile(summary_file):
        return False, 'summaryfile'
    else:
        with open(summary_file) as f:
            # Check header in summary
            header = f.readline()
            header = header.rstrip()
            exp_header = "\t".join(['species_id', 'genome_length',
                                    'covered_bases', 'fraction_covered',
                                    'mean_coverage', 'aligned_reads',
                                    'mapped_reads'])
            if header != exp_header:
                return False, 'summaryhead'

            # Make sure number of species in summary matches number
            # in species list
            nlines = count_remaining_lines(f)
            if nlines != nspecies:
                return False, 'summarynum'
        f.close()

    # Check output dir
    out_dir = ''.join([directory, '/output'])
    if not os.path.isdir(out_dir):
        return False, 'outdir'
    else:
        # Check number of files matches
        file_list = os.listdir(out_dir)
        # print(len(file_list))
        if len(file_list) != nspecies:
            return False, 'outnum'

        # Check output for every species
        for s in sp_list:
            s_file = ''.join([out_dir, '/', s, '.snps.gz'])
            # print(s_file)
            if not os.path.isfile(s_file):
                return False, 'speciesout'
            else:
                if not check_species_output_file(s_file):
                    return False, 'outsnps'

    return True, 'all'


def check_species_output_file(f):
    """Check species output file from MIDAS SNPs"""
    # print(f)
    with gzip.open(f, 'r') as fh:
        header = fh.readline()
        header = header.decode()
        header = header.rstrip()
        # print(header)
        # print("Hola")
        exp_header = "\t".join(['ref_id', 'ref_pos', 'ref_allele',
                                'depth', 'count_a', 'count_c',
                                'count_g', 'count_t'])
        check = exp_header == header
        # print(check)
    fh.close()

    return check


def count_remaining_lines(fh):
    """Counts remaining lines in a file handle"""
    # Count lines
    nlines = 0
    for l in fh:
        nlines = nlines + 1

    return nlines


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


def write_output_file(check, args):
    with open(args.outfile, 'w') as fh:
        for res in check:
            sample = os.path.basename(res[1])
            print(sample)
            line = "\t".join(res.extend())
            line = ''.join([line, "\n"])
            fh.write(line)
    fh.close()

    return 0


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
    parser.add_argument("--outfile", help=("Name of the filw to write "
                                           "output"),
                        type=str, default='midas_check_results.txt')

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
    checks = []
    for d in args.dirs:
        print("Processing directoy:\n\t{}".format(d))
        res = [d]
        res.extend(check_sample_dir(d, args))
        checks.append(res)
    #print(checks)

    write_output_file(checks, args)

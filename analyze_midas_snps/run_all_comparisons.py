#!/usr/bin/env python
# Copyright (C) 2017 Sur Herrera Paredes
# This file is based on the MarkerScanner.pl script of the AMPHORA
# pipeline by Martin Wu.

# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.

# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.

import pandas as pd
import argparse
import os
import sutilspy
import fyrd


def process_arguments():
    # Read arguments
    parser_format = argparse.ArgumentDefaultsHelpFormatter
    parser = argparse.ArgumentParser(formatter_class=parser_format)
    required = parser.add_argument_group("Required arguments")

    # Define description
    parser.description = ("Takes a directory with midas outputs and "
                          "a mapping file for samples and submits "
                          "all comparisons to MKtest.py")

    # Define required arguments
    required.add_argument("--comparisons_file", help=("File with comparisons "
                                                      "to make"),
                          type=str, required=True)
    required.add_argument("--indir", help=("Directory wher output from "
                                           "merge_midas is stored. One "
                                           "subdirectory per species"),
                          type=str, required=True)
    required.add_argument('--map_file', help=("Mapping tab-delimited file, "
                                              "must contain and ID and Groups "
                                              "columns"),
                          type=str, required=True)
    required.add_argument('--mk_bin', help=('Executable for MKtest.py'),
                          type=str, required=True)
    required.add_argument('--outdir', help=('Output directory. Will create '
                                            'a subdirectory per species'),
                          type=str, required=True)

    # Define other arguments
    parser.add_argument("--permute", help=('Number of permutations to use '
                                           'to calculate p-value. If 0 '
                                           'no permutations will be done.'),
                        type=int, default=0)
    parser.add_argument("--mode", help=('Whether to use fyrd to '
                                        'parallelize per comparison or '
                                        'submit serially'),
                        type=str, default='fyrd', choices=['fyrd', 'bash'])
    parser.add_argument("--nrows", help=('Number of rows to read. '
                                         'Useful to limit time when testing'),
                        type=float, default=float('inf'))
    parser.add_argument("--maxjobs", help=('Maximum number of fyrd jobs'),
                        type=int, default=100)
    parser.add_argument("--wait", help=('Wait for all fyrd jobs to finish'),
                        action='store_true')

    # Read arguments
    print("Reading arguments")
    args = parser.parse_args()

    # Processing goes here if needed

    return args


if __name__ == '__main__':
    args = process_arguments()

    # Read list of comparisons
    comparisons = pd.read_csv(args.comparisons_file, sep="\t")
    # comparisons.head()

    # Create overall outptu directory
    print("Creating output directory")
    try:
        os.mkdir(args.outdir)
    except FileExistsError:
        raise("Directory already exists")

    print("========ITERATIONG OVER COMPARISONS========")
    # For every comparison
    for i, r in comparisons.iterrows():
        print("Species:{} in {} vs {}".format(r['Species'], r['A'], r['B']))
        species_indir = ''.join([args.indir, '/',
                                 r['Species'],
                                 '/'])

        # Create output directory
        print("Creating species output directory")
        species_outdir = ''.join([args.outdir,
                                  '/', r['Species'],
                                  '/'])
        try:
            os.mkdir(species_outdir)
        except FileExistsError:
            print("Directory already exists")

        print("Creating MKtest.py command")
        suffix = r['A'] + '_' + r['B']
        suffix = suffix.replace(' ', '.')

        species_outfile = ''.join([species_outdir,
                                   '/', 'mk_results.',
                                   suffix, '.txt'])
        species_tables = ''.join([species_outdir, '/',
                                  'mk_tables.', suffix,
                                  '.txt'])
        species_log = ''.join([species_outdir, '/',
                               suffix, '.log'])
        species_err = ''.join([species_outdir, '/',
                               suffix, '.err'])

        group1 = ''.join(["'", r['A'], "'"])
        group2 = ''.join(["'", r['B'], "'"])

        cmd = [args.mk_bin, '--indir', species_indir,
               '--metadata_file', args.map_file,
               '--group1', group1,
               '--group2', group2,
               '--outfile', species_outfile,
               '--tables', species_tables]

        # If not equal to default, pass it
        if args.nrows != float('inf'):
            cmd.extend(['--nrows', str(args.nrows)])

        # Append output logging
        cmd.extend(['1>', species_log,
                    '2>', species_err])

        cmd = ' '.join(cmd)
        # print(cmd)

        if args.mode == 'bash':
            sutilspy.io.run_command(cmd)
        elif args.mode == 'fyrd':
            job_name = 'mktest' + r['species']
            print("Creating fyrd job {}".format(job_name))
            job = fyrd.Job(cmd,
                           clean_files=False,
                           clean_outputs=False,
                           nodes=1, cores=1,
                           time='00:20:00',
                           mem='1000mb',
                           partition='',
                           name=job_name,
                           runpath=os.getcwd(),
                           outpath='logs',
                           scriptpath='scripts')
            print("\tSubmitting job")
            job.submit(max_jobs=args.maxjobs)
        else:
            raise ValueError("Unreconized mode")

    print("========ITERATIONG OVER COMPARISONS========")

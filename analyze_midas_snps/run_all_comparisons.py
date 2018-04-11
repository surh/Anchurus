#!/usr/bin/env python
# Copyright (C) 2017-2018 Sur Herrera Paredes

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
# import sutilspy
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
    parser.add_argument("--logs", help='Directory for fyrd logs',
                        type=str, default='logs/')
    parser.add_argument("--scripts", help='Directory for fyrd scripts',
                        type=str, default='scripts')
    parser.add_argument("--pseudocount", help=("Pseudocount value to use "
                                               "in contingency tables"),
                        default=0, type=int)
    parser.add_argument("--permutations", help=("Number of permutations to "
                                                "perform to establish "
                                                "significance"),
                        type=int, default=0)
    parser.add_argument("--seed", help="Permutation seed",
                        type=int, default=None)
    parser.add_argument("--test", help=("Eventually specify test to perform."
                                        "all performs all tests. G performs "
                                        "a G test without correction."
                                        "G_Yates performs a G test with the "
                                        "Yates correction. G_Williams "
                                        "performs a G test with the Williams "
                                        "correction. hg performs the "
                                        "hypergeometric (Fisher's Exact) "
                                        "test. NI rerturns the neutrality "
                                        "index. alpha returs the Eyre-"
                                        "Walker alpha. Ratio returns the MK "
                                        "rati. DoS is the direction of "
                                        "selection statistic."),
                        default="hg", type=str,
                        choices=['all', 'G', 'G_Yates', 'G_Williamps',
                                 'hg', 'NI', 'alpha', 'ratio', 'DoS'])

    # Read arguments
    print("Reading arguments")
    args = parser.parse_args()

    # Processing goes here if needed
    if args.seed is None and args.permutations > 0:
        args.seed = np.random.randint(1000)*2 + 1
        print("Seed is ", args.seed)

    return args


def make_output_directories(args):
    """Prepare the required directories"""

    # Create overall outptut directory
    print("Creating output directory")
    try:
        os.mkdir(args.outdir)
    except FileExistsError:
        raise("Directory already exists")

    # Make fyrd directories
    if args.mode == 'fyrd':
            if os.path.isdir(args.logs):
                raise FileExistsError("Directory for fyrd logs ({}) "
                                      "already exists".format([args.logs]))
            else:
                os.mkdir(args.logs)
            if os.path.isdir(args.scripts):
                raise FileExistsError("Directory for fyrd scripts ({}) "
                                      "already exists".format([args.scripts]))
            else:
                os.mkdir(args.scripts)

    return


if __name__ == '__main__':
    args = process_arguments()

    # Read list of comparisons
    comparisons = pd.read_csv(args.comparisons_file, sep="\t")
    # comparisons.head()

    # Make output directories
    make_output_directories(args)

    print("========ITERATING OVER COMPARISONS========")
    # For every comparison
    jobs = []
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
               '--tables', species_tables,
               '--permutations', str(args.permutations),
               '--test', args.test,
               '--seed', str(args.seed),
               '--pseudocount', str(args.pseudocount)]

        # If not equal to default, pass it
        if args.nrows != float('inf'):
            cmd.extend(['--nrows', str(args.nrows)])

        # Append output logging
        cmd.extend(['1>', species_log,
                    '2>', species_err])

        cmd = ' '.join(cmd)
        # print(cmd)

        if args.mode == 'bash':
            # sutilspy.io.run_command(cmd)
            print("Running command:\n>{}".format(cmd))
            os.system(cmd)
        elif args.mode == 'fyrd':
            job_name = 'mktest' + r['Species']
            print("Creating fyrd job {}".format(job_name))
            job = fyrd.Job(cmd,
                           clean_files=False,
                           clean_outputs=False,
                           nodes=1, cores=1,
                           time='04:00:00',
                           mem='10000mb',
                           partition='',
                           name=job_name,
                           runpath=os.getcwd(),
                           outpath=args.logs,
                           scriptpath=args.scripts)
            print("\tSubmitting job")
            job.submit(max_jobs=args.maxjobs)
            job.submit(max_jobs=args.maxjobs)
            jobs.append(job)
        else:
            raise ValueError("Unreconized mode")

    print("========DONE ITERATING OVER COMPARISONS========")

    if args.mode == 'fyrd' and args.wait:
        print("========WAITING FOR MKTEST JOBS TO COMPLETE========")
        for j in jobs:
            j.wait()

        print("========DONE WAITING FOR MKTEST JOBS TO COMPLETE========")

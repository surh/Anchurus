#!/usr/bin/env python
# Copyright (C) 2018 Sur Herrera Paredes

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

import argparse
import os


def process_arguments():
    # Read arguments
    parser_format = argparse.ArgumentDefaultsHelpFormatter
    parser = argparse.ArgumentParser(formatter_class=parser_format)
    required = parser.add_argument_group("Required arguments")

    # Define description
    parser.description = ("Take output from midas_merge.py snps and "
                          "extract data from one gene.")

    # Define required arguments
    required.add_argument("--indir", help=("Midas merge directory output"),
                          required=True, type=str)
    required.add_argument("--gene", help=("Gene to extract or file with "
                                          "gene_ids"),
                          type=str, required=True)

    # Define other arguments
    parser.add_argument("--type", help=("Type of gene parameter. id or file"),
                        type=str,
                        default="id",
                        choices=['id', 'file'])
    parser.add_argument("--outdir", help=("Output directory to write files"))

    # Read arguments
    print("Reading arguments")
    args = parser.parse_args()

    # Processing goes here if needed
    # Check if file exists
    if args.type == 'file':
        if not os.path.isfile(args.gene):
            raise FileNotFoundError()

    return args


def check_midas_merge_output(dir):
    """Check if files from midas merge snps exists"""

    depth_file = '/'.join([dir, 'snps_depth.txt'])
    freq_file = '/'.join([dir, 'snps_freq.txt'])
    info_file = '/'.join([dir, 'snps_info.txt'])

    # Check if files exists
    if not os.path.isfile(depth_file):
        raise FileNotFoundError("snps_depth.txt does not exist.")
    if not os.path.isfile(freq_file):
        raise FileNotFoundError("snps_freq.txt does not exist.")
    if not os.path.isfile(info_file):
        raise FileNotFoundError("snps_info.txt does not exist.")

    return depth_file, freq_file, info_file


def gene_ids(gene, type):
    """Get list of gene ids to keep"""

    if type == "id":
        genes = [gene]
    elif type == "file":
        with open(gene, 'r') as ih:
            genes = ih.readlines().rstrip('\n')
        ih.close()
    else:
        raise ValueError("Unknown type passsed")

    return genes


if __name__ == "__main__":
    args = process_arguments()

    depth_file, freq_file, info_file = check_midas_merge_output(args.indir)
    genes = gene_ids(args.gene, args.type)

    print(genes)

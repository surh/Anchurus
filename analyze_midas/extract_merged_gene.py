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
    parser.add_argument("--outdir", help=("Output directory to write files"),
                        type=str,
                        default='out/')

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
            genes = ih.readlines()
            genes = [g.rstrip("\r\n") for g in genes]
        ih.close()
    else:
        raise ValueError("Unknown type passsed")

    return genes


def extract_genes_info(gene, depth_file, freq_file, info_file, outdir):
    """Extract info from positions for some genes"""

    # First read info file
    o_info_file = '/'.join([outdir, 'snps_info.txt'])
    keep_positions = []
    with open(info_file, 'r') as i_info, open(o_info_file, 'w') as o_info:
        header = i_info.readline().rstrip("\r\n")
        o_info.write(''.join([header, "\n"]))
        for line in i_info:
            line = line.rstrip("\r\n")
            fields = line.split("\t")
            # The following are hard-coded field numbers
            # Rely on MIDAS 1.3.1
            curr_gene = fields[12]
            pos_id = fields[0]
            if curr_gene in genes:
                keep_positions.append(pos_id)
                o_info.write(''.join([line, "\n"]))
    i_info.close()
    o_info.close()


if __name__ == "__main__":
    args = process_arguments()

    # Check input
    depth_file, freq_file, info_file = check_midas_merge_output(args.indir)

    # Get gene names
    genes = gene_ids(args.gene, args.type)

    print(genes)
    # Prepare output directory
    if os.path.isdir(args.outdir):
        raise ValueError("Output directory already exists")
    else:
        os.mkdir(args.outdir)

    # Getting data
    extract_genes_info(genes, depth_file, freq_file, info_file, args.outdir)

#!/usr/bin/env python
# Copyright (C) 2021 Sur Herrera Paredes

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
import numpy as np
import argparse


def process_arguments():
    # Read arguments
    parser_format = argparse.ArgumentDefaultsHelpFormatter
    parser = argparse.ArgumentParser(formatter_class=parser_format)
    required = parser.add_argument_group("Required arguments")

    # Define description
    parser.description = ("Read SNP files from UHGG and produce core SNP aln")

    # Define required arguments
    required.add_argument("--snps", help=("UHGG snp file"),
                          required=True, type=str)
    required.add_argument("--meta", help=("UHGG metadata"),
                          required=True, type=str)
    required.add_argument("--spec", help=("Species"),
                          required=True, type=str)

    # Define other arguments
    parser.add_argument("--output", help=("Core SNPs aln file"),
                        type=str,
                        default="core.aln.fasta")
    parser.add_argument("--presence_thres", help=("Presence threshold"),
                        type=float,
                        default=0.8)
    parser.add_argument("--min_genomes", help=("Minimum number of genomes"),
                        type=int,
                        default=5)

    # Read arguments
    print("Reading arguments")
    args = parser.parse_args()

    # Processing goes here if needed

    return args


if __name__ == "__main__":
    args = process_arguments()

    # Select genomes
    meta = pd.read_csv(args.meta, sep="\t",
                       dtype={'Genome_accession': str})
    meta = meta[meta.MGnify_accession == args.spec]
    meta = meta.loc[meta.Country.isin(["United States", "China"]), ]
    selected_genomes = meta.Genome

    # Check if enough genomes
    if not (meta.Country.value_counts() >= args.min_genomes).all():
        print("Not enough genomes")
        quit()

    with open(args.snps, 'r') as snps:
        # Read and process header
        header = snps.readline()
        header = header.split("\t")
        genome_list = np.array(header[4:])
        ii = np.intersect1d(genome_list, selected_genomes, return_indices=True)
        ii = ii[1]
        genome_list = genome_list[ii]

        # Dat = pd.DataFrame(index = genome_list)
        Dat = {g: [] for g in genome_list}

        count = 0
        while True:
            count += 1
            if (count % 10000) == 0:
                print("Analyzing SNP {}".format(count))
            # if count > 9470:
            #      break

            # Get next line from file
            line = snps.readline()
            # if line is empty
            # end of file is reached
            if not line:
                break

            line = line.split("\t")
            geno = np.array(line[4:])
            geno = geno[ii]
            geno = np.array(geno, dtype=np.float64)
            geno[geno == 255] = np.nan

            # Skip if not enough genomes
            iinan = ~np.isnan(geno)
            propgenomes = np.sum(iinan) / geno.size
            # print(propgenomes)
            if propgenomes < args.presence_thres:
                # print("no")
                continue

            # print("yes")

            # Substitute for sequence
            alleles = np.array([line[2], line[3]])
            seq = np.repeat("-", geno.shape)
            seq[geno == 0 & iinan] = alleles[0]
            seq[geno == 1 & iinan] = alleles[1]

            # Append
            # snpid = line[0] + '.' + line[1]
            # Dat = Dat.append(pd.DataFrame(seq, index = genome_list, columns = [snpid]))
            # Dat = Dat.merge(pd.DataFrame(seq, index = genome_list, columns = [snpid]),
            #                 "left", left_index = True, right_index = True)
            # Dat[snpid] = pd.DataFrame(seq, index = genome_list)
            # Dat[snpid] = seq
            for i in range(genome_list.size):
                Dat[genome_list[i]].append(seq[i])
    snps.close()

    filename = args.output
    with open(filename, 'w') as out:
        for g in Dat.keys():
            out.write('>' + g + "\n")
            out.write(''.join(Dat[g]) + "\n")
    out.close()

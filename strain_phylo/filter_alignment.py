#!/usr/bin/env python
# Copyright (C) 2019 Sur Herrera Paredes

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
import numpy as np
from Bio import AlignIO
from Bio.Alphabet import single_letter_alphabet
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Align import MultipleSeqAlignment


def process_arguments():
    # Read arguments
    parser_format = argparse.ArgumentDefaultsHelpFormatter
    parser = argparse.ArgumentParser(formatter_class=parser_format)
    required = parser.add_argument_group("Required arguments")

    # Define description
    parser.description = ("Filter alns")

    # Define required arguments
    required.add_argument("--input", help=("Input file"),
                          required=True, type=str)

    # Define other arguments
    parser.add_argument("--output", help=("Output dir"),
                        type=str,
                        default="concatenated.aln.fasta")

    # Read arguments
    print("Reading arguments")
    args = parser.parse_args()

    # Processing goes here if needed

    return args


def align2array(aln):
    """Convert multiple sequence alignment object to numpy array.
    Taken from tutorial."""

    a_array = np.array([list(rec) for rec in aln], np.character)

    return(a_array)


def array2align(arr, names, alphabet):
    """Convert numpy array to multiple sequence alignment.
    Adapted from documentation"""

    records = []

    # Iterate over array rows (i.e. records)
    for i in range(arr.shape[0]):
        seq = ''.join(np.array(arr[i], dtype=str))
        name = names[i]

        # Concatenate sequence records
        records.append(SeqRecord(Seq(seq, alphabet), id=name))

    # Convert to MSA
    new_aln = MultipleSeqAlignment(records)

    return(new_aln)


def filter_alignment_file(infile, outfile, gap_prop=0.99,
                          remove_singletons=True,
                          alphabet=single_letter_alphabet,
                          input_format='fasta',
                          output_format='fasta'):
    """Take file with a single alignment, filter it and
    write a new file"""

    print("\treading")
    aln = AlignIO.read(infile, input_format)
    print("\tfiltering")
    filtered = filter_alignment(aln=aln, gap_prop=gap_prop,
                                remove_singletons=remove_singletons,
                                alphabet=alphabet)
    print("\twriting")
    AlignIO.write(filtered, outfile, output_format)

    return(outfile)


def filter_alignment(aln, gap_prop=0.99, remove_singletons=True,
                     alphabet=single_letter_alphabet):
    """Function to filter a numpy array that represents an alignment,
    where rows are records and columns are positions. Assumes gaps are
    given by '-'"""

    # Get sequence records
    nseqs = len(aln)
    rec_names = [r.id for r in aln]

    # Convert to numpu array
    a_array = align2array(aln)

    # Prepare index. Positions to be removed will be changed
    index = np.ones(a_array.shape[1], dtype=int)

    # Iterate over columns
    # print("Iterating")
    for i in range(a_array.shape[1]):
        c = a_array[:, i]
        # print("\t", c)
        counts = np.unique(c, return_counts=True)

        # Remove constant columns
        if counts[0].shape == (1,):
            index[i] = 0
            continue

        # Count gaps
        ngaps = counts[1][b'-' == counts[0]]
        if ngaps.shape[0] == 0:
            # print("hello")
            ngaps = np.zeros(1, dtype=int)

        if ngaps / nseqs > gap_prop:
            index[i] = 0
            continue

        if remove_singletons:
            if counts[1].size == 2 and counts[1].min() == 1:
                index[i] = 0
                continue

    # Use index to slice array
    index = np.array(index, dtype=bool)
    # print(index.sum())
    filtered = a_array[:, index]

    # Convert back to alignent
    new_aln = array2align(arr=filtered, names=rec_names, alphabet=alphabet)

    return new_aln


if __name__ == "__main__":
    args = process_arguments()

    filter_alignment_file(args.input, args.output,
                          gap_prop=1, remove_singletons=True)

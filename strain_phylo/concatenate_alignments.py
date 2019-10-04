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
import os
# import align_markers as am
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
    parser.description = ("Concatenate alns")

    # Define required arguments
    required.add_argument("--indir", help=("Input dir"),
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


def concatenate_alignments(alns, alphabet=single_letter_alphabet, gap='-'):
    """Take a list of multiple sequence alignments and
    concatenate them, fill with gaps where missing sequences."""

    # Get list of species from alignments
    species = []
    species_per_aln = []
    for a in alns:
        specs = [r.id for r in a]
        species.extend(specs)
        species_per_aln.append(specs)

    species = list(set(species))

    # Create empty alignmet
    new_aln = MultipleSeqAlignment([SeqRecord(Seq('', alphabet),
                                              id=s) for s in species])

    # Iterate over each species, re-ordering when neccessary
    for i in range(len(alns)):
        # print("alginment", i)
        specs = species_per_aln[i]
        if specs != species:
            # print("\treordering")
            # new_alns.append(reorder_alignment(aln=alns[i],
            #                                   specs=specs, species=species))
            new_aln = new_aln + reorder_alignment(aln=alns[i], specs=specs,
                                                  species=species,
                                                  alphabet=alphabet,
                                                  gap=gap)
        else:
            # print("matched")
            new_aln = new_aln + alns[i]

    return new_aln


def reorder_alignment(aln, specs, species,
                      alphabet=single_letter_alphabet, gap='-'):
    """Take an alignment and reorder it acording to species list.
    Add records as gapped if missing"""

    new_aln = []
    missing_seq = ''.join([gap] * aln.get_alignment_length())
    for s in species:
        # Check if species exist in alignment
        try:
            i = specs.index(s)
        except ValueError:
            i = -1
        except:
            raise

        if i >= 0:
            new_aln.append(aln[i])
        elif i == -1:
            new_aln.append(SeqRecord(Seq(missing_seq, alphabet), id=s))

    new_aln = MultipleSeqAlignment(new_aln)

    return(new_aln)


if __name__ == "__main__":
    args = process_arguments()

    aln_files = os.listdir(args.indir)
    Alns = []
    for f in aln_files:
        a = AlignIO.read(os.path.join(args.indir, f), 'fasta')
        Alns.append(a)

    aln = concatenate_alignments(alns=Alns)
    AlignIO.write(aln, args.output, 'fasta')

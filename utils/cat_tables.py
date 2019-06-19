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

import pandas as pd
import argparse


def process_arguments():
    # Read arguments
    parser_format = argparse.ArgumentDefaultsHelpFormatter
    parser = argparse.ArgumentParser(formatter_class=parser_format)
    required = parser.add_argument_group("Required arguments")

    # Define description
    parser.description = ("Takes a list of files containing tab-delimited "
                          "files with headers, and concatenate them into "
                          "a single file.")

    # Define required arguments
    required.add_argument("FILES", help=("Files to concatenate. Each file "
                                         "must be tab-delimited, and the "
                                         "first row must contain the headers"),
                          type=str, nargs='+')
    required.add_argument("--outfile", help=("Filename to write the "
                                             "concatenated results"),
                          required=True, type=str)

    # Read arguments
    print("Reading arguments")
    args = parser.parse_args()

    # Processing goes here if needed

    return args


if __name__ == "__main__":
    args = process_arguments()

    res = pd.DataFrame()
    print("Concatenating files:")
    i = 0
    for file in args.FILES:
        print("\tProcessing {}".format(file))
        tab = pd.read_csv(file, sep="\t", index_col=False, dtype=str)
        res = res.append(tab)
        i = i + 1
    print("Processed {} files".format(str(i)))
    print("Writting {}".format(args.outfile))
    res.to_csv(args.outfile, sep="\t", na_rep="NA",
               header=True, index=False)

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


def process_arguments():
    # Read arguments
    parser_format = argparse.ArgumentDefaultsHelpFormatter
    parser = argparse.ArgumentParser(formatter_class=parser_format)
    required = parser.add_argument_group("Required arguments")

    # Define description
    parser.description = ("This script takes a list of tab delmited files "
                          "and concatenates them, while keeping only one "
                          "header and making sure the headers are consistent "
                          "accross files.")

    # Define required arguments
    required.add_argument("--infiles", help=("File paths of tab delimited "
                                             "files to concatenate."),
                          required=True, type=str, nargs='*')

    # Define other arguments
    parser.add_argument("--outfile", help=("Name of the outfile to create."),
                        type=str,
                        default="concatenated.txt")

    # Read arguments
    print("Reading arguments")
    args = parser.parse_args()

    # Processing goes here if needed

    return args


if __name__ == "__main__":
    args = process_arguments()

    with open(args.outfile, 'w') as oh:
        header = ''
        for file in args.infiles:
            with open(file) as tab:

                # Check header`
                new_header = tab.readline()
                if header == '':
                    header = new_header
                    oh.write(header)
                else:
                    if new_header != header:
                        raise ValueError("Header does not match for "
                                         "file {}".format(file))

                for line in tab:
                    oh.write(line)
            tab.close()
    oh.close()

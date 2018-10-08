# Anchurus
Analyzing SNPs in metagenomes

Anchurus was the son of king Midas of ancient Greece. The repository is named afer him
because the metagenomic SNPs are obtained through [MIDAS](https://github.com/snayfach/MIDAS)

# Directories

* **analyze_midas**: This directory contains code for analysis after midas has been used to
call SNPs.
* **mktest**: This directory contains code to perform the McDonald-Kreitman test on
gene frequencies, between two groups.
* **run_midas**: This directory contains wrappers to run midas on a series of samples. It
requires that there are two files per sample corresponding to Illumina forward and reverse
reads. Wrappers are nextflow pipelines.
* **vmwa**: This directory contains code for testing association between genetic variants
of bacteria and a host phenotype.


# Copyright & license

    (C) Copyright 2018 Sur Herrera Paredes

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.

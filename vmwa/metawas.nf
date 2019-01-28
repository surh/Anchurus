#!/usr/bin/env nextflow
// Copyright (C) 2018 Sur Herrera Paredes

// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <http://www.gnu.org/licenses/>.

// Nextflow pipeline that takes the output from midas merge of SNPs and
// and mapping file grouping samples into two groups, and performs imputation
// of mean variants with bimban, and association via LMM of each SNP with
// the grouping factor


// Params
params.genomes_file = 'genomes.txt'
params.midas_dir = ''
params.outdir = 'output'
params.map_file = 'map.txt'

// Process params
genomes_file = file(params.genomes_file)
map_file = Cahnnel.fromPath(params.map_file)

// Read list of genomes
reader = genomes_file.newReader()
GENOMES = []
while(str = reader.readLine()){
  GENOMES = GENOMES + [str]
}


process metawas{
  cpus 1
  time '2h'
  module 'R/3.5.1server'
  maxForks 4
  publishDir params.outdir, mode: 'copy', saveAs: {"${genome}_lmm.assoc.txt"}

  input:
  val genome from Channel.from(GENOMES)
  file map_file

  output:
  file "metawas/lmm/lmm.assoc.txt"

  """
  Rscript /home/sur/micropopgen/src/HMVAR/inst/scripts/bugwas.r \
    ${params.midas_dir} \
    ${genome}
  """

}

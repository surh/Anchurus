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
params.dirs = ''
params.bin = ''
params.focal_group = ''
params.outdir = 'metawas'
params.map = 'map.txt'
params.pcs = ''
params.impute = false
params.time = '5:00:00'
params.njobs = 20

// Process params
if(params.dirs == ''){
  error "No list of directories provided."
}
if(params.bin == ''){
  error "genome_metawas.r exectuable not provided."
}
if(params.focal_group == ''){
  error "A focal group must be provided."
}

dirs = file(params.dirs)
map = file(params.map)
if(params.pcs != ''){
  pcs = file(params.pcs)
}else{
  pcs = false
}



// Read list of genomes
// reader = genomes_file.newReader()
// GENOMES = []
// while(str = reader.readLine()){
//   GENOMES = GENOMES + [str]
// }
// Read list of dirs
DIRS = Channel.fromPath(dirs).
  splitCsv(sep: "\t").
  map{row -> tuple(row[0], file(row[1]))}


process metawas{
  label 'r'
  cpus 1
  time params.time
  // module 'R/3.5.1server'
  maxForks params.njobs
  publishDir params.outdir, mode: 'copy', saveAs: {"${genome}_lmm.results.txt"}
  // errorStrategy 'ignore'

  input:
  set genome, file(specdir) from Channel.from(DIRS)
  file map_file
  file pcs

  output:
  // file "metawas/lmm.results.txt"

  script:
  """
  Rscript ${params.bin} \
    ${params.midas_dir} \
    ${genome}
  """

}

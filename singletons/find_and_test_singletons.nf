#!/usr/bin/env nextflow
// Copyright (C) 2020 Sur Herrera Paredes

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

// Parameters
params.indir = ''
params.outdir = 'output'
params.map = 'map.txt'

// Process Params
map = file(params.map)
indir = file(params.indir)


// Get files
INFOS = Channel.fromPath("${indir}/**/snps_info.txt")
  .map{infofile -> tuple(infofile.getParent().name,
    file(infofile))}
ALLELES = Channel.fromPath("${indir}/**/snps_alleles.txt")
  .map{infofile -> tuple(infofile.getParent().name,
    file(infofile))}

process mktest{
  label 'r'
  tag "$spec"
  publishDir "$params.outdir/singletons", mode: 'rellink',
    pattern: "output/singletons.txt",
    saveAs: {"${spec}.txt"}
  publishDir "$params.outdir/singleton_tests", mode: 'rellink',
    pattern: "output/singleton_tests.txt",
    saveAs: {"${spec}.txt"}

  input:
  tuple spec, file(alleles), file(info) from ALLELES.join(INFOS)
  file map

  output:
  file "output/singletons.txt"
  file "output/singleton_tests.txt"

  """
  Rscript ${workflow.projectDir}/mktest.r \
    $alleles \
    $info \
    $map \
    --outdir output
  """
}

// Example nextflow.config
/*
process{
  partition = 'hbfraser,hns,owners'
  maxForks = 200
  stageInMode = 'rellink'
  errorStrategy = 'finish'
  withLabel: 'r'{
    module = 'R/3.6.1'
    memory = '4G'
    time  = '4h'
  }
}
executor{
  name = 'slurm'
  submitRateLimit = '1 sec'
  queueSize = 500
}
*/

#!/usr/bin/env nextflow
// Copyright (C) 2019 Sur Herrera Paredes

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

// Params
params.midas_dir = ''
params.maps_dir = ''
params.outdir = 'output'
params.min_snvs = 5000
params.depth_thres = 1

// Create channels
maps_dir = file(params.maps_dir)
midas_dir = file(params.midas_dir)
MAPS = Channel.fromPath("${maps_dir}/*")
  .map{map_file -> tuple(map_file.name.replaceAll('\\map\\.txt$', ""),
    file(map_file))}
MIDAS = Channel.fromPath("${midas_dir}/*", type: 'dir')
  .map{snv_dir -> tuple(snv_dir.name, file(snv_dir))}

process snv_cor{
  label 'r'
  tag "$spec"
  publishDir "${params.outdir}/${spec}",
    mode: 'rellink'

  input:
  tuple val(spec), file("map.txt"), file(snv_dir) from MAPS.join(MIDAS)
  val min_snvs from params.min_snvs
  val depth_thres from params.depth_thres

  output:
  file "cors_all.txt.gz"
  file "cors_synonymous.txt.gz"
  file "cors_nonsynonymous.txt.gz"

  """
  stitch_file.r ${workflow.projectDir}/snv_cors.r \
    $snv_dir \
    map.txt \
    $depth_thres \
    $min_snvs \
    ./bigcor/
    1
  """
}

// Example nextflow.config (for sherlock)
/*
process{
  maxForks = 300
  stageInMode = 'rellink'
  errorStrategy = 'finish'
  queue = 'hbfraser,hns'
  withLabel: 'r'{
    module = 'R/3.5.1'
    memory = '15G'
    time  = '48h'
  }
}
executor{
  name = 'slurm'
  submitRateLimit = '1 sec'
  queueSize = 800
}
*/

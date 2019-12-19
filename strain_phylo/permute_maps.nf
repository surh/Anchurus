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

// params
params.indir = ''
params.outdir = ''
params.nperms = 100
params.seed = 739

// Create input Channel
map_dir = file(params.indir)
MAPS = Channel.fromPath("${map_dir}/*")
    .map{map_file -> tuple(map_file.name.replaceAll('\\.map\\.txt$', ''),
    file(map_file))}

process permute_maps{
    label 'r'
    tag "$spec"
    publishDir params.outdir,
        pattern: "output",
        saveAs: {{"$spec"},
        mode: 'rellink'

    input:
    tuple val(spec), fle("map.txt") from inputs
    val nperms from params.nperms
    val seed from params.seed

    output:
    tuple $spec, file("output/") into PERMAPS

    """
    ${workflow.projectDir}/permute_maps.r map.txt \
        --outdir output \
        --nperms $nperms \
        --seed $seed
    """
}

// Example nextflow.config
/*
process{
  maxForks = 100
  stageInMode = 'rellink'
  errorStrategy = 'finish'
  withLabel: 'r'{
    module = 'R/3.6.1'
    memory = '1G'
    time  = '1h'
  }
}
executor{
  name = 'slurm'
  submitRateLimit = '1 sec'
  queueSize = 500
}
*/

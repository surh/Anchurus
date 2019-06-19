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

// Script that runfs MKtest.py on a set of folders.

// Script to run mktest on all genomes in a folder

// Paramteres
params.genomesdir = './'
params.map = 'map.txt'
params.group1 = ''
params.group2 = ''
params.nperm = 0
params.seed = 371
params.queue = 'hbfraser,hns'
params.njobs = 300
params.outdir = 'results/'
params.min_count = 1

// Process Params
map = file(params.map)


// Get files
GENOMES = Channel.fromPath("${params.genomesdir}/*", type: 'dir')
FREQS = Channel.fromPath("${params.genomesdir}/**/snps_freq.txt")
DEPTHS = Channel.fromPath("${params.genomesdir}/**/snps_depth.txt")
INFOS = Channel.fromPath("${params.genomesdir}/**/snps_info.txt")

process genome_mktest{
  publishDir path: params.outdir, pattern: "*_mktest.txt", mode: 'copy'
  cpus 1
  maxForks params.njobs
  queue params.queue
  errorStrategy 'retry'
  module 'fraserconda'
  time 3000.m
  memory 5.GB

  input:
  file genome from GENOMES
  file freqs from FREQS
  file depths from DEPTHS
  file infos from INFOS
  file map

  output:
  file "${genome}_mktest.txt" into MKTESTS

  """
  ${workflow.projectDir}/MKtest.py \
    --indir $genome \
    --metadata_file $map \
    --group1 ${params.group1} \
    --group2 ${params.group2} \
    --functions classes \
    --min_count ${params.min_count} \
    --permutations ${params.nperm} \
    --seed ${params.seed} \
    --test ratio \
    --min_cov 0 \
    --outfile ${genome}_mktest.txt
  """
}

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

// Pipeline for mktest via R script on list of genomes and genome DIRS
// with output from MIDAS merge.

// Paramteres
params.genomes = ''
params.bin = ''
params.map = 'map.txt'
params.focal_group = false
params.queue = 'hbfraser,hns'
params.njobs = 300
params.outdir = 'results/'
params.min_count = 1
params.freq_thres = 0.5

// Process Params
map = file(params.map)
genomes = file(params.genomes)

// Get files
GENOMES = Channel.fromPath(genomes).
  splitCsv(sep: "\t").
  map{row -> tuple(row[0], file(row[1]))}

process genome_mktest{
  label 'r'
  publishDir path: params.outdir,
    pattern: "mktest.txt",
    mode: 'copy',
    saveAs: {"${genome}_mktest.txt"}
  cpus 1
  maxForks params.njobs
  queue params.queue

  input:
  set genome, file(genomedir) from GENOMES
  file "$genomedir/snps_info.txt"
  file "$genomedir/snps_depth.txt"
  file "$genomedir/snps_freq.txt"
  file map

  output:
  file "mktest.txt" into MKTESTS

  script:
  if(params.focal_group)
    """
    ${params.bin} \
      $genomedir \
      $map \
      --depth_thres ${params.min_count} \
      --frep_thres ${params.freq_thres} \
      --focal_group ${params.focal_group} \
      --test \
      --outfile mktest.txt
    """
  else{
    """
    ${params.bin} \
      $genome \
      $map \
      --depth_thres ${params.min_count} \
      --frep_thres ${params.freq_thres} \
      --test \
      --outfile ${genome}_mktest.txt
    """
  }
}

// Example nextflow.config
/*
process {
  executor = 'slurm'
  withLabel: r {
    module = 'R/3.5.1server'
    time = '24h'
    memory = '5G'
  }
}
*/

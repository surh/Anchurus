#!/usr/bin/env nextflow
// Copyright (C) 2017 Sur Herrera Paredes

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

// Nextflow pipeline that submits sample fastq files to midas to obtain
// species profiles


// Main parameters
params.samples = 'samples.txt'
params.indir = 'samples/'
params.outdir = 'midas/'
params.sample_col = 1
params.logdir = 'logs'
params.queue = 'hbfraser,owners,bigmem,hns,normal'
params.memory = '10G'
params.time = '4:00:00'
params.cpus = 8
params.njobs = 200

// Process params
samples = file(params.samples)
sample_col = params.sample_col - 1

// Read samples file
reader = samples.newReader()
SAMPLES = []
while(str = reader.readLine()){
  // Extract sample and run IDs
  sample = str.split("\t")[sample_col]
  SAMPLES = SAMPLES + [tuple(sample,
    file("${params.indir}/${sample}_read1.fastq.bz2"),
    file("${params.indir}/${sample}_read2.fastq.bz2"))]
}


// Call run_midas.py species on every sample
process midas_species{
  cpus params.cpus
  time params.time
  memory params.memory
  maxForks params.njobs
  module 'MIDAS/1.3.1'
  queue params.queue
  publishDir params.outdir, mode: 'copy'

  input:
  set sample, f_file, r_file from SAMPLES

  output:
  set sample,
    file("${sample}/species/log.txt"),
    file("${sample}/species/readme.txt"),
    file("${sample}/species/species_profile.txt") into OUTPUTS

  """
  run_midas.py species ${sample} \
    -1 ${f_file} \
    -2 ${r_file} \
    -t ${params.cpus} \
    --remove_temp
  """
}

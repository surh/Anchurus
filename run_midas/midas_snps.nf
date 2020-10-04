#!/usr/bin/env nextflow
// Copyright (C) 2018-2020 Sur Herrera Paredes

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
// snp profiles

// Main parameters
params.samples = 'samples.txt'
params.indir = 'samples/'
params.outdir = 'midas/'
params.sample_col = 1
// params.queue = 'hbfraser,bigmem,hns,owners'
// params.memory = '10G'
// params.time = '4:00:00'
// params.cpus = 8
// params.njobs = 200
params.species_cov = 3.0
params.mapid = 94.0
params.mapq = 20
params.baseq = 30
params.readq = 30
params.aln_cov = 0.75
params.trim = 0
params.discard = false
params.baq = false
params.adjust_mq = false
// Steps argument needs to be implemented

// Process params
samples = file(params.samples)
sample_col = params.sample_col - 1
if( params.trim > 0 ) {
  trim = "--trim ${params.trim}"
}
else {
  trim = ''
}
if( params.discard == true ) {
    discard = '--discard'
}
else {
    discard = ''
}
if( params.baq == true ) {
  baq = '--baq'
}
else {
  baq = ''
}
if( params.adjust_mq == true ) {
  adjust_mq = 'adjust_mq'
}
else {
  adjust_mq = ''
}


// Read samples file
reader = samples.newReader()
SAMPLES = []
while(str = reader.readLine()){
  // Extract sample and run IDs
  sample = str.split("\t")[sample_col]
  SAMPLES = SAMPLES + [tuple(sample,
    file("${params.indir}/${sample}_read1.fastq.bz2"),
    file("${params.indir}/${sample}_read2.fastq.bz2"),
    file("${params.outdir}/${sample}/species/species_profile.txt"))]
}



// Call run_midas.py species on every sample
process midas_species{
  label 'midas'
  // cpus params.cpus
  // time params.time
  // memory params.memory
  // maxForks params.njobs
  // module 'MIDAS/1.3.1'
  // queue params.queue
  publishDir params.outdir, mode: 'copy'
  // errorStrategy 'retry'
  // maxRetries 2

  input:
  set sample, f_file, r_file, spec_profile from SAMPLES

  output:
  set sample,
    file("${sample}/snps/log.txt"),
    file("${sample}/snps/readme.txt"),
    file("${sample}/snps/species.txt"),
    file("${sample}/snps/summary.txt"),
    file("${sample}/snps/output/*.snps.gz") into OUTPUTS

  """
  mkdir ${sample}
  mkdir ${sample}/species
  cp ${spec_profile} ${sample}/species/
  run_midas.py snps ${sample} \
    -1 ${f_file} \
    -2 ${r_file} \
    -t ${params.cpus} \
    --remove_temp \
    --species_cov ${params.species_cov} \
    --mapid ${params.mapid} \
    --mapq ${params.mapq} \
    --baseq ${params.baseq} \
    --readq ${params.readq} \
    --aln_cov ${params.aln_cov} \
    -m local \
    ${trim} \
    ${discard} \
    ${baq} \
    ${adjust_mq}
  """
}


/* Example nextflow.config
process{
  queue = 'hbfraser,hns,owners,bigmem'
  maxForks = 500
  errorStrategy = 'finish'
  stageInMode = 'rellink'
  time = '5h'
  memory = '5G'
  withLabel: 'midas'{
    cpus = 4
    module="MIDAS/1.3.1"
    time = { 10.s * task.attempt }
    memory = {15.G * task.attempt }
    maxRetries = 5
    errorStrategy = { task.attempt < 5 ? 'retry' : 'finish'}
  }
}

executor{
  name = 'slurm'
  queueSize = 500
  submitRateLitmit = '30 min'
}
*/

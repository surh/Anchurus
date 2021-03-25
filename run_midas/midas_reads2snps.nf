#!/usr/bin/env nextflow
// Copyright (C) 2018-2021 Sur Herrera Paredes

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
// species profiles and snps

// Main parameters
params.indir = 'samples/'
params.outdir = 'midas/'
params.db = 'midas_db'
params.cpus = 4

// SNPs params
params.species_cov = 3.0
params.mapid = 94.0
params.mapq = 20
params.baseq = 30
params.readq = 20
params.aln_cov = 0.75
params.trim = 0
params.discard = false
params.baq = false
params.adjust_mq = false

// Process params
indir = file(params.indir)
midas_db = file(params.db)
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

// Use file pairs to create list of files
Channel
  .fromFilePairs("$indir/*_{1,2}.fq.gz")
  .into{READS1; READS2}

process midas_species{
  tag "$sample"
  label "midas"
  cpus params.cpus
  publishDir params.outdir, mode: 'rellink'

  input:
  set sample, file(reads) from READS1
  file midas_db from midas_db

  output:
  set sample,
    file("${sample}/species/log.txt"),
    file("${sample}/species/readme.txt") into SPECLOGS
  set sample, file("${sample}/species/species_profile.txt") into SPECPROFS

  """
  run_midas.py species ${sample} \
    -1 ${reads[0]} \
    -2 ${reads[1]} \
    -t ${params.cpus} \
    -d $midas_db
  """
}

process midas_snps{
  label 'midas'
  tag "$sample"
  publishDir params.outdir, mode: 'rellink'
  cpus params.cpus

  input:
  set sample, file(reads), file(spec_profile) from READS2.join(SPECPROFS)
  file midas_db from midas_db

  output:
  set sample,
    file("${sample}/snps/log.txt"),
    file("${sample}/snps/readme.txt"),
    file("${sample}/snps/species.txt"),
    file("${sample}/snps/summary.txt"),
    file("${sample}/snps/output/*.snps.gz"),
    file("${sample}/snps/temp/genomes*") into OUTPUTS

  """
  mkdir ${sample}
  mkdir ${sample}/species
  cp ${spec_profile} ${sample}/species/
  run_midas.py snps ${sample} \
    -1 ${reads[0]} \
    -2 ${reads[1]} \
    -t ${params.cpus} \
    --species_cov ${params.species_cov} \
    --mapid ${params.mapid} \
    --mapq ${params.mapq} \
    --baseq ${params.baseq} \
    --readq ${params.readq} \
    --aln_cov ${params.aln_cov} \
    -m global \
    -d $midas_db \
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
    module="MIDAS/1.3.1"
    time = { task.attempt < 3 ? '12h' : '48h' }
    memory = { 15.GB * task.attempt }
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

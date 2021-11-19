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
// snp profiles

// Main parameters
// params.samples = 'samples.txt'
params.indir = 'samples/'
params.specdir = "species/"
params.outdir = 'midas/'
params.db = 'midas_db'
// params.sample_col = 1
params.cpus = 4
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
params.paired = false
// Steps argument needs to be implemented

// Process params
// samples = file(params.samples)
midas_db = file(params.db)
// sample_col = params.sample_col - 1
indir = file(params.indir)
specdir = file(params.specdir)
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


// // Read samples file
// reader = samples.newReader()
// SAMPLES = []
// while(str = reader.readLine()){
//   // Extract sample and run IDs
//   sample = str.split("\t")[sample_col]
//   SAMPLES = SAMPLES + [tuple(sample,
//     file("${params.indir}/${sample}_read1.fastq.bz2"),
//     file("${params.indir}/${sample}_read2.fastq.bz2"),
//     file("${params.outdir}/${sample}/species/species_profile.txt"))]
// }
// Get read file
READSP = Channel.empty()
READSU = Channel.empty()
if(params.paired){
  READSP = Channel
    .fromFilePairs("$indir/*_{1,2}.fq.gz")
}else{
  READSU = Channel
    .fromPath("$indir/*.fastq.gz")
    .map{ infile -> tuple(infile.name.replaceAll(/\.fastq\.gz/, ''),
      file(infile)) }
}


// Get specfiles
Channel.fromPath("$specdir/*/species/species_profile.txt")
  .map{specfile -> tuple(specfile.getParent().getParent().name,
    file(specfile))}
  .into{SPECPROFSP; SPECPROFSU}

// Call run_midas.py species on every sample
process midas_snps_paired{
  label 'midas'
  tag "$sample"
  publishDir params.outdir, mode: 'rellink'
  cpus params.cpus

  input:
  // set sample, f_file, r_file, spec_profile from SAMPLES
  set sample, file(reads), file(spec_profile) from READSP.join(SPECPROFSP)
  file midas_db from midas_db

  output:
  set sample,
    file("${sample}/snps/log.txt"),
    file("${sample}/snps/readme.txt"),
    file("${sample}/snps/species.txt"),
    file("${sample}/snps/summary.txt"),
    file("${sample}/snps/output/*.snps.gz"),
    file("${sample}/snps/temp/genomes*") into OUTPUTSP

  when:
  params.paired

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

process midas_snps_unpaired{
  label 'midas'
  tag "$sample"
  publishDir params.outdir, mode: 'rellink'
  cpus params.cpus

  input:
  // set sample, f_file, r_file, spec_profile from SAMPLES
  set sample, file(reads), file(spec_profile) from READSU.join(SPECPROFSU)
  file midas_db from midas_db

  output:
  set sample,
    file("${sample}/snps/log.txt"),
    file("${sample}/snps/readme.txt"),
    file("${sample}/snps/species.txt"),
    file("${sample}/snps/summary.txt"),
    file("${sample}/snps/output/*.snps.gz"),
    file("${sample}/snps/temp/genomes*") into OUTPUTSU

  when:
  !params.paired

  """
  mkdir ${sample}
  mkdir ${sample}/species
  cp ${spec_profile} ${sample}/species/
  run_midas.py snps ${sample} \
    -1 ${reads[0]} \
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

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

// Nextflow pipeline that takes a matrix of SNP frequencies, a file
// with phenotypes and a matrix of covariates and performs association
// for all SNPs and testing via permutation.

// Params
params.snps = 'snps.txt'
params.phenotype = 'phenotype.txt'
params.covariates = 'covariates.txt'
params.nperm = 1000
params.bindir = '/home/sur/micropopgen/src/Anchurus/vmwa'
params.seed = 3094229
params.nrows = 20000
params.queue = 'hbfraser,hns,owners'
params.njobs = 300

// Hard-coded params
time_baseline = 1440.m
time_extra = 1440.m

// Process
if(params.nsamples <= 0){
  throw new Exception("Invalid number of samples (nsamples).")
}
snps = file(params.snps)
phenotype = file(params.phenotype)
covariates = file(params.covariates)

process permute_samples{
  cpus 1
  maxForks params.njobs
  queue params.queue
  errorStrategy 'retry'
  maxRetries 2
  module 'R'
  time 20.m
  memory 500.MB

  input:
  file phenotype
  file covariates

  output:
  file "phenoperm*.txt" into PHENOPERMS

  """
  ${params.bindir}/permute_samples.r \
    $phenotype \
    $covariates \
    --nperm ${params.nperm} \
    --prefix phenoperm
  """
}

process run_vmwa{
  cpus 1
  maxForks params.njobs
  queue params.queue
  errorStrategy 'retry'
  maxRetries params.vmwa_retries
  module 'R'
  time {time_baseline + ((task.attempt - 1 ) * time_extra)}
  memory { 4.GB + ((task.attempt - 1) * 2.GB) }

  input:
  file phenotype
  file snps
  file covariates

  output:
  file 'results.txt' into MWAS

  """
  ${params.bindir}/chunked_vmwa.r \
    $snps \
    $covariates \
    $phenotype \
    ${params.nsamples} \
    --chunk_size ${params.nrows} \
    --outfile results.txt
    --plot
  """
}

process run_vmwa_perms{
  cpus 1
  maxForks params.njobs
  queue params.queue
  errorStrategy 'retry'
  maxRetries params.vmwa_retries
  module 'R'
  time {time_baseline + ((task.attempt - 1 ) * time_extra)}
  memory { 4.GB + ((task.attempt - 1) * 2.GB) }

  input:
  file phenoperm from PHENOPERMS.flatten()
  file snps
  file covariates

  output:
  file 'results.txt' into MWASPERM

  """
  ${params.bindir}/chunked_vmwa.r \
    $snps \
    $covariates \
    $phenoperm \
    ${params.nsamples} \
    --chunk_size ${params.nrows} \
    --outfile results.txt
  """
}

process combine_permutations{
  cpus 1
  maxForks params.njobs
  queue params.queue
  errorStrategy 'retry'
  maxRetries 2
  module 'R'
  time {240.m + ((task.attempt - 1 ) * 60.m)}
  memory { 4.GB + ((task.attempt - 1) * 2.GB) }

  input:
  file 'original.txt' from MWAS
  file '*.perm' from MWASPERM.collect()

  output:
  file 'p.value_histogram.svg'
  file 'results.txt'

  """
  ${params.bindir}/process_permutations.r \
     original.txt \
    --perms *.perm
  """
}

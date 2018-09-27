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
params.nrows = 5000

// Process
snps = file(params.snps)
phenotype = file(params.phenotype)
covariates = file(params.covariates)

process permute_samples{

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

  input:
  file phenoperm from PHENOPERMS.flatten()
  file snps
  file covariates

  output:
  file 'results.txt' into MWAS

  """
  ${params.bindir}/vMWAS.r \
    $snps \
    $covariates \
    $phenoperm \
    --nrows ${params.nrows} \
    --outdir ./
  """
}

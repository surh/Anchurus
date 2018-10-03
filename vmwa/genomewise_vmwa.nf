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

// Nextflow pipeline that runs script that calcultates vMWA.r on each
// directory corresponding to a genome.

params.genomesdir = './'
params.phenotype = 'phenotype.txt'
params.covariates = 'covariates.txt'
params.nperm = 1000
params.bindir = '/home/sur/micropopgen/src/Anchurus/vmwa'
params.seed = 5743
params.queue = 'hbfraser,hns,owners'
params.njobs = 300
params.maf = 0.05
params.outdir = 'results/'

// Process Params
phenotype = file(params.phenotype)
covariates = file(params.covariates)

// Get files
GENOMES = Channel.fromPath("${params.genomesdir}/*", type: 'dir')
SNPS = Channel.fromPath("${params.genomesdir}/**/snps_freq.txt")


process genome_vmwa{
  publishDir path: params.results, pattern: "*_associations.txt", mode: 'copy'

  input:
  file genome from GENOMES
  file snps from SNPS
  file phenotype
  file covariates

  output:
  file "${genome}_associations.txt" into ASSOC
  file 'p.value_histogram.svg' into PVALSHIST
  file 'p.value_qqplot.pnf' into PVALSQQ
  file 'P_histogram.svg' optional true into PHIST
  file 'P_qqplot.pnf' optional true into PQQ

  """
  ${params.bindir}/vMWAS.r \
    $snps \
    ${params.covariates} \
    ${params.phenotype} \
    --outfile ${genome}_associations.txt \
    --maf ${params.maf} \
    --permutations ${params.nperm} \
    --lib ${params.bindir} \
    --seed ${params.seed} \
    --plot
  """
}

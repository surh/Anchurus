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

// Nextflow pipeline to produce phylogeny of the same species from
// different metagenomic samples.

// parameters
params.indir = ""
params.genomes_dir = ""
params.map_file = "map.txt"
params.min_cov = 0.8
params.outdir = "output/"


map = file(params.map_file)
genomes_dir = file(params.genomes_dir)
INDIRS = Channel.fromPath("${params.indir}/*", type: 'dir')

process alns_from_metagenomes{
  label 'r'
  publishDir "${params.outdir}/gene_coverages/",
    pattern: "output/*gene_coverage.txt",
    mode: 'rellink'
  publishDir "${params.outdir}/gene_alns/",
    pattern: "output/*aln.fasta",
    mode: 'rellink'

  input:
  file midas_dir from INDIRS
  file genomes_dir
  file map

  output:
  file "output/*.gene_coverage.txt" into COVS
  file "output/*.aln.fasta" optional true into SINGLE_ALNS

  """
  ${workflow.projectDir}/alns_from_metagenomes.r \
    $midas_dir \
    $genomes_dir \
    --min_cov ${params.min_cov} \
    --map_file $map \
    --outdir output/
  """
}

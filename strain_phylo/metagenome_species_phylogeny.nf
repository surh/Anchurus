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
params.fasttree_threads = 8


map = file(params.map_file)
genomes_dir = file(params.genomes_dir)
INDIRS = Channel.fromPath("${params.indir}/*", type: 'dir')

process alns_from_metagenomes{
  label 'r'
  publishDir "${params.outdir}/gene_coverages/",
    pattern: "output/${spec}.gene_coverage.txt",
    saveAs: {"${spec}.gene_coverage.txt"},
    mode: 'rellink'
  publishDir "${params.outdir}/gene_alns/",
    pattern: "output/*aln.fasta",
    mode: 'rellink'

  input:
  file midas_dir from INDIRS
  file genomes_dir
  file map

  output:
  file "output/${spec}.gene_coverage.txt" into COVS
  set spec, file("output/*aln.fasta") into CORE_ALNS

  script:
  spec = midas_dir.fileName
  """
  ${workflow.projectDir}/alns_from_metagenomes.r \
    $midas_dir \
    $genomes_dir \
    --min_cov ${params.min_cov} \
    --map_file $map \
    --outdir output/
  """
}

process concatenate_alns{
  // label 'py3'

  input:
  set spec, file("test/") from CORE_ALNS

  """
  ls -la
  """
}

// process concatenate_alignments{
//   label 'py3'
//   publishDir "${params.outdir}/cat_alns/", mode: 'rellink'
//
//   input:
//   set spec, file("alns/") from CORE_ALNS
//
//   output:
//   set spec, file("${spec}.concatenated.aln.fasta") into CAT_ALNS
//
//   """
//   ${workflow.projectDir}/concatenate_alignments.py \
//     --indir alns \
//     --output ${spec}.concatenated.aln.fasta
//   """
// }

// process filter_alignment{
//   label 'py3'
//   publishDir "${params.outdir}/filtered_alns/", mode: 'rellink'
//
//   input:
//   set spec, file(cat_aln) from CAT_ALNS
//
//   output:
//   set spec, file("${spec}.concatenated_filtered.aln.fasta") into FILTERED_ALNS
//
//   """
//   ${workflow.projectDir}/filter_alignment.py \
//     --input $cat_aln \
//     --output ${spec}.concatenated_filtered.aln.fasta
//   """
// }
//
// process fasttree{
//   label 'fasttree'
//   cpus params.fasttree_threads
//   publishDir "${params.outdir}/trees/", mode: 'rellink'
//
//   input:
//   set spec, file(aln) from FILTERED_ALNS
//
//   output:
//   set spec, file("${spec}.tre") into TREES
//
//
//   """
//   export OMP_NUM_THREADS=${params.fasttree_threads}
//   FastTree \
//     -nt $aln \
//     > ${spec}.tre
//   """
// }


// Example nextflow config
/*
process{
  maxForks = 100
  withLabel: 'r'{
    module = 'R/3.5.1server'
  }
  withLabel: 'py3'{
    conda = '/opt/modules/pkgs/anaconda/3.6/envs/fraserconda'
  }
  withLabel: 'fasttree'{
    module = 'FastTree/2.1.10'
  }
}
executor{
  name = 'slurm'
  submitRateLimit = '1 sec'
  queueSize = 500
}
*/

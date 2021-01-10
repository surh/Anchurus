#!/usr/bin/env nextflow
// Copyright (C) 2021 Sur Herrera Paredes

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

// Parameters
params.indir = ''
params.outdir = 'output'
params.fasttree_threads = 4

indir = file(params.indir)

// Get files
INFOS = Channel.fromPath("${indir}/**/snps_info.txt")
  .map{infofile -> tuple(infofile.getParent().name,
    file(infofile))}
ALLELES = Channel.fromPath("${indir}/**/snps_alleles.txt")
  .map{infofile -> tuple(infofile.getParent().name,
    file(infofile))}
COREGENES = Channel.fromPath("${indir}/**/core_genes.tsv")
  .map{corefile -> tuple(corefile.getParent().name,
    file(corefile))}

process core_alns{
  label 'r'
  tag "$spec"
  publishDir "$params.outdir/core_alns", mode: 'rellink',
    saveAs: {"${spec}.aln.fasta"}

  input:
  tuple spec, file(alleles), file(info), file(coregenes)
    from ALLELES.join(INFOS).join(COREGENES)

  output:
  tuple spec, file("core_aln.fasta") into COREALNS

  """
  Rscript ${workflow.projectDir}/get_core_gene_alns.r \
    $info \
    $alleles \
    $coregenes \
    core_aln.fasta
  """
}

process fasttree{
  label 'FastTree'
  tag "$spec"
  cpus params.fasttree_threads

  input:
  tuple spec, file(corealn) from COREALNS
  val threads from params.fasttree_threads

  output:
  file "core.tre"

  """
  export OMP_NUM_THREADS=$threads
  FastTree \
    -nt -gtr < $corealn > core.tre
  """
}


// Example nextflow.config
/*
process{
  partition = 'hbfraser,hns,owners'
  maxForks = 200
  stageInMode = 'rellink'
  errorStrategy = 'finish'
  withLabel: 'r'{
    module = 'R/3.6.1'
    memory = '4G'
    time  = '4h'
  }
}
executor{
  name = 'slurm'
  submitRateLimit = '1 sec'
  queueSize = 500
}
*/

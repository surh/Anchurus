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

// Command line options
params.indir = ''
params.meta = ''
params.outdir = 'output/'
params.phylocpus = 4

// Process params
indir = file(params.indir)

snps_file = "snv_catalogue/MGYG-HGUT-00022_snvs.tsv"

SNPS = Channel.fromPath("${params.indir}/*.tsv")
  .map{ snpfile -> tuple(snpfile.name.replaceAll('_snvs\\.tsv$', ""),
    file(snpfile)) }

process snps2corealns{
  label 'py3'
  tag "$spec"
  publishDir "${param.outdir}/corelans", mode: 'rellink'

  input:
  tuple spec, file(snpfile) from SNPS
  file meta from params.meta

  output:
  tuple spec, file("${spec}.core.aln.fasta") into ALNS

  """
  ${workflow.projectDir}/uhgg_snps2corealns.py \
    --snps $snpfile \
    --meta $meta \
    --spec $spec \
    --output ${spec}.core.aln.fasta \
    --presence_thres 0.8 \
    --min_genomes 5
  """
}

process fasttree{
  label 'fasttree'
  tag "$spec"
  publishDir "${param.outdir}/phylo", mode: 'rellink'
  cpus params.phylocpus

  input:
  tuple spec, file(aln) from ALNS
  val phylocpus from params.phylocpus

  output:
  file "${spec}.core.gtr.tre"

  """
  export OMP_NUM_THREADS=$phylocpus

  FastTree -nt -gtr < $aln > ${spec}.core.gtr.tre
  """
}

// Example nextflow.config
/*
process{
  maxForks = 100
  stageInMode = 'rellink'
  errorStrategy = 'finish'
  withLabel: 'py3'{
    module = 'anaconda'
    conda = '/opt/modules/pkgs/anaconda/4.8/envs/fraserconda'
  }
  withLabel: 'fasttree'{
    module = "FastTree"
    time = '10h'
  }
}
executor{
  name = 'slurm'
  submitRateLimit = '1 sec'
  queueSize = 500
}
*/

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

// Nextflow pipeline that starts with SNVs in MIDAS format and produces
// RERconverge tests for a binary phenotype.
// First step is just to get the alignments.


// parameters
params.midas_dir = ""
params.genomes_dir = ""
params.map_dir = ""
params.min_cov = 0.8
params.outdir = "output/"


map_dir = file(params.map_dir)
genomes_dir = file(params.genomes_dir)
INDIRS = Channel.fromPath("${params.midas_dir}/*", type: 'dir')
  .map{spec -> tuple(spec.fileName,
    file(spec),
    file("${map_dir}/${spec.fileName}.map.txt"))}

process alns_from_metagenomes{
  label 'r'
  publishDir "${params.outdir}/gene_alns/",
    pattern: "$spec",
    mode: 'rellink'

  input:
  set spec, file(spec_dir), file(map_file) from INDIRS.filter{it[1].exists()}
  file genomes_dir

  output:
  file "$spec" into ALNDIR

  """
  ${workflow.projectDir}/all_alns_from_metagenomes.r \
    $spec_dir \
    $genomes_dir \
    --min_cov ${params.min_cov} \
    --map_file $map_file \
    --outdir $spec \
    --type single
  """
}


// Example nextflow config
/*
process{
  maxForks = 100
  stageInMode = 'rellink'
  withLabel: 'r'{
    module = 'R/3.6.1'
    memory = '5G'
    time  = '12h'
  }
  withLabel: 'py3'{
    module = 'anaconda'
    conda = '/opt/modules/pkgs/anaconda/3.6/envs/fraserconda'
  }
  withLabel: 'fasttree'{
    module = 'FastTree/2.1.10'
    time = '12h'
  }
}
executor{
  name = 'slurm'
  submitRateLimit = '1 sec'
  queueSize = 500
}
*/

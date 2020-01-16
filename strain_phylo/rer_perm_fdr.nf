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

// Calculate FDR or RERtest from permutations

// params
params.rertest = 'rertest/'
params.rerperms = 'rerperms/'
params.outdir = 'output/'

// Match rertest and rerperms
rertest_dirs = file(params.rertest)
rerperms_dirs = file(params.rerperms)
RERTEST = Channel.fromPath("$rertest_dirs/*", type:'dir', maxDepth:1)
  .map{dir -> spec = dir.name;
    tuple(spec, file("$rertest_dirs/$spec/${spec}.cors.txt"))}
RERPERMS = Channel.fromPath("$rerperms_dirs/*", type:'dir', maxDepth:1)
  .map{dir -> spec = dir.name;
    tuple(spec, file(dir))}
IN = RERTEST.join(RERPERMS)

process rer_fdr{
  label 'r'
  tag "$spec"
  publishDir params.outdir

  input:
  tuple val(spec), file(rertest), file(rerperms) from IN

  output:
  tuple spec, file("${spec}.rer.fdr.txt")

  """
  stitch_file,r ${workflow.projectDir}/rer_perm_fdr.r \
    $rerperms \
    $rertest \
    ${spec}.rer.fdr.txt
  """
}

// Example nextflow.config
/*
process{
  maxForks = 100
  stageInMode = 'rellink'
  errorStrategy = 'finish'
  withLabel: 'r'{
    module = 'R/3.6.1'
    memory = '4G'
    time  = '2h'
  }
}
executor{
  name = 'slurm'
  submitRateLimit = '1 sec'
  queueSize = 500
}
*/

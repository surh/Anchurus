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

// Nextflow pipeline that runs benchmark_imputation.r script from HMVAR
// on a list of directories.

// Parameters
params.dirs = ''
params.bin = ''
params.map = 'map.txt'
params.seed = 76543
params.njobs = 10
param.outdir = "output/"
params.memory = "20G"
params.time = "10:00:00"
params.queue = 'hbfraser,hns'


// Process params
if(params.bin == ''){
  error "You must provide a path to the benchmark_imputation.r executable."
}
dirs = file(dirs)
map = file(params.map)

// Read list of dirs
DIRS = Channel.fromPath(dirs).
  splitCsv(sep: "\t").
  map{row -> tuple(file(row[0]))}

process benchmark_imputation{
  label 'r'
  memory params.memory
  time params.time
  makForks params.njobs
  queue params.queue
  cpus 1

  input:
  set file(specdir) from DIRS
  file map

  output:
  file 'benchmark_imputation/imputation_results.txt' into IMPRES
  file 'benchmark_imputation/summary_stats.txt' into IMPSUM

  """
  ${params.bin} \
    $specdir \
    --map_file $map \
    --outdir benchmark_imputation \
    --m 5 \
    --seed ${params.seed}
  """
}

process collect_results{
  label 'py3'
  cpus 1
  memory '5G'
  time '1:00:00'
  queue params.queue
  publishDir params.outdir

  input:
  file "*.txt" from IMPRES.collect()

  output:
  file "imputation_results.txt"

  """
  ${workflow.projectDir}/../utils/cat_tables.py \
    *.txt \
    --outfile imputation_results.txt
  """
}

process collect_summaries{
  label 'py3'
  cpus 1
  memory '1G'
  time '00:30:00'
  queue params.queue
  publishDir params.outdir

  input:
  file "*.txt" from IMPSUM.collect()

  output:
  file "imputation_summary.txt"

  """
  ${workflow.projectDir}/../utils/cat_tables.py \
    *.txt \
    --outfile imputation_summary.txt
  """
}

// Example nextflow.config
/*
process {
  executor = 'slurm'
  withLabel: py3 {
    module = 'fraserconda'
  }
  withLabel: r {
    module = 'R/3.5.1server'
  }
}
*/

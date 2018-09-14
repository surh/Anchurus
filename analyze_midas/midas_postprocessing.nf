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

// Nextflow pipeline that processes the different genomes after running
// 'midas_merge.py snps'. It selects samples of interest, and renames
// SNPs so they have unique names across bacterial genomes.

// Params
params.indir = 'midas.merge/'     // 'midas_merge.py snps' output
params.samples = 'samples.txt'    // Must NOT have header row
params.genomes = 'genome_ids.txt' // Must have header row
params.outdir = './'
params.bindir = '/home/sur/micropopgen/src/Anchurus/analyze_midas/'
params.queue = 'hbfraser,hns,owners'
params.memory = '2.GB'
params.time = '20.m'
params.njobs = 200

// Process params
samples = file(params.samples)
genomes = file(params.genomes)

// Read genomes file
reader = genomes.newReader()
GENOMES = []
i = 0
while(str = reader.readLine()){

  row = sample = str.split("\t")

  // Check and skip header
  if(i == 0){
    if(row[0] != 'species' || row[1] != 'id'){
      throw new Exception("Wrong header")
    }

    i = i + 1
    continue
  }

  // input specdir
  dir = file("${params.indir}/${row[0]}/")
  id = row[1]
  freq_file = file("${dir}/snps_freq.txt")
  depth_file = file("${dir}/snps_depth.txt")
  info_file = file("${dir}/snps_info.txt")

  GENOMES = GENOMES + [tuple(id, dir, freq_file, depth_file, info_file)]
}


process post_process_midas_snps{
  module 'R'
  cpus 1
  // time {params.time + ((task.attempt - 1 ) * 10.m)}
  // memory { params.memory + ((task.attempt - 1) * 2.GB) }
  time {20.m + ((task.attempt - 1 ) * 10.m)}
  memory { 2.GB + ((task.attempt - 1) * 2.GB) }
  maxForks params.njobs
  queue params.queue
  publishDir params.outdir, mode: 'copy'
  errorStrategy 'retry'
  maxRetries 2

  input:
  set id, dir, freq, depth, info from GENOMES
  file samples
  // file genomes

  output:
  file 'snps_freq.txt' optional true into FREQS
  file 'snps_depth.txt' optional true into DEPTHS
  file 'snps_info.txt' optional true into INFOS

  """
  ${params.bindir}/midas_postprocessing.r \
    --outdir ./ \
    --genome_ids ${id} \
    --overwrite \
    specdir \
    ${dir} \
    ${samples}
  """
}

process concatenate_freqs{
  module 'fraserconda'
  publishDir params.outdir, mode: 'copy'
  cpus 1
  time '2:00:00'
  memory params.memory
  maxForks params.njobs
  queue params.queue
  publishDir params.outdir, mode: 'copy'
  errorStrategy 'retry'
  maxRetries 2


  input:
  file '*.txt' from FREQS.collect()

  output:
  file 'freqs.txt'

  """
  ${params.bindir}/concatenate_check_header.py \
    --infiles *.txt \
    --outfile freqs.txt
  """
}

process concatenate_depth{
  module 'fraserconda'
  publishDir params.outdir, mode: 'copy'
  cpus 1
  time '2:00:00'
  memory params.memory
  maxForks params.njobs
  queue params.queue
  publishDir params.outdir, mode: 'copy'
  errorStrategy 'retry'
  maxRetries 2

  input:
  file '*.txt' from DEPTHS.collect()

  output:
  file 'depth.txt'

  """
  ${params.bindir}/concatenate_check_header.py \
    --infiles *.txt \
    --outfile depth.txt
  """
}

process concatenate_info{
  module 'fraserconda'
  publishDir params.outdir, mode: 'copy'
  cpus 1
  time '2:00:00'
  memory params.memory
  maxForks params.njobs
  queue params.queue
  publishDir params.outdir, mode: 'copy'
  errorStrategy 'retry'
  maxRetries 2

  input:
  file '*.txt' from INFOS.collect()

  output:
  file 'info.txt'

  """
  ${params.bindir}/concatenate_check_header.py \
    --infiles *.txt \
    --outfile info.txt
  """
}

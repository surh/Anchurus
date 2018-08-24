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

// Nextflow pipeline that submits sample fastq files to midas to obtain
// snp profiles

// Main parameters
params.samples = 'samples.txt'
params.indir = 'samples/'
params.outdir = 'midas/'
params.sample_col = 1
params.queue = 'hbfraser,bigmem,hns,owners'
params.memory = '10G'
params.time = '4:00:00'
params.cpus = 8
params.njobs = 200
params.species_cov = 3.0
params.mapid = 94.0
params.mapq = 20
params.baseq = 30
params.readq = 30
params.trim = 0
params.discard = false
params.baq = false
params.adjust_mq = false
params.steps = 'all'

// Process params
samples = file(params.samples)
sample_col = params.sample_col - 1

// Read samples file
reader = samples.newReader()
SAMPLES = []
while(str = reader.readLine()){
  // Extract sample and run IDs
  sample = str.split("\t")[sample_col]
  SAMPLES = SAMPLES + [tuple(sample,
    file("${params.indir}/${sample}_read1.fastq.bz2"),
    file("${params.indir}/${sample}_read2.fastq.bz2"))]
}

#!/usr/bin/env nextflow
// Copyright (C) 2017 Sur Herrera Paredes

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
// species profiles


// Main parameters
params.samples = 'samples.tst'
params.indir = 'samples/'
params.outdir = 'midas/'
params.sample_col = 1
params.logdir = 'logs'
params.queue = 'hbfraser,owners,bigmem,hns,normal'
params.memory = '10G'
params.time = '4:00:00'

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

  GENOMES = GENOMES + [tuple(id, dir)]
}

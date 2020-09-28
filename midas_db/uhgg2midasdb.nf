#!/usr/bin/env nextflow
// Copyright (C) 2020 Sur Herrera Paredes

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

// Command line parameters
params.indir = "input/"
params.outdir = "output/"

// Process parameters
indir = file(params.indir)

GENOMEDIRS = Channel.fromPath("$indir/*/*", type: 'dir')
  .map{specdir -> tuple(specdir.name, file("${specdir}/genome/"))}

process prepare{
  tag "$spec"
  label "bedtools"
  publishDir params.outdir, mode: 'copy'

  input:
  tuple spec, file(genomedir) from GENOMEDIRS

  output:
  tuple spec, file("$spec/") into SPECDIRS

  """
  # Create outdir
  mkdir $spec

  # Copy assembly fna
  cp $genomedir/${spec}.fna $spec/

  # Copy CDS faa
  cp $genomedir/${spec}.faa $spec/

  # GFF to BED
  cat $genomedir/${spec}.gff | \
    cut -f 1,4,5,7,9 | \
    sed 's/;/\t/' | \
    cut -f 1,2,3,4,5 | \
    sed 's/ID=//' | \
    awk '{print \$1 "\t" \$2 - 1  "\t" \$3 "\t" \$5 "\t.\t" \$4}' > ${spec}.bed

  # Extract gene fasta (ffn)
  bedtools getfasta \
    -fi $genomedir/${spec}.fna \
    -bed ${spec}.bed \
    -s \
    -name | \
    sed 's/([+-])\$//' > $spec/${spec}.ffn

  # Genes table
  echo -e "gene_id\tscaffold_id\tstart\tend\tstrand\tgene_type" > \
    $spec/${spec}.genes
  cat $genomedir/${spec}.gff | \
    cut -f 1,3,4,5,7,9 | \
    sed 's/;/\t/' | \
    cut -f 1,2,3,4,5,6 | \
    sed 's/ID=//' | \
    awk '{print \$6 "\t" \$1 "\t" \$3 "\t" \$4 "\t" \$5 "\t" \$2}' >> $spec/${spec}.genes
  """
}

// Example nextflow.config
/*
process{
  queue = 'hbfraser,hns, owners'
  maxForks = 500
  errorStrategy = 'finish'
  stageInMode = 'rellink'
  time = '1h'
  memory = '1G'
  withLabel: 'bedtools'{
    module = 'bedtools/2.28.0'
  }
}

executor{
  name = 'slurm'
  queueSize = 500
  submitRateLitmit = '1 sec'
}
*/

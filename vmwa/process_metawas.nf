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

// Takes output of metawas.nf and performs the next steps. Requires HMVAR,
//  stitch_file.r in PATH, and the location of HMVAR scripts. And bedtools

// parameters
params.midas_db = ''
params.scripts = ''
params.metawas = 'lmmres/'
params.outdir = 'output/'
params.genomes_file = 'genomes.txt'
params.metawas_pval = 1e-6
params.metawas_suffix = '_lmm.results.txt'
params.filter = false

// process Params
metawas_dir = file(params.metawas)
genomes_file = file(params.genomes_file)

 // Read genomes file
reader = genomes_file.newReader()
GENOMES = []
while( line = reader.readLine() ) {
  GENOMES = GENOMES + [tuple(line,
    file("${metawas_dir}/${line}${params.metawas_suffix}"),
    file("${params.midas_db}/rep_genomes/${line}/genome.features.gz"))]
}


process manhattans{
  publishDir "${params.outdir}/manhattans/", mode: 'rellink'
  label 'r'

  input:
  file metawas_dir

  output:
  file "manhattans/*.manhattan.lmm.png"
  file "manhattans/*.manhattan.lmmpcs.png"

  """
  stitch_file.r ${params.scripts}/manhattan_type.r \
    $metawas_dir
    'manhattans/'
    $params.metawas_pval
    $params.metawas_suffix
  """
}

// Copied from micropopgen snps_to_genes.nf
process snps_to_genes{
  publishDir "params.outdir/closest/", mode: 'rellink', pattern: "*.closest"
  maxForks 20

  input:
  set genome, lmm_file, feat_file from GENOMES

  output:
  file 'snps.bed' optional true
  file 'genome.features.bed'
  file "${genome}.closest" optional true

  script:
  if(params.filter)
    """
    # Convert snps to BED
    #awk '(\$6 <= ${params.metawas_pval}){print \$1 "\\t" \$2 "\\t" \$2}' \
    #  ${lmm_file} | sort -k1,1 -k2,2n > snps.bed
    awk '(\$14 != "none" && \$3 != "ps"){print \$1 "\\t" \$3 "\\t" \$3}' \
    ${lmm_file} | sort -k1,1 -k2,2n > snps.bed

    # Convert features to BED
    zcat ${feat_file} | awk '{print \$2 "\t" \$3 "\t" \$4 "\t" \$1}' | \
      grep -v scaffold_id | sort -k1,1 -k2,2n > genome.features.bed

    # Find closest
    closestBed -D a -a snps.bed -b genome.features.bed > ${genome}.closest
    """
  else
    """
    # Convert snps to BED
    awk '!(\$1 == "chr" && \$2 == "rs" && \$3 == "ps"){print \$1 "\\t" \$3 "\\t" \$3}' \
    ${lmm_file} | sort -k1,1 -k2,2n > snps.bed

    # Convert features to BED
    zcat ${feat_file} | awk '{print \$2 "\t" \$3 "\t" \$4 "\t" \$1}' | \
      grep -v scaffold_id | sort -k1,1 -k2,2n > genome.features.bed

    # Find closest
    closestBed -D a -a snps.bed -b genome.features.bed > ${genome}.closest
    """
}

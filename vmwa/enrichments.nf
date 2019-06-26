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

params.files = ''
params.dist_thres = 500
params.count_thres = 3
params.outdir = 'enrichments/'
params.suffix = '_lmm.results.txt'
params.score_column = 'p_lrt.lmmpcs'

// Process
files = file(params.files)

// Read list of dirs
FILES = Channel.fromPath(files).
  splitCsv(sep: "\t").
  map{row -> tuple(row[0], file(row[1]), file(row[2]), file(row[3]))}.
  into{FILES_GO; FILES_KO; FILES_OG}


process go_enrichments{
  label 'r'
  publishDir "${params.outdir}/GO/", mode: 'rellink'

  input:
  set spec, file(lmm), file(closest), file(annots) from FILES_GO

  output:
  file "*.enrichments.txt"

  """
  Rscript ~/micropopgen/src/HMVAR/inst/bin/annotation_enrichments.r \
    $lmm \
    ./ \
    --suffix ${params.suffix} \
    --closest $closest \
    --annotations $annots \
    --dist_thres ${params.dist_thres} \
    --min_size ${params.count_thres} \
    --annot_column GO_terms \
    --score_column ${params.score_column} \
    --gene_score min \
    --alternative less \
    --method gsea
  """
}

process ko_enrichments{
  label 'r'
  publishDir "${params.outdir}/KO/", mode: 'rellink'

  input:
  set spec, file(lmm), file(closest), file(annots) from FILES_KO

  output:
  file "*.enrichments.txt"

  """
  Rscript ~/micropopgen/src/HMVAR/inst/bin/annotation_enrichments.r \
    $lmm \
    ./ \
    --suffix ${params.suffix} \
    --closest $closest \
    --annotations $annots \
    --dist_thres ${params.dist_thres} \
    --min_size ${params.count_thres} \
    --annot_column KEGG_KOs \
    --score_column ${params.score_column} \
    --gene_score min \
    --alternative less \
    --method gsea
  """
}

process og_enrichments{
  label 'r'
  publishDir "${params.outdir}/OG/", mode: 'rellink'

  input:
  set spec, file(lmm), file(closest), file(annots) from FILES_OG

  output:
  file "*.enrichments.txt"

  """
  Rscript ~/micropopgen/src/HMVAR/inst/bin/annotation_enrichments.r \
    $lmm \
    ./ \
    --suffix ${params.suffix} \
    --closest $closest \
    --annotations $annots \
    --dist_thres ${params.dist_thres} \
    --min_size ${params.count_thres} \
    --annot_column OGs \
    --score_column ${params.score_column} \
    --gene_score min \
    --alternative less \
    --method gsea
  """
}

// example nextflow.config
/*
process {
  executor = slurm
  maxForks = 40
  label: r{
    module = 'R/3.5.1server'
  }
}
*/

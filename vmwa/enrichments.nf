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

// Process
files = file(params.files)

// Read list of dirs
FILES = Channel.fromPath(files).
  splitCsv(sep: "\t").
  map{row -> tuple(row[0], file(row[1]), file(row[2]), file(row[3]))}


process genome_enrichments{
  publishDir params.outdir, mode: 'copy'
  makForks = 4

  input:
  set spec, file(lmm), file(closest), file(annots) from FILES

  output:
  file "enrichments/*"

  """
  Rscript ~/micropopgen/src/HMVAR/inst/scripts/metawas_enrichments.r \
    $lmm \
    $closest \
    $annots \
    --dist_thres ${params.dist_thres} \
    --count_thres ${params.count_thres} \
    --outdir enrichments/ \
    --prefix $spec
  """
}

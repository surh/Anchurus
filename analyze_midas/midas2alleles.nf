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

// Takes output from midas_merge.py snps and assigns genotypes, and filter
// genotypes for downstream analysis

// Parameters
params.indir = ""
params.outdir = "output/"

// QP Parameters
params.min_depth = 5
params.min_snv_prop = 0.8
params.maf_thres = 0.2

// Filter params
params.map_file = "map.txt"
params.min_sample_per_group = 5
params.min_snv_prop_per_sample = 0.5
params.min_core_gene_prev = 0.8
params.min_core_gene_cov = 0.8
params.min_core_genes = 0.8

// Proess params
map_file = file(params.map_file)
indir = file(params.indir)
outdir = file(params.outdir)

println(indir)

SNVDIRS = Channel.fromPath("$indir/*", type: 'dir')
  .map{snvdir -> tuple(snvdir.name,
    file(snvdir))}

process qp_genotypes{
  label 'r'
  tag "$spec"
  publishDir "$outdir/qp_alleles", mode: 'rellink',
    pattern: "snps_alleles.txt",
    saveAs: {"$spec/snps_alleles.txt"}

  input:
  tuple spec, file(snvdir) from SNVDIRS
  val min_depth from params.min_depth
  val min_snv_prop from params.min_snv_prop
  val maf_thres from params.maf_thres

  output:
  tuple spec, file("snps_alleles.txt") into QPALLELES
  tuple spec, file("$indir/snps_info.txt") into SNVINFOS

  """
  Rscript ${workflow.projectDir}/qp_snvs.r \
    $indir \
    snps_alleles.txt \
    --min_depth $min_depth \
    --min_svn_prop $min_snv_prop \
    --maf_thres $maf_thres
  """
}

process filter_genotypes{
  label 'r'
  tag "$spec"
  publishDir "$outdir/fltered_alleles", mode: 'rellink',
    pattern: "output",
    saveAs: {"$spec"}

  input:
  tuple spec, file("snps_alleles.txt"),
    file("snps_info.txt") from QPALLELES.join(SNVINFOS)
  val map_file from params.map_file
  val min_sample_per_group from params.min_sample_per_group
  val min_snv_prop_per_sample from params.min_snv_prop_per_sample
  val min_core_gene_prev from params.min_core_gene_prev
  val min_core_gene_cov from params.min_core_gene_cov
  val min_core_genes from params.min_core_genes

  output:
  tuple spec, file("output") optional true into FILTERED

  """
  Rscript ${workflow.projectDir}/filter_genotypes.r \
    snps_info.txt \
    $map_file \
    output \
    --min_sample_per_group $min_sample_per_group \
    --min_snv_prop_per_sample $min_snv_prop_per_sample \
    --min_core_gene_prev $min_core_gene_prev \
    --min_core_gene_cov $min_core_gene_cov \
    --min_core_genes $min_core_genes \
    --process_info yes
  """
}

// Example nextflow.config
/*
process{
  partition = 'hbfraser,hns,owners'
  maxForks = 200
  stageInMode = 'rellink'
  errorStrategy = 'finish'
  withLabel: 'r'{
    module = 'R/3.6.1'
    memory = '4G'
    time  = '4h'
  }
}
executor{
  name = 'slurm'
  submitRateLimit = '1 sec'
  queueSize = 500
}
*/

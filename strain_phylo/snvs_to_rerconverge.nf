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

// Nextflow pipeline that starts with SNVs in MIDAS format and produces
// RERconverge tests for a binary phenotype.
// First step is just to get the alignments.

/*
--midas_dir (optional)
MIDAS merge output. Must match species file. Starts pipeline from the
beginning.
--genomes_dir
Representative genomes dir. One genome per species with name matching species
fileName
--map_dir
Directory with one map per species. Should include only samples to keep.
Filename mut match <species name>.map.txt
--alns_dir (optional)
Directory one directory per species containing all alignments for that
species. Directories must match species name. Skips alns_from_metagenomes.
--gene_trees_dir (optional)
Directory with one directory per species containing all trees for that
species. Directories must match species names.
--master_trees_dir
Directory with core phylogeny of each species. Files must be named
<species name>.tre
--cov_dir
Directory with gene coverage matrices per species. Files must be named
<species name>.gene_coverage.txt
--focal_phenotype
Value of phenotype to compare against all other samples.
--snvs
Indicates which SNVs to keep. Either all ('all'), non-synonymous ('ns')
or synonymous ('s')
--min_cov
Minimum coverage for genes in a given strain to be analyzed.
--outdir
Directory where to place output
--baseml_threads
Number of threads to split all baseml jobs. Each job uses just one thread.
Default: 1.
*/

// parameters
params.midas_dir = ""
params.genomes_dir = ""
params.map_dir = ""
params.master_trees_dir = ""
params.cov_dir = ""
params.focal_phenotype = "USA"
params.snvs = 'all'
params.min_cov = 0.8
params.outdir = "output/"
params.baseml_threads = 1

// Optional parameters
params.alns_dir = ""
params.gene_trees_dir = ""


map_dir = file(params.map_dir)
genomes_dir = file(params.genomes_dir)
INDIRS = (params.midas_dir == ""
  ? Channel.empty()
  : Channel.fromPath("${params.midas_dir}/*", type: 'dir')
    .map{spec -> tuple(spec.fileName,
      file(spec),
      file("${map_dir}/${spec.fileName}.map.txt"))})

// Create channel with gene level alignments
ALNDIR = (params.alns_dir == ""
  ? Channel.empty()
  : Channel.fromPath("${params.alns_dir}/*", type: 'dir')
      .map{spec -> tuple(spec.name, file(spec))})

// Channel with master trees
MASTERTREE = Channel.fromPath("${params.master_trees_dir}/*.tre")
  .map{filename -> tuple(filename.name.replace('.tre', ''), file(filename))}
  .into{MT_BASEML; MT_RER}
// MT_BASEML.subscribe{println it}

// Channel with gene coverages
COV = Channel.fromPath("${params.cov_dir}/*.gene_coverage.txt")
  .map{filename -> tuple(filename.name.replace('.gene_coverage.txt', ''), file(filename))}

// Channel with gene level trees
// Create channel with gene level alignments
GENETREESDIR = (params.gene_trees_dir == ""
  ? Channel.empty()
  : Channel.fromPath("${params.gene_trees_dir}/*", type: 'dir')
      .map{spec -> tuple(spec.name, file(spec))})

SPECMAPS = Channel.fromPath("${map_dir}/*")
  .map{filename -> tuple(filename.name.replace('.map.txt', ''), file(filename))}

process alns_from_metagenomes{
  label 'r'
  tag "$spec"
  publishDir "${params.outdir}/gene_alns/",
    pattern: "output",
    saveAs: {"${spec}/"},
    mode: 'rellink'

  input:
  tuple val(spec), file(midas_dir), file(map_file) from INDIRS
  file genomes_dir
  val snvs from params.snvs

  output:
  tuple val(spec), file("output") optional true into MIDAS2ALNS

  when:
  map_file.exists()

  """
  ${workflow.projectDir}/all_alns_from_metagenomes.r \
    $midas_dir \
    $genomes_dir \
    --min_cov ${params.min_cov} \
    --map_file $map_file \
    --outdir output/ \
    --type single \
    --snvs $snvs
  """
}
println "============="
// ALNDIR.mix(MIDAS2ALNS).subscribe{println it}
// ALNDIR.mix(MIDAS2ALNS).join(MT_BASEML).subscribe{println it}
MIDAS2ALNS.mix(ALNDIR).join(MT_BASEML).subscribe{println it}
// MT_BASEML.subscribe{println it}
// join(MT_BASEML).join(COV)
println "#############"

// process baseml{
//   label 'baseml'
//   tag "$spec"
//   cpus params.baseml_threads
//   publishDir "${params.outdir}/gene_trees/",
//     pattern: "output/gene_trees",
//     saveAs: {"${spec}/"},
//     mode: 'rellink'
//
//   input:
//   tuple spec, file("alns_dir"), file(master_tree), file(cov) from ALNDIR.mix(MIDAS2ALNS).join(MT_BASEML).join(COV)
//
//   output:
//   tuple val(spec), file("output/gene_trees/") into ALNS2BASEML
//
//   """
//   ${workflow.projectDir}/baseml_all_genes.py \
//     --aln_dir alns_dir/ \
//     --cov_file $cov \
//     --master_tree $master_tree \
//     --outdir output/ \
//     --min_cov ${params.min_cov} \
//     --baseml baseml \
//     --cpus ${params.baseml_threads} \
//     --resume
//   """
//
// }
//
// // println "============="
// // GENETREESDIR.mix(ALNS2BASEML).subscribe{println it}
// // GENETREESDIR.subscribe{println it}
// // TEST = GENETREESDIR.mix(ALNS2BASEML)
//
// process trees2tab{
//   tag "$spec"
//   publishDir "${params.outdir}/tree_tabs",
//     pattern: 'trees_tab.txt',
//     saveAs: {"${spec}.trees.txt"},
//     mode: 'rellink'
//
//   input:
//   tuple val(spec), file("trees") from GENETREESDIR.mix(ALNS2BASEML)
//   // tuple val(spec), file("trees") from TEST
//
//   output:
//   tuple val(spec), file("trees_tab.txt") into TREETABS
//
//   """
//   for f in trees/*.tre; \
//     do echo -e `basename \$f`"\\t"`cat \$f`; \
//     done | sed 's/\\.baseml\\.tre//' > trees_tab.txt
//   """
// }
//
// process rertest{
//   tag "$spec"
//   label 'r'
//   publishDir "${params.outdir}/rertest/",
//     pattern: "output",
//     saveAs: {"${spec}/"},
//     mode: 'rellink'
//
//   input:
//   tuple val(spec),
//     file("trees_tab.txt"),
//     file("master_tree.tre"),
//     file("map.txt") from TREETABS.join(MT_RER).join(SPECMAPS)
//   val pheno from params.focal_phenotype
//
//   output:
//   path "output"
//   tuple val(spec), file("output/${spec}.cors.txt") into RERCORS
//   tuple val(spec), file("output/${spec}.rerw.dat") into RERWS
//   tuple val(spec), file("output/${spec}.Trees.dat") into RERTREES
//
//   """
//   ${workflow.projectDir}/rertest.r \
//     trees_tab.txt \
//     master_tree.tre \
//     --map_file map.txt \
//     --outdir output \
//     --focal_phenotype $pheno \
//     --spec $spec
//   """
// }

// Example nextflow.config
/*
process{
  maxForks = 100
  stageInMode = 'rellink'
  errorStrategy = 'finish'
  withLabel: 'r'{
    module = 'R/3.6.1'
    memory = '5G'
    time  = '12h'
  }
  withLabel: 'py3'{
    module = 'anaconda'
    conda = '/opt/modules/pkgs/anaconda/3.6/envs/fraserconda'
  }
  withLabel: 'fasttree'{
    module = 'FastTree/2.1.10'
    time = '12h'
  }
  withLabel: 'baseml'{
    module = "anaconda:paml/4.9i"
    conda = '/opt/modules/pkgs/anaconda/3.6/envs/fraserconda'
    time = '200h'
  }
}
executor{
  name = 'slurm'
  submitRateLimit = '1 sec'
  queueSize = 500
}
*/

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

// Nextflow pipeline that runs RERConverge testst on a permuted set
// of mapping files

/* Options:
--maps_dirs
A directory with permuted map directories. There must be one directory per
species named after the species, and each species directory must contain
all permuted maps from tha species. Permuted maps must be named map_<##>.txt.
--tree_tabs_dir
A directory with gene trees in tab format. Files must be named
<species>.trees.txt
--master_trees_dir
A directory with master species trees. Files must be named <species>.tre
--outdir
Directory for output
--focal_phenotype
Phenotype group that will be converted to 1 for RER test, everything else
will be 0.
*/

// parameters
params.maps_dirs = ''
params.tree_tabs_dir = ''
params.master_trees_dir = ''
params.rertest_dirs = ''
params.outdir = 'output/'
params.focal_phenotype = "USA"

// Process inputs
maps_dir = file(params.maps_dirs)
tree_tabs_dir = file(params.tree_tabs_dir)
master_trees_dir = file(params.master_trees_dir)
rertest_dirs = file(params.rertest_dirs)

INPUTS = Channel.fromPath("$maps_dir/**", type: 'file', maxDepth: 2)
  .map{map_file -> tuple(map_file.getParent().name,
    map_file.name.replaceAll("^map_","").replaceAll('\\.txt$', ""),
    file(map_file))}
  .map{spec, perm, map_file -> tuple(spec, perm, map_file,
    file("$tree_tabs_dir/${spec}.trees.txt"),
    file("$master_trees_dir/${spec}.tre"))}
// INPUTS.subscribe{println it}

// TREETABS = Channel.fromPath("$tree_tabs_dir/*", type: 'file')
//   .map{trees_file -> tuple(trees_file.name.replaceAll('\\.trees\\.txt$',''),
//     file(trees_file))}
// TREETABS.subscribe{println it}
// There is no left join!!!
// MAPS.join(TREETABS, remainder: true).subscribe{println it}
// MAPS.join(TREETABS, remainder: true).filter{items -> items[1] != null}.subscribe{println it}

// Results from RERConverge on real labels
RERTEST = Channel.fromPath("$rertest_dirs/*", type:'dir', maxDepth:1)
  .map{dir -> spec = dir.name;
    tuple(spec, file("$rertest_dirs/$spec/${spec}.cors.txt"))}

process rertest{
  tag "${spec}.perm_$nperm"
  label 'r'
  publishDir "${params.outdir}/rertest/${spec}/",
    pattern: "output/perm_${nperm}.cors.txt",
    saveAs: {"perm_${nperm}.cors.txt"},
    mode: 'rellink'

  input:
  tuple val(spec),
    val(nperm),
    file("map.txt"),
    file("trees_tab.txt"),
    file("master_tree.tre") from INPUTS
  val pheno from params.focal_phenotype

  output:
  path "output"
  tuple val(spec), file("output/perm_${nperm}.cors.txt") into RERCORS
  tuple val(spec), file("output/perm_${nperm}.rerw.dat") into RERWS
  tuple val(spec), file("output/perm_${nperm}.Trees.dat") into RERTREES

  """
  Rscript ${workflow.projectDir}/rertest.r \
    trees_tab.txt \
    master_tree.tre \
    --map_file map.txt \
    --outdir output \
    --focal_phenotype $pheno \
    --spec perm_$nperm
  """
}



// RERPERMS = Channel.fromPath("$rerperms_dirs/*", type:'dir', maxDepth:1)
//   .map{dir -> spec = dir.name;
//     tuple(spec, file(dir))}

process rer_fdr{
  label 'r'
  tag "$spec"
  publishDir params.outdir

  input:
  tuple spec, file(rertest),
    file("rerperms/") from RERTEST.join(RERCORS.groupTuple())

  output:
  tuple spec, file("${spec}.rer.fdr.txt")

  """
  Rscript ${workflow.projectDir}/rer_perm_fdr.r \
    "rerperms/" \
    $rertest \
    ${spec}.rer.fdr.txt
  """
}



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
}
executor{
  name = 'slurm'
  submitRateLimit = '1 sec'
  queueSize = 500
}
*/

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


// parameters
params.midas_dir = ""
params.genomes_dir = ""
params.map_dir = ""
params.alns_dir = ""
params.master_trees_dir = ""
params.cov_dir = ""
params.min_cov = 0.8
params.outdir = "output/"


map_dir = file(params.map_dir)
genomes_dir = file(params.genomes_dir)
INDIRS = (params.midas_dir == ""
  ? Channel.empty()
  : Channel.fromPath("${params.midas_dir}/*", type: 'dir')
    .map{spec -> tuple(spec.fileName,
      file(spec),
      file("${map_dir}/${spec.fileName}.map.txt"))})
// INDIRS = Channel.fromPath("${params.midas_dir}/*", type: 'dir')
//   .map{spec -> tuple(spec.fileName,
//     file(spec),
//     file("${map_dir}/${spec.fileName}.map.txt"))}

// Create channel with gene level alignments
ALNDIR = (params.alns_dir == ""
  ? Channel.empty()
  : Channel.fromPath(params.alns_dir)
      .map{spec -> tuple(spec.fileName, file(spec))})

process alns_from_metagenomes{
  label 'r'
  tag "$spec"
  publishDir "${params.outdir}/gene_alns/",
    pattern: "output",
    saveAs: {"${spec}/"},
    mode: 'rellink'

  input:
  set spec, file(midas_dir), file(map_file) from INDIRS
  file genomes_dir

  output:
  set spec, file("output") into MIDAS2ALNS

  when:
  map_file.exists()

  """
  ${workflow.projectDir}/all_alns_from_metagenomes.r \
    $midas_dir \
    $genomes_dir \
    --min_cov ${params.min_cov} \
    --map_file $map_file \
    --outdir output/ \
    --type single
  """
}

process baseml{
  label 'baseml'
  tag "$spec"
  publishDir "${params.outdir}/gene_trees/",
    pattern: "output",
    saveAs: {"${spec}/"},
    mode: 'rellink'

  input:
  set spec, file("alns_dir") from ALNDIR.mix(MIDAS2ALNS)
  path master_tree from "${workflow.launchDir}/${params.master_trees_dir}/${spec}.tre"
  path cov from "${workflow.launchDir}/${params.cov_dir}/${spec}.gene_coverage.txt"

  output:
  set spec, file("output") into ALNS2BASEML

  """
  ${workflow.projectDir}/baseml_all_genes.py
    --aln_dir alns_dir/ \
    --cov_file $cov \
    --master_tree $master_tree \
    --outdir output/ \
    --min_cov ${params.min_cov} \
    --baseml baseml
  """

}


process trees2tab{
  tag "$spec"
  publishDir "${params.outdir}/tree_tabs",
    pattern: 'trees_tab.txt',
    saveAs: {"${spec}.trees.txt"},
    mode: 'rellink'

  input:
  set spec, file("trees") from ALNS2BASEML

  output:
  set spec, file("trees_tab.txt")

  """
  for f in trees/*.tre; \
    do echo "\$f\\t"`cat \$f`; \
    done | sed 's/\\.baseml\\.tre//' > trees_tab.txt
  """
}

// Example nextflow config
/*
process{
  maxForks = 100
  stageInMode = 'rellink'
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
    module = "anaconda:paml/4.91"
    conda = '/opt/modules/pkgs/anaconda/3.6/envs/fraserconda'
    time = '48h'
  }
}
executor{
  name = 'slurm'
  submitRateLimit = '1 sec'
  queueSize = 500
}
*/

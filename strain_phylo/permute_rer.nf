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
--tree_tabs_dir
*/


// parameters
params.maps_dirs = ''
params.tree_tabs_dir = ''

// Process inputs
maps_dir = file(params.maps_dirs)
tree_tabs_dir = file(params.tree_tabs_dir)

MAPS = Channel.fromPath("$maps_dir/**", type: 'file', maxDepth: 2)
  .map{map_file -> tuple(map_file.getParent().name,
    map_file.name.replaceAll("^map_","").replaceAll('\\.txt$', ""),
    file(map_file))}
// MAPS.subscribe{println it}

TREETABS = Channel.fromPath("$tree_tabs_dir/*", type: 'file')
  .map{trees_file -> tuple(trees_file.name.replaceAll('\\.trees\\.txt$',''),
    file(trees_file))}
// TREETABS.subscribe{println it}

MAPS.cross(TREETABS).subscribe{println it}

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
--gene_trees_dir
*/


// parameters
params.maps_dirs = ''
params.gene_trees_dir = ''


// Process inputs
maps_dir = file(params.maps_dir)
gene_trees_dir = file(gene_trees_dir)

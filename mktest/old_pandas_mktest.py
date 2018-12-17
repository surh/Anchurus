#!/usr/bin/env python
# Copyright (C) 2017-2018 Sur Herrera Paredes

# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.

# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.

import pandas as pd
import numpy as np


def read_and_process_data(map_file, info_file, depth_file, freqs_file,
                          groups, cov_thres=1, nrows=float('inf')):
    """Reads MIDAS output files, selects gene sites and samples above threshold,
    determines mutation effect (s or n) and makes sure files are consistent
    with each other"""

    print("\tReading data...")
    # Read data
    info = pd.read_csv(info_file, sep="\t")
    depth = pd.read_csv(depth_file, sep="\t")
    freq = pd.read_csv(freqs_file, sep="\t")

    # Remove non gene sites
    ii = ~info.gene_id.isnull()
    info = info.loc[ii, :]
    depth = depth.loc[ii, :]
    freq = freq.loc[ii, :]

    # Remove site_id columns
    depth = depth.drop(axis=1, labels='site_id')
    freq = freq.drop(axis=1, labels='site_id')
    # print(depth.shape)
    # print(freq.shape)

    # subset for tests
    if nrows < float('inf'):
        nrows = int(nrows)
        info = info.head(nrows)
        depth = depth.head(nrows)
        freq = freq.head(nrows)

    print("\tClassifying sites into synonymous and non-synonymous")
    # Determine effect of sites (this is constant and indepentent of samples)
    info['Effect'] = info.apply(determine_mutation_effect, axis=1)

    # Check that sample names match between freq and depth
    if not all(freq.columns == depth.columns):
        raise ValueError("Columns don't match between freq and depth files")

    print("\tFiltering samples by group and coverage")
    # Read map file and select groups
    map = pd.read_csv(map_file, sep="\t")
    map.index = map.ID
    map = map.loc[map.Group.isin(groups), :].copy()

    # Remove samples from other groups
    ci = depth.columns.isin(map.ID)
    depth = depth.loc[:, ci]
    freq = freq.loc[:, ci]
    # print(depth.shape)
    # print(freq.shape)

    # Reorder map
    map = map.loc[depth.columns, :]

    # Calculate coverage in sites
    map['coverage'] = depth.mean(axis=0)

    # Remove samples below coverage
    ci = map.coverage >= cov_thres
    map = map.loc[ci, :]
    depth = depth.loc[:, map.index]
    freq = freq.loc[:, map.index]
    # print(depth.shape)
    # print(freq.shape)

    return map, freq, info, depth


def calculate_mk_oddsratio(map, info, depth, freq, depth_thres=1):
    """Determine if sites are fixed of polymorphic, calculate MK
    contingency table per gene, and calculate odds ratio"""

    if not all(map.index == depth.columns):
        raise ValueError("Samples in map and depth don't match")
    if not all(map.index == freq.columns):
        raise ValueError("Samples in map and freq don't match")
    if not all(freq.columns == depth.columns):
        raise ValueError("Samples in freq and depth don't match")

    print("\tDetermining if sites are fixed or polymorphic")
    # Determine type of mutation
    info['Type'] = determine_site_dist(map=map, depth=depth, freq=freq,
                                       info=info, depth_thres=depth_thres)
    # print(info.shape)
    # print(freq.shape)
    # print(depth.shape)
    # print(info.head(45))

    print("\tCalculate MK contingency table per gene")
    # Calculate MK contingency table per gene
    Genes = pd.DataFrame(columns=['Gene', 'Dn', 'Ds', 'Pn', 'Ps'])
    for g in info.gene_id.unique():
        dat = info.loc[info.gene_id == g, :].copy()
        # print(g)
        # print(dat.shape)
        tab = pd.crosstab(dat.Effect, dat.Type,
                          rownames=['Effect'],
                          colnames=['Type'])
        # print(tab)
        tab = tab.reindex(index=pd.Index(['n', 's']),
                          columns=pd.Index(['fixed', 'polymorphic']),
                          fill_value=0)
        s = pd.Series(g, index=['Gene']).append(
            tab.fixed).append(tab.polymorphic)
        Genes = Genes.append(pd.DataFrame([list(s)],
                                          columns=Genes.columns),
                             ignore_index=True)

    print("\tCalculate statistic")
    # Calculate ratio
    np.seterr(divide='ignore', invalid='ignore')
    Genes['ratio'] = pd.to_numeric(
        Genes.Dn * Genes.Ps) / pd.to_numeric(Genes.Ds * Genes.Pn)
    np.seterr(divide='raise', invalid='raise')
    Genes.replace(np.inf, np.nan, inplace=True)
    # Genes['hg.pval'] = Genes.apply(mktest_fisher_exact, axis=1)
    # Genes.head()

    return Genes, info


def determine_mutation_effect(r):
    """Mini function for apply, takes a series and checks if the mutation
    is synonymous (s) or non-synonymopus (n)"""

    # print(r)
    # ii = r.loc[['count_a', 'count_c', 'count_g', 'count_t']] > 0
    aa = pd.Series(r.amino_acids.split(sep=','),
                   index=['A', 'C', 'G', 'T'])
    # aa = np.array(r.amino_acids.split(sep=','))

    if aa[r.major_allele] == aa[r.minor_allele]:
        effect = 's'
    else:
        effect = 'n'

    return effect




def determine_site_dist(map, depth, freq, info, depth_thres=1):
    """For all sites, determine if they are fixed or polymorphic"""

    group1, group2 = map.Group.unique()
    Dist = []
    for i in range(info.shape[0]):
        # Add samples IDs as map header and match samples
        # in map with samples in depth

        # Create site data frame
        site = map.copy()
        site['depth'] = depth.loc[depth.index[i], map.index]
        site['freq'] = freq.loc[freq.index[i], map.index]

        # Remove samples without information for site
        site = site[site.depth >= depth_thres]
        # print(site)

        # Determine if it is polymorphic or fixed
        # site_crosstab = pd.crosstab(site.freq < 0.5, site.Group)
        site_crosstab = np.array([[sum((site.Group == group1) & (site.freq >= 0.5)),
                                   sum((site.Group == group2) & (site.freq >= 0.5))],
                                  [sum((site.Group == group1) & (site.freq < 0.5)),
                                   sum((site.Group == group2) & (site.freq < 0.5))]])
        # print(site_crosstab)
        if((site_crosstab.diagonal() == [0, 0]).all() or
           (np.fliplr(site_crosstab).diagonal() == [0, 0]).all()):
            mutation_type = 'fixed'
        else:
            mutation_type = 'polymorphic'

        Dist.append(mutation_type)
        # print(info.site_id.iloc[i], mutation_type)

    return(Dist)

if __name__ == "__main__":
    args = process_arguments()

    # Create file names
    info_file = ''.join([args.indir, '/snps_info.txt'])
    depth_file = ''.join([args.indir, '/snps_depth.txt'])
    freqs_file = ''.join([args.indir, '/snps_freq.txt'])

    # Prepare list of groups
    groups = [args.group1, args.group2]

    map, freq, info, depth = read_and_process_data(map_file=args.metadata_file,
                                                   info_file=info_file,
                                                   depth_file=depth_file,
                                                   freqs_file=freqs_file,
                                                   groups=groups,
                                                   cov_thres=args.min_cov,
                                                   nrows=args.nrows)
    Genes, info = calculate_mk_oddsratio(map=map, info=info,
                                         depth=depth, freq=freq,
                                         depth_thres=args.min_count)

    # Perform permutations if needed
    if args.permutations > 0:
        Perms = np.empty((Genes.shape[0], args.permutations + 1))
        Perms[:] = np.nan
        Perms[:, 0] = Genes.ratio
        np.random.seed(args.seed)
        print("Seed for permutations is {}".format(args.seed))
        for i in range(1, args.permutations + 1):
            # print(i)
            map_i = map.copy()
            map_i['Group'] = np.random.permutation(map_i.Group)
            Genes_i, info_i = calculate_mk_oddsratio(map=map_i,
                                                     info=info,
                                                     depth=depth,
                                                     freq=freq,
                                                     depth_thres=args.min_count)
            Perms[:, i] = Genes_i.ratio
        # Perms

        # Calculate permutation p-values
        Genes['nperm'] = args.permutations + \
            1 - pd.isnull(Perms).sum(axis=1)
        np.seterr(invalid='ignore', divide='ignore')
        Genes['P'] = np.greater_equal(Perms[:, np.repeat(
            0, args.permutations + 1)], Perms).sum(axis=1) / Genes['nperm']
        np.seterr(invalid='raise', divide='raise')

    Genes.to_csv(args.outfile, sep="\t", na_rep='NA', index=False)

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

# Imports
import os
import csv
import numpy as np
import pandas as pd
import scipy.stats as stats
import argparse


class GenomeSite:
    """A class for represintinc sites in genome that have potential SNPS"""

    def __init__(self, site_id, contig, position, ref_allele='',
                 major_allele='',
                 minor_allele='', locus_type='', gene_id='',
                 aminoacid_A='',
                 aminoacid_C='', aminoacid_G='', aminoacid_T=''):
        self.id = site_id
        self.contig = contig
        self.position = position
        self.ref_allele = ref_allele
        self.major_allele = major_allele
        self.minor_allele = minor_allele
        self.locus_type = locus_type
        self.gene_id = gene_id
        self.aminoA = aminoacid_A
        self.aminoC = aminoacid_C
        self.aminoG = aminoacid_G
        self.aminoT = aminoacid_T

    def codon_aminoacid(self, base):
        """This function returns the aminoacid that would be coded by
        the specified base in the site"""

        if base in ['A', 'a']:
            return(self.aminoA)
        elif base in ['C', 'c']:
            return(self.aminoC)
        elif base in ['G', 'g']:
            return(self.aminoG)
        elif base in ['T', 't']:
            return(self.aminoT)
        else:
            raise ValueError(("base must be one of the four canonical "
                              "nucleotides"))

    def substitution_type(self):
        """This function returns the type of subsitution encoded
        by the two alleles in the genomic site"""

        substitution_type = ''
        aa1 = self.codon_aminoacid(base=self.major_allele)
        aa2 = self.codon_aminoacid(base=self.minor_allele)
        if aa1 == aa2:
            substitution_type = 'synonymous'
        else:
            substitution_type = 'non-synonymous'

        return(substitution_type)


class Gene:
    """A class for representing a gene"""

    def __init__(self, gene_id, contig, start, end, strand=''):
        if(start > end):
            raise ValueError("Start cannot be greater than end")
        self.id = gene_id
        self.contig = contig
        self.start = int(start)
        self.end = int(end)
        self.strand = strand

    def extend(self, pos):
        pos = int(pos)
        if pos > self.end:
            self.end = pos
        elif pos < self.start:
            self.start = pos

    def info(self):
        print("===Gene===")
        print(">Gene id: {}".format(self.id))
        print(">Gene contig: {}".format(self.contig))
        print(">Gene start: {}".format(str(self.start)))
        print(">Gene end: {}".format(str(self.end)))


class MKtest:
    """A class for holding the McDonald-Kreitmant test"""

    def __init__(self, name, Ds=0, Dn=0, Ps=0, Pn=0):
        self.name = name
        self.Dn = Dn
        self.Ds = Ds
        self.Ps = Ps
        self.Pn = Pn

    def update(self, Ds=0, Dn=0, Ps=0, Pn=0):
        """Update the contigency matrix"""

        self.Dn += Dn
        self.Ds += Ds
        self.Ps += Ps
        self.Pn += Pn

    def mk_ratio(self, pseudocount=0):
        """Calculate the McDonald Kreitman ratio (Dn/Ds)/(Pn/Ps)"""

        num = (self.Dn + pseudocount) * (self.Ps + pseudocount)
        denom = (self.Ds + pseudocount) * (self.Pn + pseudocount)
        ratio = num / denom
        return ratio

    def alpha(self, pseudocount=0):
        """Calculate the Smith & Eyre-Walker alpha 1 - """
        ni = self.neutrality_index(pseudocount=pseudocount, log=False)
        alpha = 1 - ni
        return alpha

    def hg_test(self, pseudocount=0):
        """Hypergeometric (Fisher's exact) test"""

        res = stats.fisher_exact([[self.Ds + pseudocount,
                                   self.Ps + pseudocount],
                                  [self.Dn + pseudocount,
                                   self.Pn + pseudocount]])
        return res

    def g_test(self, correction, pseudocount=0):
        """G-test for independence. Original McDonald & Kreitman 1991
        suggestion"""

        # Create 2x2 contingency matrix
        mat = np.matrix([[self.Ds + pseudocount, self.Ps + pseudocount],
                         [self.Dn + pseudocount, self.Pn + pseudocount]])

        if correction == 'none':
            res = stats.chi2_contingency(observed=mat,
                                         lambda_="log-likelihood",
                                         correction=False)
        elif correction == 'yates':
            # apply yates correction, the default and only option
            # on scipy.stats
            res = stats.chi2_contingency(observed=mat,
                                         lambda_="log-likelihood",
                                         correction=True)
        elif correction == "williams":
            # Original correction used by McDonald & Kreitman (1991).
            # According to McDonald (same as above) biostat handbook,
            # it doesn't make much difference
            # (http://www.biostathandbook.com/small.html)
            g, p, df, e = stats.chi2_contingency(observed=mat,
                                                 lambda_="log-likelihood",
                                                 correction=False)

            # Calculate q correction. Only for 2 x 2 table
            n = mat.sum()
            q1 = (n * (1 / mat.sum(axis=1)).sum() - 1)
            q2 = (n * (1 / mat.sum(axis=0)).sum() - 1)
            q = 1 + q1 * q2 / (6 * n)

            # correct g and recalculate p-value
            g = g / q
            p = 1 - stats.chi2.cdf(g, df)

            # combine results
            res = [g, p, df, e]

        else:
            raise ValueError(("Correction must be one of 'none', 'yates' "
                              "or 'williams'"))

        return(res)

    def DoS(self, pseudocount):
        """Estimate Direction of Selection (DoS) from Stoletzki & Eyre-Walker
        2010"""
        first = self.Dn / (self.Dn + self.Ds)
        second = self.Pn / (self.Pn + self.Ps)
        DoS = first - second

        return DoS

    def neutrality_index(self, pseudocount=1, log=True):
        """Calculate neutrality index (Pn/Dn)/(Ps/Ds).
        Following Li et al. (2008), we add a psedocount and
        return the -log10(NI). Also referred to as haldane
        estimation."""

        num = (self.Pn + pseudocount) * (self.Ds + pseudocount)
        denom = (self.Ps + pseudocount) * (self.Dn + pseudocount)
        ni = num / denom

        if log:
            ni = -np.log10(ni)

        return(ni)


def calculate_contingency_tables(Map, args):
    """Take metadata dictionaries and location of MIDAS merge
    files, and calculate MK contingency tables. Ideally run
    after confirming existence of MIDAS files."""

    # Find total number of sites and genes present. Define an object
    # per site and gene.
    print("\tRead snps_info.txt")
    Genes, Sites = process_snp_info_file(args)
    # print("Number of sites: {}".format(str(len(Sites))))
    # print("Number of genes: {}".format(str(len(Genes))))

    # Remove sites that do not have at least one sample per group at
    # minimum depth. And find which samples are to be used per site (Counts).
    print("\tChose sites based on depth in groups to compare")
    Counts = process_snps_depth_file(args, Map, Sites)
    # print("Number of sites: {}".format(str(len(Sites))))
    # print("Number of genes: {}".format(str(len(Genes))))
    # print("Sites with counts: {}".format(str(len(Counts))))

    # Identify alleles and calculate contingency table
    print("\tRead frequencies and calculate")
    MK = process_snp_freq_file(args, Counts, Map, Sites)
    # print("Number of sites: {}".format(str(len(Sites))))
    # print("Number of genes: {}".format(str(len(Genes))))
    # print("Sites with counts: {}".format(str(len(Counts))))
    # print("Genes with MK: {}".format(str(len(MK))))

    return MK, Genes


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


def calculate_statistic(mk, test, pseudocount=0):
    """Takes an MK object and returns a dictionary
    with the statistics asked"""

    tests = dict()
    # Calculate neutrality index
    if 'NI' in test:
        tests['NI.pval'] = float('nan')
        try:
            tests['NI'] = mk.neutrality_index(log=True,
                                              pseudocount=pseudocount)
        except ZeroDivisionError:
            tests['NI'] = float('nan')

    # Calculate ratio
    if 'ratio' in test:
        tests['ratio.pval'] = float('nan')
        try:
            tests['ratio'] = mk.mk_ratio(pseudocount=pseudocount)
        except ZeroDivisionError:
            tests['ratio'] = float('nan')

    # Hypergeometric test
    if 'hg' in test:
        tests['hg'], tests['hg.pval'] = mk.hg_test(pseudocount=pseudocount)

    # G test of indenpendece try multiple corrections
    if 'G' in test:
        try:
            g, pval, df, E = mk.g_test(correction='none',
                                       pseudocount=pseudocount)
            tests['G'] = g
            tests['G.pval'] = pval
            tests['G.df'] = df
            tests['G.E'] = E
        except ValueError:
            tests['G'] = float('nan')
            tests['G.pval'] = float('nan')
            tests['G.df'] = float('nan')
            tests['G.E'] = float('nan')

    if 'G_Yates' in test:
        try:
            g, pval, df, E = mk.g_test(correction='yates',
                                       pseudocount=pseudocount)
            tests['G_Yates'] = g
            tests['G_Yates.pval'] = pval
            tests['G_Yates.df'] = df
            tests['G_Yates.E'] = E
        except ValueError:
            tests['G_Yates'] = float('nan')
            tests['G_Yates.pval'] = float('nan')
            tests['G_Yates.df'] = float('nan')
            tests['G_Yates.E'] = float('nan')

    if 'G_Williams' in test:
        try:
            g, pval, df, E = mk.g_test(correction='williams',
                                       pseudocount=pseudocount)
            tests['G_Williams'] = g
            tests['G_Williams.pval'] = pval
            tests['G_Williams.df'] = df
            tests['G_Williams.E'] = E
        except ValueError:
            tests['G_Williams'] = float('nan')
            tests['G_Williams.pval'] = float('nan')
            tests['G_Williams.df'] = float('nan')
            tests['G_Williams.E'] = float('nan')

    # Eyre-Walker alpha
    if 'alpha' in test:
        tests['alpha.pval'] = float('nan')
        try:
            tests['alpha'] = mk.alpha(pseudocount=pseudocount)
        except ZeroDivisionError:
            tests['alpha'] = float('nan')

    if 'DoS' in test:
        tests['DoS.pval'] = float('nan')
        try:
            tests['DoS'] = mk.DoS(pseudocount=pseudocount)
        except ZeroDivisionError:
            tests['DoS'] = float('nan')

    return tests


def confirm_midas_merge_files(args):
    """Confirm files are present. No integrity check"""

    # Check files exist in input directory
    file_list = os.listdir(args.indir)
    if 'snps_freq.txt' not in file_list:
        msg = "Could not find snps_freq.txt at {}".format(args.indir)
        raise FileNotFoundError(msg)
    if 'snps_info.txt' not in file_list:
        msg = "Could not find snps_info.txt at {}".format(args.indir)
        raise FileNotFoundError(msg)
    if 'snps_depth.txt' not in file_list:
        msg = "Could not find snps_depth.txt at {}".format(args.indir)
        raise FileNotFoundError(msg)
    if not os.path.isfile(args.metadata_file):
        msg = "Could not find metadata file {}".format(args.metadata_file)
        raise FileNotFoundError(msg)

    print("\tAll files found")
    return


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
        if( (site_crosstab.diagonal() == [0,0]).all() or (np.fliplr(site_crosstab).diagonal() == [0,0]).all()):
            mutation_type = 'fixed'
        else:
            mutation_type = 'polymorphic'

        Dist.append(mutation_type)
        # print(info.site_id.iloc[i], mutation_type)

    return(Dist)


def mktest_fisher_exact(g):
    """Per perform the fisher's exact test on a gene MK
    contingency table."""

    tab = np.array([[g.Dn, g.Pn], [g.Ds, g.Ps]])
    oddsratio, pval = stats.fisher_exact(tab, alternative='two-sided')
    # oddsratio, pval = stats.fisher_exact(tab, alternative='greater')

    return pval


def process_arguments():
    # Read arguments
    parser_format = argparse.ArgumentDefaultsHelpFormatter
    parser = argparse.ArgumentParser(formatter_class=parser_format)
    required = parser.add_argument_group("Required arguments")

    # Define description
    parser.description = ("Script to perform McDonald-Kreitman test "
                          "from a midas output file of snps.")

    # Define required arguments
    required.add_argument("--indir", help="Input directory",
                          type=str,
                          required=True)
    required.add_argument("--metadata_file", help="Mapping file for samples",
                          type=str,
                          required=True)
    required.add_argument("--group1", help="Group1 of comparison",
                          type=str,
                          required=True)
    required.add_argument("--group2", help="Group2 of comparison",
                          type=str,
                          required=True)

    # Define other arguments
    parser.add_argument("--functions", help=("Which set of functions to use, "
                                             "pandas or classes"),
                        default='pandas', type=str,
                        choices=['pandas', 'classes'])
    parser.add_argument("--min_count", help=("Min number of reads (depth) "
                                             " at a position in a sample to "
                                             "consider that sample for "
                                             "that position."),
                        default=1, type=int)
    parser.add_argument("--min_cov", help=("min metagenomic coverage across "
                                           "genome in a sample to keep that "
                                           "sample for that genome"),
                        default=1.0, type=float)
    parser.add_argument("--nrows", help="Number of gene positions to read",
                        default=float('inf'), type=float)
    parser.add_argument("--outfile", help="Output file with results",
                        default="mk_results.txt", type=str)
    parser.add_argument("--permutations", help=("Number of permutations to "
                                                "perform to establish "
                                                "significance"),
                        type=int, default=0)
    parser.add_argument("--pseudocount", help=("Pseudocount value to use "
                                               "in contingency tables"),
                        default=0, type=int)
    parser.add_argument("--seed", help="Permutation seed",
                        type=int, default=None)
    parser.add_argument("--test", help=("Eventually specify test to perform."
                                        "all performs all tests. G performs "
                                        "a G test without correction."
                                        "G_Yates performs a G test with the "
                                        "Yates correction. G_Williams "
                                        "performs a G test with the Williams "
                                        "correction. hg performs the "
                                        "hypergeometric (Fisher's Exact) "
                                        "test. NI rerturns the neutrality "
                                        "index. alpha returs the Eyre-"
                                        "Walker alpha. Ratio returns the MK "
                                        "rati. DoS is the direction of "
                                        "selection statistic."),
                        default="hg", type=str,
                        choices=['all', 'G', 'G_Yates', 'G_Williamps',
                                 'hg', 'NI', 'alpha', 'ratio', 'DoS'])

    # Read arguments
    print("Reading arguments")
    args = parser.parse_args()

    # Processing goes here if needed
    if args.seed is None and args.permutations > 0:
        args.seed = np.random.randint(1000)*2 + 1

    return args


def process_metadata_file(mapfile, permute=False):
    """Process metadata file and permute if needed.
    Creates dictionaries of ID => group and grop => IDs.
    Uses pandas."""

    map = pd.read_csv(mapfile, sep='\t')
    map.index = map.ID
    # print(map.head())

    # if permute:
    #     map['Group'] = np.random.permutation(map.Group)
    #
    # # Create dictionary sampleID => group
    # Samples = {map.ID[i]: [map.Group[i]] for i in range(len(map))}
    #
    # # Create dictionary group => sampleIDs
    # Groups = dict()
    # for g in set(map.Group):
    #     samples = list(map.ID[map.Group == g])
    #     Groups[g] = samples

    # NOTE: ADD CHECKS FOR COLUMN NAMES AND VALUES

    return map


def process_snps_depth_file(args, Map, Sites):
    """Use depth to decide which samples to keep for each site.

    It modifies Sites by removing sites that do not have samples at
    minimum depth for both groups. It returns the dictionary Counts,
    that provides presence absence vector per site indicating which
    samples passed the threshold for each site.

    Might not need indices anymore."""

    Counts = {}
    with open(args.indir + '/snps_depth.txt') as depth_fh:
        # Read header
        header = depth_fh.readline()
        header = header.rstrip()
        header = header.split('\t')

        # Get samples
        samples = np.array(header[1:])

        # Match map to samples in file
        Map_present = Map.loc[samples]

        # Create index of for samples to keep based on Map
        s_ii = Map_present.Group.notnull()

        # Read all lines after header
        depth_reader = csv.reader(depth_fh, delimiter='\t')
        i = 0
        for row in depth_reader:
            i += 1
            # Break if max number of rows passed.
            if i > args.nrows:
                break

            # Get site ID and check if it is in Sites (for MK this is
            # equivalent to check if this a gene)
            site_id = row[0]
            if not (site_id in Sites):
                continue

            # Get all counts, convert to integer and decide which are present
            # based on min depth threshold
            counts = row[1:]
            counts = np.array(counts, dtype='int')
            counts = counts >= args.min_count

            # Get counts per group
            samples1 = (Map_present.Group[counts & s_ii] == args.group1).sum()
            samples2 = (Map_present.Group[counts & s_ii] == args.group2).sum()

            # Keep sites with at least one sample per group
            if not (samples1 > 1 and samples2 > 1):
                # Delete sites that don't pass the tresshold
                if site_id in Sites:
                    del Sites[site_id]
            else:
                # If site is kept, create dictionary entry with index of
                # samples tothat pass threshold
                # NOTE: ASSUMING SAME ORDER IN SAMPLES BETWEEN SITES
                Counts[site_id] = counts

    depth_fh.close()

    return Counts


def process_snp_freq_file(args, Counts, Map, Sites):
    """Process snp_freq.txt from MIDAS. Produces MK table"""

    print("Processing snp_freq.txt")
    MK = {}
    with open(args.indir + '/snps_freq.txt') as freqs_fh:
        # Read header
        header = freqs_fh.readline()
        header = header.rstrip()
        header = header.split('\t')

        # Get samples
        samples = np.array(header[1:])

        # Match map to samples in file
        Map_present = Map.loc[samples]

        # Create index of for samples to keep based on Map
        s_ii = Map_present.Group.notnull()

        # Read all lines after header
        freqs_reader = csv.reader(freqs_fh, delimiter='\t')
        i = 0
        for row in freqs_reader:
            i += 1
            # Break if max number of rows passed.
            if i > args.nrows:
                break

            # Skip is site not selected (because not in gene
            # or not enough depth)
            site_id = row[0]
            if not (site_id in Sites):
                continue

            # Get info from GenomeSite and Gene objects
            gene = Sites[site_id].gene_id
            s_type = Sites[site_id].substitution_type()
            present_index = Counts[site_id]

            # Initiallize MKtest object if SNP is in novel gene
            if not (gene in MK):
                MK[gene] = MKtest(name=gene)

            # find allele per sample
            allele = np.array(row[1:], dtype='float') < 0.5

            # Count major allele on each group
            group1_count = (present_index &
                            s_ii &
                            allele &
                            (Map_present.Group == args.group1)).sum()
            group2_count = (present_index &
                            s_ii &
                            allele &
                            (Map_present.Group == args.group2)).sum()

            # Classify variants based on distribution
            if group1_count > 0 and group2_count > 0:
                fixed = False
            elif group1_count == 0 and group2_count == 0:
                fixed = False
            else:
                fixed = True

            # Classify variants based on effect on aminoacid
            if s_type == 'synonymous':
                if fixed:
                    MK[gene].update(Ds=1)
                else:
                    MK[gene].update(Ps=1)
            elif s_type == 'non-synonymous':
                if fixed:
                    MK[gene].update(Dn=1)
                else:
                    MK[gene].update(Pn=1)
            else:
                raise ValueError("Invalid substitution type")

            # print("\t".join([gene, s_type, str(fixed), site_id]))
    freqs_fh.close()

    return MK


def process_snp_info_file(args):
    """Process the snps_info.txt file from MIDAS. Assumes column
    info based on hard-coded values. Developped for MIDAS 1.3.1.abs

    Returns dictionaries of Genes (class Gene) and
    Sites (class GenineSite)"""

    Genes = {}
    Sites = {}
    with open(args.indir + '/snps_info.txt') as info_fh:
        # Get headers
        header = info_fh.readline()
        header = header.split('\t')

        # Create reader for following lines
        info_reader = csv.reader(info_fh, delimiter='\t')
        i = 0

        # Set columns. These are hard-coded values that specify which
        # information is on each column. Based on MIDAS 1.3.1
        site_id_col = 0
        contig_col = 1
        pos_col = 2
        ref_allele_col = 3
        major_allele_col = 4
        minor_allele_col = 5
        locus_type_col = 11
        gene_id_col = 12
        aminoacids_col = 15

        # print("============HEADERs============")
        # print(">Site id: {}".format(header[site_id_col]))
        # print(">Contig: {}".format(header[contig_col]))
        # print(">Position: {}".format(header[pos_col]))
        # print(">Ref allele: {}".format(header[ref_allele_col]))
        # print(">Major allele: {}".format(header[major_allele_col]))
        # print(">Minor allele: {}".format(header[minor_allele_col]))
        # print(">Locus type: {}".format(header[locus_type_col]))
        # print(">Gene id: {}".format(header[gene_id_col]))
        # print(">Aminoacids: {}".format(header[aminoacids_col]))

        # Read each line after header
        for row in info_reader:
            i += 1

            # Break if max number of rows given
            if i > args.nrows:
                break

            # Get info
            gene = row[gene_id_col]
            site_id = row[site_id_col]
            aminoacids = row[aminoacids_col]

            # Skip intergenig regions
            if gene == 'NA':
                continue

            # Get aminoacid per variant
            aa = aminoacids.split(',')
            # print(site_id)

            # Define GenomeSite object
            Sites[site_id] = GenomeSite(site_id=site_id,
                                        contig=row[contig_col],
                                        position=row[pos_col],
                                        ref_allele=row[ref_allele_col],
                                        major_allele=row[major_allele_col],
                                        minor_allele=row[minor_allele_col],
                                        locus_type=row[locus_type_col],
                                        gene_id=gene,
                                        aminoacid_A=aa[0],
                                        aminoacid_C=aa[1],
                                        aminoacid_G=aa[2],
                                        aminoacid_T=aa[3])

            # Define Gene object
            # Check if site is in previously defined gene
            if gene in Genes:
                # Update gene if already present by extending position.
                Genes[gene].extend(row[pos_col])

            else:
                # Create new Gene object
                Genes[gene] = Gene(gene_id=gene,
                                   contig=row[contig_col],
                                   start=row[pos_col],
                                   end=row[pos_col])

    info_fh.close()

    return Genes, Sites


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


def test_and_write_results(MK, Genes, outfile,
                           test='hg', pseudocount=0,
                           permutations=0):
    """Take MK results, perform test and write outfile with
    results."""

    # Get list of tests to perform
    supported_tests = ['NI', 'ratio', 'hg',
                       'G', 'G_Yates', 'G_Williams',
                       'alpha', 'DoS']
    if test == 'all':
        test = supported_tests
    elif test not in supported_tests:
        raise ValueError("Test not supported")
    else:
        test = [test]

    # Create header
    header_base = ['gene', 'contig', 'start', 'end',
                   'Dn', 'Ds', 'Pn', 'Ps']
    pval_list = [''.join([t, '.pval']) for t in test]
    if permutations > 0:
        perm_list = [''.join([t, '.perm']) for t in test]
    else:
        perm_list = []

    header = header_base + test + pval_list + perm_list

    # Open files for output
    with open(outfile, mode='w') as fh:
        # Write header as first line of results
        fh.write("\t".join(header) + "\n")

        # Iterate over every MK element
        # print(MK[0])
        for gene, mk in MK[0].items():
            # print(gene)
            # print(mk.Dn, mk.Ds, mk.Pn, mk.Ps)

            if permutations == 0:
                # Calculate statistics
                tests = calculate_statistic(mk, test, pseudocount)

                # prepare res
                res = [str(tests[t]) for t in test + pval_list]
                res = [gene, Genes[gene].contig,
                       str(Genes[gene].start),
                       str(Genes[gene].end),
                       str(mk.Dn), str(mk.Ds),
                       str(mk.Pn), str(mk.Ps)] + res

                # res = [gene, Genes[gene].contig, str(Genes[gene].start),
                #        str(Genes[gene].end),
                #        str(mk.Dn), str(mk.Ds), str(mk.Pn), str(mk.Ps),
                #        str(ni), str(ratio), str(ratio_pseudo),
                #        str(hg_odds), str(hg_p), str(hg_odds_pseudo),
                #        str(hg_p_pseudo),
                #        str(g_none_p), str(g_yates_p),str(g_williams_p),
                #        str(g_none_p_pseudo), str(g_yates_p_pseudo),
                #        str(g_williams_p_pseudo),
                #        str(alpha), str(alpha_pseudo)]

                # th.write(str(res) + "\n")
                fh.write("\t".join(res) + "\n")
                # alpha = mk.alpha()
                # print("MK ratio is: {}".format(str(ratio)))
                # print("MK alpha is: {}".format(str(alpha)))
            elif permutations > 0:
                res = test_by_permutation(gene, MK, permutations,
                                          test, pval_list, pseudocount)
                # print(res)
                res = [str(res[k]) for k in res]
                res = [gene, Genes[gene].contig,
                       str(Genes[gene].start),
                       str(Genes[gene].end),
                       str(mk.Dn), str(mk.Ds),
                       str(mk.Pn), str(mk.Ps)] + res
                fh.write("\t".join(res) + "\n")
            else:
                raise ValueError("Invalid permutations")

    fh.close()
    # th.close()


def test_by_permutation(gene, MK, permutations, test, pval_list, pseudocount):
    nperm = int(permutations + 1)
    perm_table = np.full(shape=(nperm, len(test)), fill_value=np.nan)
    row = 0
    for p in MK:
        if gene in p:
            p_stat = calculate_statistic(p[gene], test,
                                         pseudocount=pseudocount)
            p_res = [p_stat[t] for t in test]
            perm_table[row] = p_res

        row = row + 1

    # Pvalues
    nperms = nperm - np.isnan(perm_table).sum(axis=0)
    np.seterr(invalid='ignore', divide='ignore')
    perm_pvals = (perm_table >= perm_table[0]).sum(axis=0) / nperms
    np.seterr(invalid='raise', divide='raise')
    nperm_names = [''.join([t, '.nperm']) for t in test]

    keys = np.concatenate((test, pval_list, nperm_names))
    vals = np.concatenate((perm_table[0], perm_pvals, nperms))
    # vals = np.array(vals, dtype=np.character)
    res = dict(zip(keys, vals))

    return res


if __name__ == "__main__":
    args = process_arguments()

    # Check midas files
    print("Checking MIDAS files exist")
    confirm_midas_merge_files(args)

    if args.functions == 'classes':
        # Create dictionaries that have all the samples per group (Groups),
        # and the group to which each sample belongs (Samples)
        # Read DataFrame
        # NOTE: need to check names
        print("Read metadata")
        Map = process_metadata_file(args.metadata_file)

        print("Calculate MK contingency tables")
        MK, Genes = calculate_contingency_tables(Map, args)
        MK = [MK]
        if args.permutations > 0:
            print("Permuting")
            print("Seed is {}".format(str(args.seed)))
            np.random.seed(args.seed)
            for i in range(args.permutations):
                Sp, Gp = process_metadata_file(args.metadata_file,
                                               permute=True)
                mk, genes = calculate_contingency_tables(Sp, Gp, args)
                MK.append(mk)

        print("Testing and writing")
        test_and_write_results(MK, Genes, args.outfile,
                               test=args.test,
                               pseudocount=args.pseudocount,
                               permutations=args.permutations)
    elif args.functions == 'pandas':
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

        # Write results
        Genes.to_csv(args.outfile, sep="\t", na_rep='NA', index=False)

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
import sutilspy
import csv
import numpy as np
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
        if base in ['A','a']:
            return(self.aminoA)
        elif base in ['C','c']:
            return(self.aminoC)
        elif base in ['G','g']:
            return(self.aminoG)
        elif base in ['T','t']:
            return(self.aminoT)
        else:
            raise ValueError("base must be one of the four canonical nucleoties")

    def substitution_type(self):
        substitution_type = ''
        if self.codon_aminoacid(base = self.major_allele) == self.codon_aminoacid(base = self.minor_allele):
            substitution_type = 'synonymous'
        else:
            substitution_type = 'non-synonymous'

        return(substitution_type)

class Gene:
    """A class for representing a gene"""

    def __init__(self, gene_id,contig,start,end, strand = ''):
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

    def __init__(self, name, Ds = 0, Dn = 0, Ps = 0, Pn = 0):
        self.name = name
        self.Dn = Dn
        self.Ds = Ds
        self.Ps = Ps
        self.Pn = Pn

    def update(self, Ds = 0, Dn = 0, Ps = 0, Pn = 0):
        """Update the contigency matrix"""
        self.Dn += Dn
        self.Ds += Ds
        self.Ps += Ps
        self.Pn += Pn

    def mk_ratio(self, pseudocount = 0):
        """Calculate the McDonald Kreitman ratio (Dn/Ds)/(Pn/Ps)"""
        ratio = ((self.Dn + pseudocount) / (self.Ds + pseudocount)) / ((self.Pn + pseudocount) / (self.Ps + pseudocount))
        return(ratio)

    def alpha(self, pseudocount = 0):
        """Calculate the Smith & Eyre-Walker alpha 1 - """
        ni = self.neutrality_index(pseudocount = pseudocount, log = False)
        alpha = 1 - ni
        return(alpha)

    def hg_test(self, pseudocount = 0):
        """Hypergeometric (Fisher's exact) test"""

        res = stats.fisher_exact([[self.Ds + pseudocount,self.Ps + pseudocount],
                                  [self.Dn + pseudocount,self.Pn + pseudocount]])
        return(res)

    def g_test(self, correction, pseudocount = 0):
        """G-test for independence. Original McDonald & Kreitman 1991 suggestion"""

        # Create 2x2 contingency matrix
        mat = np.matrix([[self.Ds + pseudocount,self.Ps + pseudocount],
                         [self.Dn + pseudocount,self.Pn + pseudocount]])

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
            # it doesn't make much difference (http://www.biostathandbook.com/small.html)
            g, p, df, e = stats.chi2_contingency(observed=mat,
                                         lambda_="log-likelihood",
                                         correction=False)

            # Calculate q correction. Only for 2 x 2 table
            n = mat.sum()
            q = 1 + (n * (1 / mat.sum(axis = 1)).sum() - 1) * (n * (1 / mat.sum(axis = 0)).sum() - 1) / (6 * n)

            # correct g and recalculate p-value
            g = g / q
            p = 1 - stats.chi2.cdf(g, df)

            # combine results
            res = [g, p , df, e]

        else:
            raise ValueError("Correction must be one of 'none', 'yates' or 'williams'")

        return(res)
    def neutrality_index(self, pseudocount = 1, log = True):
        """Calculate neutrality index (Pn/Dn)/(Ps/Ds). Following Li et al. (2008), we add a psedocount and return the -log10(NI)"""

        ni = ((self.Pn + pseudocount) / (self.Dn + pseudocount)) / ((self.Ps + pseudocount) / (self.Ds + pseudocount))

        if log:
            ni = -np.log10(ni)


        return(ni)

def process_snp_info_file(args):
    """Process the snps_info.txt file from MIDAS"""

    Genes = {}
    Sites = {}
    with open(args.indir + '/snps_info.txt') as info_fh:
        header = info_fh.readline()
        header = header.split('\t')
        print(header)
        info_reader = csv.reader(info_fh, delimiter = '\t')
        i = 0

        # Set columns
        site_id_col = 0
        contig_col = 1
        pos_col = 2
        ref_allele_col = 3
        major_allele_col = 4
        minor_allele_col = 5
        locus_type_col = 11
        gene_id_col = 12
        aminoacids_col = 15

        print("============HEADERs============")
        print(">Site id: {}".format(header[site_id_col]))
        print(">Contig: {}".format(header[contig_col]))
        print(">Position: {}".format(header[pos_col]))
        print(">Ref allele: {}".format(header[ref_allele_col]))
        print(">Major allele: {}".format(header[major_allele_col]))
        print(">Minor allele: {}".format(header[minor_allele_col]))
        print(">Locus type: {}".format(header[locus_type_col]))
        print(">Gene id: {}".format(header[gene_id_col]))
        print(">Aminoacids: {}".format(header[aminoacids_col]))

        #
        for row in info_reader:
            i += 1
            if i > args.nrows:
                break
            #print(row)
            #print(row[gene_id_col], row[site_id_col])
            #print(row[aminoacids_col])
            gene = row[gene_id_col]
            site_id = row[site_id_col]
            aminoacids = row[aminoacids_col]
            #print(aminoacids)
            #print(site_id)

            if gene == 'NA':
                # skip intergenig regions
                continue

            #print("\tgene")
            # Get aminoacid per position
            aa = aminoacids.split(',')
            #print(aa)

            # Define site
            #print(site_id)
            Sites[site_id] = GenomeSite(site_id = site_id,
                                        contig = row[contig_col],
                                        position = row[pos_col],
                                        ref_allele = row[ref_allele_col],
                                        major_allele = row[major_allele_col],
                                        minor_allele = row[minor_allele_col],
                                        locus_type = row[locus_type_col],
                                        gene_id = gene, aminoacid_A = aa[0],
                                        aminoacid_C = aa[1],
                                        aminoacid_G = aa[2],
                                        aminoacid_T = aa[3])

            # For genes
            if gene in Genes:
                # update genes
                Genes[gene].extend(row[pos_col])
                #print(gene)
                #print(Genes[gene])
                #Genes[gene].info()

            else:
                # Define gene
                Genes[gene] = Gene(gene_id=gene, contig = row[contig_col],
                                   start = row[pos_col], end = row[pos_col])
                #Genes[gene].info()
                #print(Genes[gene])


    info_fh.close()
    #print(Groups)
    print("Number of sites: {}".format(str(len(Sites))))
    print("Number of genes: {}".format(str(len(Genes))))

    return Genes, Sites

def process_snps_depth_file(args,Groups,Sites):
    """Use depth to decide which samples to keep. It modifies Sites and returns Counts"""

    Counts = {}
    with open(args.indir + '/snps_depth.txt') as depth_fh:
        header = depth_fh.readline()
        header = header.rstrip()
        header = header.split('\t')

        # Get sample and column indices
        samples = header[1:]
        indices = {}
        for s in samples:
            indices[s] = header.index(s)
        print(indices)

        depth_reader = csv.reader(depth_fh, delimiter = '\t')
        i = 0
        for row in depth_reader:
            i += 1
            if i > args.nrows:
                break
            #print(row)

            # Get site ID and check if it is in Sites (for MK this is
            # equivalent to check if this a gene)
            site_id = row[0]
            #print(site_id)
            if not site_id in Sites:
                continue

            # Get all counts and convert to integer
            counts = row[1:]
            counts = list(map(int,counts))
            #print(counts)

            # Convert count to presence/absence vector based on
            # threshold of number of counts to use position in sample
            counts = [int(c >= args.min_count) for c in counts]

            # Get counts per group
            # GLITCH: Here it fails if map has extra samples not present in files
            #print(set(Groups[args.group1]) & set(indices.keys()))
            #print(args.group1)
            #print(Groups[args.group1])
            #print(indices.keys())
            #print(set(indices.keys()))

            samples1 = [int(counts[ indices[l] - 1 ]) for l in set(Groups[args.group1]) & set(indices.keys())]
            samples2 = [int(counts[ indices[l] - 1 ]) for l in set(Groups[args.group2]) & set(indices.keys())]
            samples1 = sum(samples1)
            samples2 = sum(samples2)
            #print(samples1)
            #print(samples2)
            if not (samples1 > 1 and samples2 > 1):
                #print("\t====Group1:{},Group2:{},SiteID:{}====".format(samples1,samples2,site_id))
                # delete
                #print(site_id)
                if site_id in Sites:
                    del Sites[site_id]
            else:
                # NOTE: ASSUMING SAME ORDER IN SAMPLES BETWEEN SITES
                Counts[site_id] = counts



    depth_fh.close()
    print("Number of sites: {}".format(str(len(Sites))))
    print("Number of genes: {}".format(str(len(Genes))))
    print("Sites with counts: {}".format(str(len(Counts))))

    return Counts

def process_snp_freq_file(args,Counts,Groups,Samples):
    """Process snp_freq.txt from MIDAS. Produces MK table"""

    print(Groups)
    MK = {}
    with open(args.indir + '/snps_freq.txt') as freqs_fh:
        header = freqs_fh.readline()
        header = header.rstrip()
        header = header.split('\t')

        # Get sample and column indices
        samples = header[1:]
        indices = {}
        for s in samples:
            indices[s] = header.index(s)
        print(indices)
        print(header)

        freqs_reader = csv.reader(freqs_fh, delimiter = '\t')
        i = 0
        for row in freqs_reader:
            i += 1
            if i > args.nrows:
                break

            # Check if site was selected based on sites
            site_id = row[0]
            if not site_id in Sites:
                #print("==Skipping")
                continue

            gene = Sites[site_id].gene_id
            s_type = Sites[site_id].substitution_type()
            present_index = np.array(Counts[site_id])
            group_index = np.array([Samples[s][0] for s in samples])
    #         if site_id == '77719':
    #             print("==========================")
    #             print(row)
    #             print(site_id)
    #             print("Major Allele: {}".format(Sites[site_id].major_allele))
    #             print("Minor Allele: {}".format(Sites[site_id].minor_allele))
    #             print("Substitution type: {}".format(s_type))
    #             print("Gene: {}".format(gene))
    #             print(present_index)
    #             print(group_index)

            # Create MKtest if needed
            if gene not in MK:
                MK[gene]= MKtest(name=gene)

            # find allele per sample
            allele_freqs = np.array([int(float(f) < 0.5) for f in row[1:]])
            #print(allele_freqs)

            # Remove non covered positions
            ii = np.where(present_index)
            group_index = group_index[ii]
            allele_freqs = allele_freqs[ii]
            #print(group_index)
            #print(allele_freqs)

            # Count alleles per group
            group1_count = allele_freqs[np.where(group_index == args.group1)].sum()
            group2_count = allele_freqs[np.where(group_index == args.group2)].sum()
            #print(group1_count)
            #print(group2_count)

            if group1_count > 0 and group2_count > 0:
                fixed = False
            elif group1_count > 0 or group2_count > 0:
                fixed = True

            if s_type == 'synonymous':
                if fixed:
                    MK[gene].update(Ds = 1)
                else:
                    MK[gene].update(Ps = 1)
            elif s_type == 'non-synonymous':
                if fixed:
                    MK[gene].update(Dn = 1)
                else:
                    MK[gene].update(Pn = 1)
            else:
                raise ValueError("Invalid substitution type")

            #print("==========================")

    freqs_fh.close()
    print("Number of sites: {}".format(str(len(Sites))))
    print("Number of genes: {}".format(str(len(Genes))))
    print("Sites with counts: {}".format(str(len(Counts))))
    print("Genes with MK: {}".format(str(len(MK))))

    return MK

def confirm_midas_merge_files(args):
    """Confirm files are present. No integrity check"""

    # Check files exist in input directory
    file_list = os.listdir(args.indir)
    if 'snps_freq.txt' not in file_list:
        raise FileNotFoundError("Could not find snps_freq.txt at {}".format(args.indir))
    if 'snps_info.txt' not in file_list:
        raise FileNotFoundError("Could not find snps_info.txt at {}".format(args.indir))
    if 'snps_depth.txt' not in file_list:
        raise FileNotFoundError("Could not find snps_depth.txt at {}".format(args.indir))
    if not os.path.isfile(args.metadata_file):
        raise FileNotFoundError("Could not find metadata file {}".format(args.metadata_file))

if __name__ == "__main__":

    # Argparse
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    required = parser.add_argument_group("Required arguments")
    required.add_argument("--indir", help = "Input directory", type = str,
                          required = True)
    required.add_argument("--metadata_file", help = "Mapping file for samples", type = str,
                          required = True)
    required.add_argument("--group1", help = "Group1 of comparison", type = str,
                          required = True)
    required.add_argument("--group2", help = "Group2 of comparison",
                          required = True)

    parser.add_argument("--test", help = "Eventually specify test to perform",
                        default = "G", type = str)
    parser.add_argument("--outfile", help = "Output file with results",
                   default = "mk_results.txt", type = str)
    parser.add_argument("--min_count", help = "min depth at a position in a sample to consider that sample in that position",
                        default = 1, type = int)
    parser.add_argument("--nrows", help = "Number of gene positions to read",
                        default = float('inf'), type = float)
    parser.add_argument("--tables", help = "Output file for contingency tables",
                        default = "mk_tables.txt", type = str)
    parser.add_argument("--pseudocount", help = "Pseudocount value to use in contingency tables",
                        default = 1, type = int)

    args = parser.parse_args()

    ######## Check files #################
    confirm_midas_merge_files(args)

    #### Read metadata ####
    Groups = sutilspy.io.process_run_list(args.metadata_file,
                                          1, 0, header = True)
    Samples = sutilspy.io.process_run_list(args.metadata_file,
                                           0, 1, header = True)

    ######## Read info #######
    Genes, Sites = process_snp_info_file(args)

    ###### Chose sites based on depth in groups to compare #######
    Counts = process_snps_depth_file(args, Groups, Sites)

    ####### Read frequencies and calculate #########
    MK = process_snp_freq_file(args, Counts, Groups, Samples)

    ################ Test and results ########
    with open(args.outfile,mode='w') as fh, open(args.tables,mode='w') as th:
        header = ['gene','contig','start','end',
                  'Dn','Ds','Pn','Ps',
                  'ni', 'ratio','ratio_pseudo','hg_odds','hg_p','hg_odds_pseudo','hg_p_pseudo',
                  'g_none_p','g_yates_p','g_williams_p',
                  'g_none_p_pseudo','g_yates_p_pseudo','g_williams_p_pseudo',
                  'alpha','alpha_pseudo']
        fh.write("\t".join(header) + "\n")
        for gene,mk in MK.items():
            th.write("=============================================\n")
            th.write(gene)
            th.write("\t\tFixed\tPolymorphic\n\tSynonymous\t{}\t{}\n\tnon-synonymous\t{}\t{}\n".format(mk.Ds,mk.Ps,mk.Dn,mk.Pn))

            # Calculate neutrality index
            try:
                ni = mk.neutrality_index(log=True, pseudocount = args.pseudocount)
            except ZeroDivisionError:
                ni = float('nan')

            # Calculate ratio with and without pseudocount
            try:
                ratio = mk.mk_ratio(pseudocount=0)
            except ZeroDivisionError:
                ratio = float('nan')
            try:
                ratio_pseudo = mk.mk_ratio(pseudocount=args.pseudocount)
            except ZeroDivisionError:
                ratio = float('nan')

            # Hypergeometric test
            hg_odds, hg_p = mk.hg_test(pseudocount = 0)
            hg_odds_pseudo, hg_p_pseudo = mk.hg_test(pseudocount = args.pseudocount)

            # G test of indenpendece try multiple corrections
            try:
                g_none, g_none_p, g_none_df, g_none_E = mk.g_test(correction='none',
                                                                  pseudocount=0)
            except ValueError:
                g_none = float('nan')
                g_none_p = float('nan')
                g_none_df = float('nan')
                g_none_E = float('nan')

            try:
                g_yates, g_yates_p, g_yates_df, g_yates_E = mk.g_test(correction='yates',
                                                                      pseudocount=0)
            except ValueError:
                g_yates = float('nan')
                g_yates_p = float('nan')
                g_yates_df = float('nan')
                g_yates_E = float('nan')

            try:
                g_williams, g_williams_p, g_williams_df, g_williams_E = mk.g_test(correction='williams',
                                                                                  pseudocount=0)
            except ValueError:
                g_williams = float('nan')
                g_williams_p = float('nan')
                g_williams_df = float('nan')
                g_williams_E = float('nan')

            # G test for independence with pseududocounts
            try:
                g_none_pseudo, g_none_p_pseudo, g_none_df_pseudo, g_none_E_pseudo = mk.g_test(correction='none',
                                                                                              pseudocount=args.pseudocount)
            except ValueError:
                g_none_pseudo = float('nan')
                g_none_p_pseudo = float('nan')
                g_none_df_pseudo = float('nan')
                g_none_E_pseudo = float('nan')

            try:
                g_yates_pseudo, g_yates_p_pseudo, g_yates_df_pseudo, g_yates_E_pseudo = mk.g_test(correction='yates',
                                                                                                  pseudocount=args.pseudocount)
            except ValueError:
                g_yates_pseudo = float('nan')
                g_yates_p_pseudo = float('nan')
                g_yates_df_pseudo = float('nan')
                g_yates_E_pseudo = float('nan')

            try:
                g_williams_pseudo, g_williams_p_pseudo, g_williams_df_pseudo, g_williams_E_pseudo = mk.g_test(correction='williams',
                                                                                                              pseudocount=args.pseudocount)
            except ValueError:
                g_williams_pseudo = float('nan')
                g_williams_p_pseudo = float('nan')
                g_williams_df_pseudo = float('nan')
                g_williams_E_pseudo = float('nan')


            # Eyre-Walker alpha
            try:
                alpha = mk.alpha(pseudocount=0)
            except ZeroDivisionError:
                alpha = float('nan')

            alpha_pseudo = mk.alpha(pseudocount=args.pseudocount)

            # prepare res
            res = [gene, Genes[gene].contig, str(Genes[gene].start), str(Genes[gene].end),
                   str(mk.Dn), str(mk.Ds), str(mk.Pn), str(mk.Ps),
                   str(ni), str(ratio), str(ratio_pseudo),
                   str(hg_odds), str(hg_p), str(hg_odds_pseudo),str(hg_p_pseudo),
                   str(g_none_p), str(g_yates_p),str(g_williams_p),
                   str(g_none_p_pseudo), str(g_yates_p_pseudo),str(g_williams_p_pseudo),
                   str(alpha), str(alpha_pseudo)]

            th.write(str(res) + "\n")
            fh.write("\t".join(res) + "\n")
            #alpha = mk.alpha()
            #print("MK ratio is: {}".format(str(ratio)))
            #print("MK alpha is: {}".format(str(alpha)))
    fh.close()
    th.close()

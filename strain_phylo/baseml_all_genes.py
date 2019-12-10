#!/usr/bin/env python
# Copyright (C) 2019 Sur Herrera Paredes

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


from Bio import AlignIO
from skbio import TreeNode
from Bio import Align
from Bio.Phylo.PAML import baseml
import os
import io
import pandas as pd
import numpy as np
import argparse
import multiprocessing as mp


def process_arguments():
    # Read arguments
    parser_format = argparse.ArgumentDefaultsHelpFormatter
    parser = argparse.ArgumentParser(formatter_class=parser_format)
    required = parser.add_argument_group("Required arguments")

    # Define description
    parser.description = ("Script that uses baseml to generate gene trees "
                          "for all gene alignments in an input")

    # Define required arguments
    required.add_argument("--aln_dir", help=("Path to directory containing "
                                             "alignments in fasta format"),
                          required=True, type=str)
    required.add_argument("--cov_file",
                          help=("A tab-delimited file having gene coverage "
                                "per strain."),
                          required=True, type=str)
    required.add_argument("--master_tree",
                          help=("Path to master species tree."),
                          required=True, type=str)

    # Define other arguments
    parser.add_argument("--outdir",
                        help=("Path to directory to store outputs."),
                        type=str,
                        default="output/")
    parser.add_argument("--min_cov",
                        help=("Minimum coverage of a gene in a given strain "
                              "for the strain to be included in that gene"
                              "tree"),
                        type=float, default=0.8)
    parser.add_argument("--baseml",
                        help=("Path to baseml executable."),
                        default='baseml',
                        type=str)
    parser.add_argument("--resume",
                        help=("Flag that indicates if this is resuming an "
                              "stopped run. If passed, the script will "
                              "assume that genes for which there is already "
                              "an existing output tree file are already "
                              "done and it will skip them. Otherwise, the "
                              "scprit will re-run baseml on every gene."),
                        default=False,
                        action="store_true")

    # Read arguments
    print("Reading arguments")
    args = parser.parse_args()

    # Processing goes here if needed
    # Adding extra options
    args.min_sp = 5

    return args


def subset_aln(infile, outfile, to_keep={},
               informat='fasta', outformat='fasta'):
    """Creates alignment file that contains only sequences in to_keep set."""

    aln = AlignIO.read(infile, informat)

    # Filter aln
    new_aln = Align.MultipleSeqAlignment([])
    for s in aln:
        if s.name in to_keep:
            new_aln.append(s)

    AlignIO.write(alignments=new_aln, handle=outfile, format=outformat)

    return len(new_aln)


def run_baseml(aln_file, tre_file, outdir="output/",
               model=7,
               clock=0,
               tree_format='newick',
               aln_format='fasta',
               baseml_bin='baseml'):
    """Estimate rates from dna alignment and phylogenetic tree
    using PAML's baseml program."""

    # Read original tree and aln
    aln = AlignIO.read(aln_file, aln_format)
    tre = TreeNode.read(tre_file, tree_format)

    # Homogenize tree and aln files
    aln_seqs = {s.id for s in aln}
    tre_tips = {n.name for n in tre.tips()}
    to_keep = aln_seqs & tre_tips
    # Filter aln
    new_aln = Align.MultipleSeqAlignment([])
    for s in aln:
        if s.name in to_keep:
            new_aln.append(s)
    # Filter tre
    new_tre = tre.shear(to_keep)
    new_tre.prune()

    # Write new files (eg. baseml input)
    new_tre_file = os.path.join(outdir, 'tree.newick')
    new_aln_file = os.path.join(outdir, 'aln.phylip')
    TreeNode.write(new_tre, new_tre_file, 'newick')
    # Custom phylip writer
    with open(new_aln_file, 'w') as oh:
        header = ' '.join([str(len(new_aln)),
                           str(len(new_aln[0, :]))]) + "\n"
        oh.write(header)
        for s in new_aln:
            # Need TWO spaces between sequence name and sequence
            line = '  '.join([s.id, str(s.seq)]) + "\n"
            oh.write(line)
    oh.close()

    # run basml
    bml = baseml.Baseml(alignment=new_aln_file, tree=new_tre_file,
                        out_file=os.path.join(outdir, "baseml.out"),
                        working_dir=outdir)
    bml.set_options(model=model, runmode=0, clock=clock)
    res = bml.run(verbose=True, parse=True,
                  command=baseml_bin)

    return(res)


def baseml_aln(aln_file, tre_file, baseml_dir,
               outfile, baseml_bin='baseml'):
    """Single gene aln baseml"""

    # is this redundant with run_baseml?

    # aln_file = os.path.join(aln_dir, g + '.aln.fasta')
    # subset_aln_file = os.path.join(outdir, "gene_alns",
    #                                prefix + '.aln.fasta')
    # gene_baseml_dir = os.path.join(outdir, "baseml", g)
    # gene_tree_file = os.path.join(gene_trees_dir,
    #                               prefix + ".baseml.tre")

    # Error if no alignment exists
    if not os.path.isfile(aln_file):
        raise FileNotFoundError
    # Error if output tree file already exists
    if os.path.isfile(outfile):
        raise FileExistsError

    # Run
    os.mkdir(baseml_dir)
    res = run_baseml(aln_file=aln_file,
                     tre_file=tre_file,
                     outdir=baseml_dir,
                     baseml_bin=baseml_bin)

    tre = TreeNode.read(io.StringIO(res.get('tree')))
    TreeNode.write(tre, file=outfile)

    return outfile


def baseml_all_genes(cov_file, aln_dir, tre_file, outdir="./output/",
                     cov_thres=0.8, n_threshold=5, baseml_bin="baseml",
                     resume=False, cpus=1):
    """Run baseml on all genes with only samples above
    certain coverage threshold."""

    # Prepare output directory structure
    gene_alns_dir = os.path.join(outdir, "gene_alns")
    baseml_dir = os.path.join(outdir, "baseml")
    gene_trees_dir = os.path.join(outdir, "gene_trees")
    if resume:
        # If option resume is true, check first if output dirs exist.
        if not os.path.isdir(outdir):
            os.mkdir(outdir)
        if not os.path.isdir(gene_alns_dir):
            os.mkdir(gene_alns_dir)
        if not os.path.isdir(baseml_dir):
            os.mkdir(baseml_dir)
        if not os.path.isdir(gene_trees_dir):
            os.mkdir(gene_trees_dir)
    else:
        # If option resume is false, it will exit with an error if
        # output files already exist.
        os.mkdir(outdir)
        os.mkdir(gene_alns_dir)
        os.mkdir(baseml_dir)
        os.mkdir(gene_trees_dir)

    # Read coverage information
    covs = pd.read_csv(cov_file, sep="\t", dtype={'gene': np.character})
    covs = covs.set_index('gene')

    # Parallelize with multiprocessing
    with mp.Pool(cpus) as p:
        # Split coverage file into the number of processes.
        for c in np.array_split(covs, cpus):
            # Run baseml on every gene on every chunk
            # _baseml_iterate(covs=c,
            #                 aln_dir=aln_dir,
            #                 outdir=outdir,
            #                 gene_trees_dir=gene_trees_dir,
            #                 tre_file=tre_file,
            #                 cov_thres=cov_thres,
            #                 n_threshold=n_threshold,
            #                 resume=resume,
            #                 baseml_bin=baseml_bin)
            p.apply_async(_baseml_iterate,
                          kwds={'covs': c,
                                'aln_dir': aln_dir,
                                'outdir': outdir,
                                'gene_trees_dir': gene_trees_dir,
                                'tre_file': tre_file,
                                'cov_thres': cov_thres,
                                'n_threshold': n_threshold,
                                'resume': resume,
                                'baseml_bin': baseml_bin})
        # for g, c in covs.iterrows():
        #     # Create file names
        #     aln_file = os.path.join(aln_dir, g + '.aln.fasta')
        #     subset_aln_file = os.path.join(outdir, "gene_alns", g + '.aln.fasta')
        #     gene_baseml_dir = os.path.join(outdir, "baseml", g)
        #     gene_tree_file = os.path.join(gene_trees_dir, g + ".baseml.tre")
        #
        #     # Skip if alignment does not exist
        #     if not os.path.isfile(aln_file):
        #         continue
        #     # Skip if output tree file already exists
        #     if resume and os.path.isfile(gene_tree_file):
        #         continue
        #
        #     # Find samples to keep
        #     to_keep = set(c.index[c >= cov_thres])
        #     # print(g, len(to_keep))
        #     n_samples = subset_aln(infile=aln_file,
        #                            outfile=subset_aln_file,
        #                            to_keep=to_keep)
        #
        #     if n_samples < n_threshold:
        #         continue
        #
        #     # Run baseml
        #     os.mkdir(gene_baseml_dir)
        #     res = run_baseml(aln_file=subset_aln_file, tre_file=tre_file,
        #                      outdir=gene_baseml_dir,
        #                      baseml_bin=baseml_bin)
        #
        #     tre = TreeNode.read(io.StringIO(res.get('tree')))
        #     TreeNode.write(tre, file=gene_tree_file)
        p.close()
        print("Submitted {} processes".format(str(cpus)))
        p.join()


def _baseml_iterate(covs, aln_dir, outdir, gene_trees_dir, tre_file,
                    cov_thres=0.8,
                    n_threshold=5,
                    resume=False,
                    baseml_bin='baseml'):
    """Takes a pandas data frame with gene coverages, and
    runs baseml iteratively in all of the genes.

    No checks are performed. Intended to be called indirectly."""

    # Run baseml on every gene
    for g, c in covs.iterrows():
        # Create file names
        aln_file = os.path.join(aln_dir, g + '.aln.fasta')
        subset_aln_file = os.path.join(outdir, "gene_alns", g + '.aln.fasta')
        gene_baseml_dir = os.path.join(outdir, "baseml", g)
        gene_tree_file = os.path.join(gene_trees_dir, g + ".baseml.tre")

        # Skip if alignment does not exist
        if not os.path.isfile(aln_file):
            continue
        # Skip if output tree file already exists
        if resume and os.path.isfile(gene_tree_file):
            continue

        # Find samples to keep
        to_keep = set(c.index[c >= cov_thres])
        # print(g, len(to_keep))
        n_samples = subset_aln(infile=aln_file,
                               outfile=subset_aln_file,
                               to_keep=to_keep)

        if n_samples < n_threshold:
            continue

        # Run baseml
        if not os.path.isdir(gene_baseml_dir):
            os.mkdir(gene_baseml_dir)
        res = run_baseml(aln_file=subset_aln_file, tre_file=tre_file,
                         outdir=gene_baseml_dir,
                         baseml_bin=baseml_bin)

        tre = TreeNode.read(io.StringIO(res.get('tree')))
        TreeNode.write(tre, file=gene_tree_file)


if __name__ == "__main__":
    args = process_arguments()

    baseml_all_genes(cov_file=args.cov_file,
                     aln_dir=args.aln_dir,
                     tre_file=args.master_tree,
                     outdir=args.outdir,
                     cov_thres=args.min_cov,
                     n_threshold=args.min_sp,
                     baseml_bin=args.baseml,
                     resume=args.resume)

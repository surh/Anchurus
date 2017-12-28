#!/usr/bin/env python
# Copyright (C) 2017 Sur Herrera Paredes
import sutilspy
import os
import fyrd
import argparse

def build_midas_command(sample,read1,read2,bin,args):
    """Build run_midas.py snps command

    Builds a command line command for MIDAS run_midas.py  script,
    in the snps mode

    Args:
        sample (str): Sample ID. typically starts with SRS
        read1 (str): File path to read1 files.
        read2 (str): File path to read2 files
        bin (str): Executable path to run_midas.py
        args (Namespace): Resuts from argparse ArgumentParser.parse_args method

    Returns: A string corresponding to a run_midas.py command
    """

    # Build MIDAS comand
    midas_command = [bin,"snps",args.outdir + "/" + sample,
                         "-1", read1,
                         "-2", read2,
                         "-t","8",
                         "--species_cov",str(args.species_cov),
                         "--mapid", str(args.mapid),
                         "--mapq", str(args.mapq),
                         "--baseq", str(args.baseq),
                         "--readq", str(args.readq),
                         "--aln_cov", "0.75",
                         "-m", 'local']
                         #"--species_id","Haemophilus_parainfluenzae_62356"]
    if args.trim > 0:
        midas_command.extend(["--trim",str(args.trim)])
    if args.discard:
        midas_command.append("--discard")
    if args.baq:
        midas_command.append("--baq")
    if args.adjust_mq:
        midas_command.append("--adjust_mq")
    if args.steps == 'align':
        midas_command.append("--build_db")
        midas_command.append("--align")
    elif args.steps == 'call':
        midas_command.append("--pileup")
    else:
        raise ValueError("Incorrect --steps option ({})".format(args.steps))

    midas_command = " ".join(midas_command)
    # print("#######")
    # print(midas_command)
    # print("#######")

    return(midas_command)


def process_arguments():
    # Read arguments
    parser_format = argparse.ArgumentDefaultsHelpFormatter
    parser = argparse.ArgumentParser(formatter_class=parser_format)
    required = parser.add_argument_group("Required arguments")

    # Define description
    parser.description = ("This script takes a list of sample IDs, "
                          "and submits all those samples to MIDAS "
                          "script <run_midas.py species>. It must "
                          "be run after <run_midas.py species> "
                          "since it requires the output of this script.\n\n"
                          "It creates a job per sample and it can utilize "
                          "fyrd, slurm or torque.")

    # Define required arguments
    required = parser.add_argument_group("Required arguments")
    required.add_argument("--samples",
                          help=("File with samples (one per line) to "
                                "be processed."),
                          type=str, required=True)
    required.add_argument("--indir",
                          help=("Directory where the raw reads files per "
                                "sample are located."),
                          type=str, required=True)
    required.add_argument("--outdir",
                          help=("Directory where to save output. It must "
                                "be also the directory where the output from "
                                "<run_midas.py species> is located."),
                          type=str, required=True)

    # Define other arguments
    # Optional arguments
    parser.add_argument("--sample_col",help = "Column where sample id is located in --samples",
                        default = 1, type = int)
    parser.add_argument("--method", help = "Method to use for submissions. qsub are slurm are kept for legacy reasons, it is strongly recommended to use fyrd always",
                        type = str, default = 'fyrd', choices = ['qsub','slurm','fyrd'])
    parser.add_argument("--logdir", help = "If method is cluster-based, where to store the logfiles",
                         type = str, default = "logs")
    parser.add_argument("--submissions_dir", help = "Directory where to store submission dirs",
                        type = str, default = "submissions")
    parser.add_argument("--queue",
                        help="If method is  'slurm, the partition to use",
                        default="all",
                        choices=['hbfraser','owners','batch','bigmem','all','hns'],
                        type=str)
    parser.add_argument("--memory", help = "Amount of memory to request",
                        default = "10G", type = str)
    parser.add_argument("--time", help = "If method is slurm, amount of time to reserve",
                        type = str, default = "4:00:00")
    parser.add_argument("--species_cov", help = "Include species with greated coverage than identified",
                        type = float, default = 3.0)
    parser.add_argument("--mapid", help = "Minimum mapping identity", type = float,
                        default = 94.0)
    parser.add_argument("--mapq", help = "Minumum mappping quality",
                        type = int, default = 20)
    parser.add_argument("--baseq", help = "Discard bases with quality under the specified value",
                        type = int, default = 30)
    parser.add_argument("--readq", help = "Minimum read quality",
                        type = int, default = 30)
    parser.add_argument("--trim", help = "Trim N base-pairs from 3' end of the read",
                        type = int, default = 0)
    parser.add_argument("--discard", help = "Flag to discard discordant read pairs",
                        action = "store_true")
    parser.add_argument("--baq", help = "Flag to enable per-base alignment quality (BAQ)",
                        action = "store_true")
    parser.add_argument("--adjust_mq", help = "Adjust MAPQ",
                        action = "store_true")
    parser.add_argument("--steps", help = "Steps to perform for <run_midas.py snps. Either build the database and align. Or call SNPs",
                        default = "align", choices = ['align','call'])

    # Read arguments
    print("Reading arguments")
    args = parser.parse_args()

    # Processing goes here if needed
    if args.queue == 'all':
        args.queue = 'hbfraser,owners,hns,normal,bigmem'

    return args


if __name__ == "__main__":
    args = process_arguments()

    # Read samples
    samples = sutilspy.io.return_column(infile = args.samples,
                                        col = args.sample_col,
                                        separator = '\t',
                                        header = False)
    # Prepare directories
    if args.method in ['qsub','slurm','fyrd']:
        if not os.path.isdir(args.logdir):
            print("Creating directory {}".format(args.logdir))
            os.mkdir(args.logdir)
        if not os.path.isdir(args.submissions_dir):
            print("Creating directory {}".format(args.submissions_dir))
            os.mkdir(args.submissions_dir)


    #print(samples)

    pre_commands = []

    # Add module dependencies
    #pre_commands.append("module load MIDAS/1.2.1")
    pre_commands.append("module load MIDAS/1.3.1")
    pre_commands.append("echo MIDAS database is $MIDAS_DB")
    bin = "run_midas.py"

    for sample in samples:
        print("== Processing sample {}".format(sample))
        sample_file_base = args.indir + "/" + sample
        read1 = sample_file_base + "_read1.fastq.bz2"
        read2 = sample_file_base + "_read2.fastq.bz2"

        # We need to check if the read files exist
        if not os.path.isfile(read1):
            raise FileNotFoundError("File {} not found".format(read1))
        if not os.path.isfile(read2):
            raise FileNotFoundError("File {} not found".format(read2))

        # We also need to check if the species coverage file exists
        species_file = args.outdir + "/" + sample + "/species/species_profile.txt"
        if not os.path.isfile(species_file):
            raise FileNotFoundError("File {} not found".format(species_file))

        midas_command = build_midas_command(sample, read1, read2, bin, args)

        # Final list of commands
        commands = pre_commands[:]
        commands.append(midas_command)

        job_name = sample + ".midas"
        logfile = args.logdir + "/midas.snps." + sample + ".log"
        errorfile = args.logdir + "/midas.snps." + sample + ".err"

        submission_file = args.submissions_dir + "/midas.snps." + sample + ".bash"

        # Create submission file or job if method is fyrd
        if args.method == 'qsub':
            with open(submission_file, 'w') as fh:
                # memory = "16000mb"
                nodes = "nodes=1:ppn=8"
                sutilspy.io.write_qsub_submission(fh=fh, commands=commands,
                                                  name=job_name,
                                                  memory=args.memory,
                                                  logfile=logfile,
                                                  errorfile=errorfile,
                                                  nodes=nodes)
            fh.close()
            os.chmod(submission_file, 0o744)
        elif args.method == 'slurm':
            with open(submission_file, 'w') as fh:
                # memory = "16G"
                nodes = "1"
                cpus = "8"
                sutilspy.io.write_slurm_submission(fh=fh,
                                                   commands=commands,
                                                   name=job_name,
                                                   memory=args.memory,
                                                   logfile=logfile,
                                                   errorfile=errorfile,
                                                   queue=args.queue,
                                                   nodes='1',
                                                   cpus='8',
                                                   time=args.time)
            fh.close()
            os.chmod(submission_file, 0o744)
        elif args.method == 'fyrd':
            print("\tCreating fyrd.Job")
            midas_job = fyrd.Job(midas_command,runpath = os.getcwd(),outpath = args.logdir,
                                  scriptpath = args.submissions_dir, clean_files = False,
                                  clean_outputs = False, mem = args.memory, name = job_name,
                                  outfile = "midas.snps." + sample + ".log",
                                  errfile = "midas.snps." + sample + ".err",
                                  partition = args.queue,
                                  nodes = 1, cores = 8, time = args.time,
                                  modules = "MIDAS/1.3.1")

        else:
            raise ValueError("Invalid method {}".format(args.method))



        # Submit submission file
        if args.method == 'qsub':
            #print(submission_file)
            sutilspy.io.qsub_submissions([submission_file],args.logdir)
        elif args.method == 'slurm':
            #print(submission_file)
            sutilspy.io.sbatch_submissions([submission_file], args.logdir)
        elif args.method == 'bash':
            sutilspy.io.run_command(submission_file)
        elif args.method == 'fyrd':
            #midas_job.write(overwrite = True)
            #print("\tWriting submission scripts")
            print("\tSubmitting job")
            midas_job.submit(max_jobs = 1000)
        else:
            raise ValueError("Incorrect method supplie ({})".format(args.method))

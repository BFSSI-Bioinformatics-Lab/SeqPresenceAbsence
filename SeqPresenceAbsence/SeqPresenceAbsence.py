import os
import click
import shutil
import logging
import pandas as pd
from tqdm import tqdm
from pathlib import Path
from copy import deepcopy
from multiprocessing import Pool, cpu_count
from dataclasses import dataclass
from subprocess import Popen, PIPE, DEVNULL

__version__ = "0.3.0"
__author__ = "Forest Dussault"
__email__ = "forest.dussault@canada.ca"

DEPENDENCIES = [
    'blastn',
    'blastx',
    'makeblastdb',
    'muscle',
    'perl'
]

ROOT_DIR = Path(__file__).parent


def print_version(ctx, param, value):
    if not value or ctx.resilient_parsing:
        return
    print(f"Version: {__version__}")
    print(f"Author: {__author__}")
    print(f"Email: {__email__}")
    quit()


def convert_to_path(ctx, param, value):
    if not value or ctx.resilient_parsing:
        return
    value = Path(value)
    return value


@dataclass
class QueryObject:
    """
    Dataclass to store metadata for a sample
    """

    # Must be instantiated with these attributes
    sample_name: str
    fasta_path: Path

    # Updated later in the lifecycle
    fasta_name: str = None
    blastn_path: Path = None
    filtered_blastn_path: Path = None
    blastn_df: pd.DataFrame = None
    blastn_df_processed: pd.DataFrame = None
    target_dict_: dict = None

    def init_target_dict(self, target_list):
        header_list = get_fasta_headers(target_list)
        self.target_dict_ = {header: 0 for header in header_list}

    @staticmethod
    def get_fasta_headers(fasta: Path) -> list:
        """
        Pulls headers any fasta file (e.g. lines starting with >) and returns them as a single list
        """
        fasta_headers = []
        with open(str(fasta)) as f:
            for line in f.readlines():
                if line.startswith(">"):
                    line = line[1:]
                    line = line.strip()
                    fasta_headers.append(line)
        return fasta_headers


@click.command(
    help=f"seqPresenceAbsence is a simple script for querying an input nucleotide FASTA file against a database of "
    f"sequences. Will return an .xlsx and .csv report of presence/absence of the sequences. Version: {__version__}.")
@click.option("-i", "--indir",
              type=click.Path(exists=True),
              required=True,
              help='Path to directory containing FASTA files you want to query',
              callback=convert_to_path)
@click.option("-t", "--targets",
              type=click.Path(exists=True),
              required=True,
              help='Path to multi-FASTA containing targets of interest',
              callback=convert_to_path)
@click.option('-o', '--outdir',
              type=click.Path(exists=False),
              required=True,
              default=None,
              help='Root directory to store all output files',
              callback=convert_to_path)
@click.option('-p', '--perc_identity',
              type=click.FLOAT,
              required=False,
              default=95.00,
              help='Equivalent to the -perc_identity argument in blastn. Defaults to 95.00.')
@click.option('-v', '--verbose',
              is_flag=True,
              default=False,
              help='Set this flag to enable more verbose logging.')
@click.option('--version',
              help='Specify this flag to print the version and exit.',
              is_flag=True,
              is_eager=True,
              callback=print_version,
              expose_value=False)
def cli(indir, targets, outdir, perc_identity, verbose):
    if verbose:
        logging.basicConfig(
            format='\033[92m \033[1m %(asctime)s \033[0m %(message)s ',
            level=logging.DEBUG,
            datefmt='%Y-%m-%d %H:%M:%S')
    else:
        logging.basicConfig(
            format='\033[92m \033[1m %(asctime)s \033[0m %(message)s ',
            level=logging.INFO,
            datefmt='%Y-%m-%d %H:%M:%S')
    logging.info(f"Started seqPresenceAbsence")

    # Check dependencies, fail if necessary
    check_all_dependencies()

    # Set number of concurrent processes for any multiprocessing steps
    n_processes = int(cpu_count() / 8)

    # Validation
    if not targets.suffix == '.fasta':
        logging.error(f"Suffix for --targets argument must be '.fasta', please try again. Your file: {targets}")
        quit()

    # Setup
    os.makedirs(str(outdir), exist_ok=True)

    database = call_makeblastdb(db_file=targets)
    logging.debug(f"Created BLAST database at {database}")

    sample_name_dict = get_sample_name_dict(indir=indir)
    logging.debug(f"Detected {len(sample_name_dict)} samples")

    # Call blast on every sample against the provided database
    query_object_list = []
    for sample_name, infile in tqdm(sample_name_dict.items()):
        # Create QueryObject
        query_object = QueryObject(sample_name=sample_name, fasta_path=infile)

        # Call blastn against sample and parse results
        query_object.blastn_path = call_blast(infile=query_object.fasta_path, database=database, outdir=outdir)
        query_object.blastn_df = parse_blastn(blastn_file=query_object.blastn_path, perc_identity=perc_identity)
        query_object.filtered_blastn_path = export_df(query_object.blastn_df,
                                                      outfile=query_object.blastn_path.with_suffix(
                                                          ".BLASTn_filtered"))
        os.remove(str(query_object.blastn_path))
        query_object.blastn_path = None

        # Initialize the target dictionary
        query_object.init_target_dict(targets)
        query_object.fasta_name = get_fasta_headers(query_object.fasta_path)[0].replace(" ", "_").replace(",", "")
        query_object_list.append(query_object)
    query_object_list = generate_final_report(query_object_list, outdir=outdir)

    # Generate a multifasta per locus containing the sequence for each sample
    loci_dir = outdir / "loci"
    if loci_dir.exists():
        shutil.rmtree(loci_dir)
    loci_dir.mkdir(exist_ok=False)

    # Get master locus list
    master_locus_list = []
    for query_object in query_object_list:
        master_locus_list = list(set(master_locus_list + list(query_object.target_dict_.keys())))

    # Write strand aware sequences to new fasta files
    for locus in tqdm(master_locus_list):
        locus_file = loci_dir / (locus + ".fasta")
        with open(str(locus_file), "w") as f:
            reference_sequence = ""
            for query_object in query_object_list:
                header = f">{query_object.sample_name}\n"
                df = query_object.blastn_df
                try:
                    sequence = str(df[df['sseqid'] == locus]['qseq_strand_aware'].values[0]) + "\n"
                except IndexError:
                    continue
                reference_sequence = str(df[df['sseqid'] == locus]['sseq_strand_aware'].values[0]) + "\n"
                f.write(header)
                f.write(sequence)
            reference_header = f">REFERENCE_TARGET\n"
            if reference_sequence is not "":
                f.write(reference_header)
                f.write(reference_sequence)

    logging.info("Removing empty files from loci dir")
    remove_empty_files_from_dir(in_dir=loci_dir)

    # Create aligned versions of each marker multifasta with muscle
    logging.info("Aligning fasta files in loci dir with MUSCLE")
    aligned_dir = Path(loci_dir / 'aligned')
    aligned_dir.mkdir(exist_ok=False)
    for f in tqdm(list(loci_dir.glob("*.fasta"))):
        call_muscle(infile=f, outfile=(aligned_dir / f.with_suffix(".align.fasta").name))

    logging.info(f"Concatenating all sequences in {aligned_dir}")
    sample_ids = [query_object.sample_name for query_object in query_object_list]
    concatenated_sequences = concatenate_sequence_directory(sample_ids=sample_ids, sequence_directory=aligned_dir,
                                                            n_processes=n_processes,
                                                            outdir=outdir)
    logging.info(f"Concatenated sequences available at {concatenated_sequences}")
    logging.info("Script Complete!")


def remove_empty_files_from_dir(in_dir: Path):
    file_list = list(in_dir.glob("*"))
    file_list = [f for f in file_list if f.is_file()]
    for f in file_list:
        if f.lstat().st_size == 0:
            f.unlink()


def get_sample_name_dict(indir: Path) -> dict:
    fasta_files = list(indir.glob("*.fna"))
    fasta_files += list(indir.glob("*.fasta"))
    fasta_files += list(indir.glob("*.fa"))
    fasta_files += list(indir.glob("*.ffn"))
    fasta_files += list(indir.glob("*.fas"))

    sample_name_dict = {}
    for f in fasta_files:
        sample_name = f.with_suffix("").name
        sample_name_dict[sample_name] = f
    return sample_name_dict


def export_df(df: pd.DataFrame, outfile: Path) -> Path:
    df.to_csv(str(outfile), sep="\t", index=None)
    return outfile


def get_fasta_headers(fasta: Path) -> list:
    """
    Pulls headers any fasta file (e.g. lines starting with >) and returns them as a single list
    """
    fasta_headers = []
    with open(str(fasta)) as f:
        for line in f.readlines():
            if line.startswith(">"):
                line = line[1:]
                line = line.strip()
                header = line.split(" ")[0]

                # TODO: Bug - contigs that have "|" in them will break this fix, whereas targets with "|" need this fix?
                # This is really sketchy
                # if '|' in header:
                #     header = header.split("|")[1]
                fasta_headers.append(header)
    return fasta_headers


def parse_blastn(blastn_file: Path, perc_identity: float, header: list = None) -> pd.DataFrame:
    """
    Parses *.BLASTn file generated by call_blastn(), then returns the df
    :param blastn_file: file path to *.BLASTn file
    :param perc_identity: Equivalent of -perc_identity in blastn
    :param header: List of values of expected headers in blastn file
    :return: DataFrame contents of *.BLASTn file
    """
    if header is None:
        header = ["qseqid", "stitle", "sseqid", "slen", "length", "qstart", "qend",
                  "pident", "score", "sstrand", "bitscore", "qseq", "sseq"]
    df = pd.read_csv(blastn_file, delimiter="\t", names=header, index_col=None)
    df['lratio'] = df['length'] / df['slen']
    df['qseq_strand_aware'] = df.apply(get_reverse_complement_qseq, axis=1)
    df['sseq_strand_aware'] = df.apply(get_reverse_complement_sseq, axis=1)

    df = df.query(f"lratio >= 0.90 & pident >= {perc_identity} & lratio <=1.10")
    df = df.sort_values(by=['bitscore'], ascending=False)
    return df


def get_highest_bitscore_hit(df: pd.DataFrame) -> pd.DataFrame:
    df = df.sort_values(by=['bitscore']).reset_index()
    return df.head(1)


def call_blast(infile: Path, database: Path, outdir: Path):
    out_blast = outdir / infile.with_suffix(".BLASTn").name
    blast_params = "6 qseqid stitle sseqid slen length qstart qend pident score sstrand bitscore qseq sseq"
    blast_cmd = f"blastn -db {database} -query {infile} -word_size 15 -outfmt '{blast_params}' > {out_blast}"
    run_subprocess(blast_cmd)
    return out_blast


def combine_dataframes(dfs: [pd.DataFrame]) -> pd.DataFrame:
    """
    Receives a list of DataFrames and concatenates them. They must all have the same header.
    :param dfs: List of DataFrames
    :return: Single concatenated DataFrame
    """
    df = pd.concat(dfs, sort=False)
    return df


def export_multifasta(blastn_df: pd.DataFrame, sample_name: str, outdir: Path):
    outname = outdir / (sample_name + ".fasta")
    with open(str(outname), "w") as f:
        for index, row in blastn_df.iterrows():
            header = ">" + sample_name + "_" + row['sseqid']
            sequence = row['qseq_strand_aware']
            f.write(header + "\n")
            f.write(sequence + "\n")


def filter_empty_hits(df: pd.DataFrame):
    df = df.reset_index()
    rows = df.iterrows()
    next(rows)
    next(rows)
    rows_to_remove = []
    for i, row in rows:
        presabs_list = list(row.values[1:])
        if 1 not in presabs_list:
            rows_to_remove.append(i)
    df = df.drop(df.index[tuple([rows_to_remove])])
    return df


def generate_final_report(query_object_list: [QueryObject], outdir: Path) -> [QueryObject]:
    df_master_list = []
    for query_object in query_object_list:
        query_object.blastn_df = pd.read_csv(query_object.filtered_blastn_path, sep='\t')

        export_multifasta(blastn_df=query_object.blastn_df, sample_name=query_object.sample_name, outdir=outdir)

        sseqid_series = list(query_object.blastn_df['sseqid'])
        for target in query_object.target_dict_.keys():
            if target in sseqid_series:
                query_object.target_dict_[target] += 1

        query_object.blastn_df_processed = pd.DataFrame(list(query_object.target_dict_.items()),
                                                        columns=['target', 'count'])
        query_object.blastn_df_processed['sample_name'] = query_object.sample_name
        query_object.blastn_df_processed['fasta_name'] = query_object.fasta_name
        df_master_list.append(query_object.blastn_df_processed)

    df_combined = combine_dataframes(df_master_list)
    df_pivot = df_combined.pivot_table(index=['sample_name', 'fasta_name'],
                                       columns='target',
                                       values='count').reset_index()

    csv_path = outdir / f'TargetOutput.tsv'
    csv_path_t = outdir / f'TargetOutput_transposed.tsv'
    xlsx_path = outdir / f'TargetOutput.xlsx'
    xlsx_path_t = outdir / f'TargetOutput_transposed.xlsx'
    xlsx_path_t_filtered = outdir / f'TargetOutput_transposed_filtered.xlsx'

    # Excel report
    writer = pd.ExcelWriter(str(xlsx_path), engine='xlsxwriter')
    df_pivot.to_excel(writer, sheet_name='TargetOutput', index=None)
    worksheet = writer.sheets['TargetOutput']
    worksheet.conditional_format('B2:AMJ4000', {'type': '2_color_scale'})
    writer.save()

    # Transposed Excel report
    # df_pivot.set_index('sample_name', inplace=True)
    df_transposed = df_pivot.transpose()
    writer = pd.ExcelWriter(str(xlsx_path_t), engine='xlsxwriter')
    df_transposed.to_excel(writer, sheet_name='TargetOutput', header=None)
    worksheet = writer.sheets['TargetOutput']
    worksheet.conditional_format('A2:AMJ4000', {'type': '2_color_scale'})
    writer.save()

    # Transposed filtered Excel report
    # Filter out rows without any hits
    df_transposed_filtered = filter_empty_hits(df=df_transposed)
    writer = pd.ExcelWriter(str(xlsx_path_t_filtered), engine='xlsxwriter')
    df_transposed_filtered.to_excel(writer, sheet_name='TargetOutput', header=None, index=None)
    worksheet = writer.sheets['TargetOutput']
    worksheet.conditional_format('A2:AMJ4000', {'type': '2_color_scale'})
    writer.save()

    # CSV files
    df_pivot.to_csv(csv_path, sep='\t', index=None)
    df_transposed.to_csv(csv_path_t, sep='\t', header=None)

    logging.info(f"Created .csv of count data at {csv_path}")
    logging.info(f"Created .xlsx of count data at {xlsx_path}")

    return query_object_list


def dependency_check(dependency: str) -> bool:
    """
    Checks if a given program is present in the user's $PATH
    :param dependency: String of program name
    :return: True if program is in $PATH, False if not
    """
    check = shutil.which(dependency)
    if check is not None:
        return True
    else:
        return False


def check_all_dependencies():
    # Dependency check
    logging.info("Conducting dependency check...")
    dependency_dict = dict()
    for dependency in DEPENDENCIES:
        dependency_dict[dependency] = dependency_check(dependency)
    if False in dependency_dict.values():
        logging.error("ERROR: Cannot locate some dependencies in $PATH...")
        for key, value in dependency_dict.items():
            if not value:
                logging.error(f"Dependency missing: {key}")
        quit()
    else:
        for key, value in dependency_dict.items():
            logging.debug(f"Dependency {key}: {value}")
    logging.info("Dependencies OK")


def call_makeblastdb(db_file: Path) -> Path:
    """
    Makes a system call to makeblastdb on a given database file. Can handle *.gz, *.fasta, or no suffix.
    """
    db_name = db_file.with_suffix(".blastDB")
    db_type = 'nucl'
    if db_file.suffix == ".gz":
        cmd = f"gunzip -c {db_file} | makeblastdb -in - -parse_seqids -out {db_name} -title {db_name} "
        cmd += f"-dbtype {db_type}"
        run_subprocess(cmd, get_stdout=True)
    elif db_file.suffix == ".fasta":
        cmd = f"makeblastdb -in {db_file} -parse_seqids -out {db_name} -title {db_name} "
        cmd += f"-dbtype {db_type}"
        run_subprocess(cmd, get_stdout=True)
    elif db_file.suffix == "":
        os.rename(str(db_file), str(db_file.with_suffix(".fasta")))
        cmd = f"makeblastdb -in {db_file} -parse_seqids -out {db_name} -title {db_name} "
        cmd += f"-dbtype {db_type}"
        run_subprocess(cmd, get_stdout=True)
    else:
        logging.debug("Invalid file format provided to call_makeblastdb()")
    return db_name


def call_muscle(infile: Path, outfile: Path = None):
    if outfile is None:
        outfile = infile.with_suffix(".align.fasta")
    cmd = f"muscle -in {infile} -out {outfile} -maxiters 1"
    run_subprocess(cmd, get_stdout=True)


def run_subprocess(cmd: str, get_stdout: bool = False) -> str:
    """
    Will run an external bash command via shell. Set get_stdout to True to retrieve contents of stdout.
    """
    if get_stdout:
        p = Popen(cmd, shell=True, stdout=PIPE, stderr=PIPE)
        out, err = p.communicate()
        out = out.decode().strip()
        err = err.decode().strip()
        if out != "":
            return out
        elif err != "":
            return err
        else:
            return ""
    else:
        p = Popen(cmd, shell=True)
        p.wait()


def get_reverse_complement_qseq(row) -> str:
    """
    Takes DataFrame row and returns the reverse complement of the sequence if the strand is 'minus'
    :param row:
    :return:
    """
    complement_dict = {'A': 'T',
                       'C': 'G',
                       'G': 'C',
                       'T': 'A'}
    sequence = row['qseq'].upper()
    if row['sstrand'] == 'minus':
        reverse_complement = "".join(complement_dict.get(base, base) for base in reversed(sequence))
        return reverse_complement
    else:
        return sequence


def get_reverse_complement_sseq(row) -> str:
    """
    TODO: DRY. This is basically the exact same fn as get_reverse_complement_qseq.
    Takes DataFrame row and returns the reverse complement of the sequence if the strand is 'minus'
    :param row:
    :return:
    """
    complement_dict = {'A': 'T',
                       'C': 'G',
                       'G': 'C',
                       'T': 'A'}
    sequence = row['sseq'].upper()
    if row['sstrand'] == 'minus':
        reverse_complement = "".join(complement_dict.get(base, base) for base in reversed(sequence))
        return reverse_complement
    else:
        return sequence


def concatenate_sequence_directory(sample_ids: [str], sequence_directory: Path, n_processes: int, outdir: Path) -> Path:
    """
    Given a sequence directory containing aligned multi-FASTA files, will attempt to concatenate all of the sequences
    into a single file. Minimal example showing input and output (concatenated_sequences.fasta):

    files in directory:
        gene_1.fasta
        gene_2.fasta

    contents of gene1.fasta:
        >sample_1
        ATCG
        >sample2
        ATTT

    contents of gene2.fasta:
        >sample_1
        GGGGAGGGGGTCA

    concatenated_sequences.fasta:
        >sample_1
        ATCGGGGGAGGGGGTCA
        >sample_2
        ATTTNNNNNNNNNNNNN

    - Provided sample_ids list must contain exact matches to headers in FASTA files
    - All sequences in a single FASTA file must be the exact same length; this is done through an aligner like MUSCLE
    - Expectation is that every FASTA file will have the same headers, or a subset of headers present in sample_ids
    """
    concat_seqs_dict = generate_concat_seqs_dict(sample_ids=sample_ids, indir=sequence_directory,
                                                 n_processes=n_processes)
    outfile = write_concat_seqs_dict(concat_seqs_dict=concat_seqs_dict, outdir=outdir)
    return outfile


def write_concat_seqs_dict(concat_seqs_dict: dict, outdir: Path) -> Path:
    outfile = outdir / 'concatenated_sequences.fasta'
    if outfile.exists():
        outfile.unlink()
    outfile_ = open(str(outfile), 'a+')
    for sample_id, sequence in concat_seqs_dict.items():
        outfile_.write(f">{sample_id}\n")
        outfile_.write(f"{sequence}\n")
    outfile_.close()
    return outfile


def generate_concat_seqs_dict(sample_ids: set, indir: Path, n_processes=4) -> dict:
    # Potentially takes a lot of RAM. Stores the concatenated sequences for each sample.
    sequence_storage = {sample_id: "" for sample_id in sample_ids}
    sequence_storage['REFERENCE_TARGET'] = ''
    fasta_files = sorted(list(indir.glob("*.fasta")))

    # Set # of concurrent processes to run
    pool = Pool(processes=n_processes)
    cluster_dicts = [
        pool.apply_async(populate_template_dict, args=(sequence_storage, f, sample_ids)) for f in fasta_files
    ]
    cluster_dicts = [result.get() for result in cluster_dicts]

    # Merge all of the dictionaries into one
    for d in cluster_dicts:
        for cluster, sequence in d.items():
            sequence_storage[cluster] += sequence
    return sequence_storage


def populate_template_dict(template_dict: dict, cluster_file: Path, sample_ids: list):
    with open(str(cluster_file), 'r') as f:
        cluster_dict = deepcopy(template_dict)
        lines = f.readlines()
        cluster_samples = []
        seq_length = 0
        seq_lengths = []
        for line in lines:
            line = line.strip()
            if line.startswith(">"):
                # Add the sequence length to the ongoing tracking list
                if seq_length != 0:
                    seq_lengths.append(seq_length)
                seq_length = 0
                # Remove '>'
                header = line.replace(">", "")
                # Add this header to cluster_samples to keep track of which samples we have sequence for
                cluster_samples.append(header)
            else:
                cluster_dict[header] += line
                seq_length += len(line)
        try:
            assert len(set(seq_lengths)) == 1
        except AssertionError:
            print(f"ERROR: Varying sequence lengths detected in {f.name}")
            print(seq_lengths)
            quit()
        cluster_samples = set(cluster_samples)
        for s in sample_ids:
            if s not in cluster_samples:
                cluster_dict[s] += ('N' * seq_lengths[0])
        return cluster_dict


if __name__ == "__main__":
    cli()

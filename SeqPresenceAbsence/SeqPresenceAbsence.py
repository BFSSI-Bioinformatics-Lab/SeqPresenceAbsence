import os
import click
import shutil
import logging
import pandas as pd
from tqdm import tqdm
from pathlib import Path
from dataclasses import dataclass
from subprocess import Popen, PIPE, DEVNULL

__version__ = "0.1.3"
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
FASCONCAT = ROOT_DIR / 'FASconCAT-G_v1.04.pl'


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
@click.option('-k', '--keep_db_seqs',
              is_flag=True,
              default=False,
              help='Set this flag to keep the target sequence in addition to the query sequence from BLAST.')
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
def cli(indir, targets, outdir, perc_identity, keep_db_seqs, verbose):
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

    # TODO: Dynamically detect and adjust to sequence type for targets/query (nucl vs prot)
    # if protein:
    #     logging.debug("Using blastx for alignments")
    # else:
    #     logging.debug("Using blastn for alignments")

    check_all_dependencies()

    if not targets.suffix == '.fasta':
        logging.error(f"Suffix for --targets argument must be '.fasta', please try again. Your file: {targets}")
        quit()

    os.makedirs(str(outdir), exist_ok=True)

    database = call_makeblastdb(db_file=targets)
    logging.debug(f"Created BLAST database at {database}")

    sample_name_dict = get_sample_name_dict(indir=indir)
    logging.debug(f"Detected {len(sample_name_dict)} samples")

    query_object_list = []
    for sample_name, infile in tqdm(sample_name_dict.items()):
        # Create QueryObject
        query_object = QueryObject(sample_name=sample_name, fasta_path=infile)

        # Call blastn against sample and parse results
        query_object.blastn_path = call_blast(infile=query_object.fasta_path, database=database, outdir=outdir,
                                              keep_db_seqs=keep_db_seqs)
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
    os.makedirs(str(loci_dir), exist_ok=True)

    # Get master locus list
    master_locus_list = []
    for query_object in query_object_list:
        master_locus_list = list(set(master_locus_list + list(query_object.target_dict_.keys())))

    for locus in tqdm(master_locus_list):
        locus_file = loci_dir / (locus + ".fas")
        with open(str(locus_file), "w") as f:
            for query_object in query_object_list:
                header = f">{query_object.sample_name}\n"
                df = query_object.blastn_df
                try:
                    sequence = str(df[df['sseqid'] == locus]['qseq_strand_aware'].values[0]) + "\n"
                except IndexError:
                    continue
                f.write(header)
                f.write(sequence)

    logging.info("Removing empty files from loci dir")
    remove_empty_files_from_dir(in_dir=loci_dir)

    # Create alignmed versions of each marker multifasta with muscle
    logging.info("Aligning fasta files in loci dir with MUSCLE")
    aligned_dir = Path(loci_dir / 'aligned')
    aligned_dir.mkdir(exist_ok=True)
    for f in tqdm(list(loci_dir.glob("*.fas"))):
        call_muscle(infile=f, outfile=(aligned_dir / f.with_suffix(".align.fas").name))

    logging.info(f"Calling FASconCAT on contents of {aligned_dir}")
    call_fasconcat(target_dir=aligned_dir, fasconcat_exec=FASCONCAT)

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

    sample_name_dict = {}
    for f in fasta_files:
        sample_name = f.with_suffix("").name
        sample_name_dict[sample_name] = f
    return sample_name_dict


def export_df(df: pd.DataFrame, outfile: Path) -> Path:
    df.to_csv(str(outfile), sep="\t", index=None)
    return outfile


def create_target_dict(targets: Path) -> dict:
    header_list = get_fasta_headers(targets)
    return {header: 0 for header in header_list}


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
                # This is really sketchy
                if '|' in header:
                    header = header.split("|")[1]
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
    df['qseq_strand_aware'] = df.apply(get_reverse_complement_row, axis=1)
    df = df.query(f"lratio >= 0.90 & pident >= {perc_identity} & lratio <=1.10")
    df = df.sort_values(by=['bitscore'], ascending=False)
    return df


def get_highest_bitscore_hit(df: pd.DataFrame) -> pd.DataFrame:
    df = df.sort_values(by=['bitscore']).reset_index()
    return df.head(1)


def call_blast(infile: Path, database: Path, outdir: Path, keep_db_seqs: bool):
    out_blast = outdir / infile.with_suffix(".BLASTn").name
    blast_params = "6 qseqid stitle sseqid slen length qstart qend pident score sstrand bitscore qseq "
    if keep_db_seqs:
        blast_params += 'sseq'

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

    # TODO: Split df_pivot into those with NO hits and those with some

    csv_path = outdir / f'TargetOutput.tsv'
    csv_path_t = outdir / f'TargetOutput_transposed.tsv'
    xlsx_path = outdir / f'TargetOutput.xlsx'
    xlsx_path_t = outdir / f'TargetOutput_transposed.xlsx'

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


def call_raxml_ng(fasta: Path):
    """
    raxml-ng --msa /mnt/QuizBoy/probeDetectionTesting/TreeTest/out/loci/aligned/FcC_supermatrix.fas  --model GTR+G --threads 2
    """
    cmd = f"raxml-ng --msa {fasta} --model GTR+G --threads 2"


def call_raxml(fasta: Path, n_cpu: int):
    """
    raxmlHPC-PTHREADS-SSE3 -s FcC_smatrix.phy.reduced -n VDB_Tree_Galapagos -p 12345 -x 12345 -N 1000 -m GTRGAMMA -T 50 -f a
    """
    tree_name = fasta.with_suffix(".tree").name
    cmd = f"raxmlHPC-PTHREADS-SSE3 -s {fasta} -n {tree_name} -p 12345 -N 1000 -m GTRGAMMA -T {n_cpu} -f a"
    run_subprocess(cmd)


def call_muscle(infile: Path, outfile: Path = None):
    if outfile is None:
        outfile = infile.with_suffix(".align.fasta")
    cmd = f"muscle -in {infile} -out {outfile} -maxiters 1"
    run_subprocess(cmd, get_stdout=True)


def call_fasconcat(target_dir: Path, fasconcat_exec: Path):
    cmd = f"perl {fasconcat_exec} -s -p"
    p = Popen(cmd, shell=True, cwd=str(target_dir), stdout=DEVNULL, stderr=DEVNULL)
    p.wait()


def run_subprocess(cmd: str, get_stdout: bool = False) -> str:
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


def get_reverse_complement_row(row) -> str:
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


if __name__ == "__main__":
    cli()

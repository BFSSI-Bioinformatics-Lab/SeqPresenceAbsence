import os
import click
import shutil
import logging
import pandas as pd
from pathlib import Path
from dataclasses import dataclass
from subprocess import Popen, PIPE

DEPENDENCIES = [
    'blastn',
    'makeblastdb'
]


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
    help="SequencePresenceAbsence is a simple script for querying an input FASTA file against a database of probes. "
         "Will return a visualization of presence/absence of the probes.")
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
@click.option('-v', '--verbose',
              is_flag=True,
              default=False,  # Set this to false eventually
              help='Set this flag to enable more verbose logging.')
def cli(indir, targets, outdir, verbose):
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

    logging.info(f"Started {Path(__file__).name}")

    check_all_dependencies()

    os.makedirs(str(outdir), exist_ok=True)

    database = call_makeblastdb(targets)
    logging.debug(f"Created BLAST database at {database}")

    sample_name_dict = get_sample_name_dict(indir=indir)
    logging.debug(f"Detected {len(sample_name_dict)} samples")

    query_object_list = []
    with click.progressbar(sample_name_dict.items(), length=len(sample_name_dict), label="BLASTn") as bar:
        for sample_name, infile in bar:
            query_object = QueryObject(sample_name=sample_name, fasta_path=infile)
            query_object.blastn_path = call_blastn(query_object.fasta_path, database, outdir)
            query_object.blastn_df = parse_blastn(query_object.blastn_path)
            query_object.filtered_blastn_path = export_df(query_object.blastn_df,
                                                          outfile=query_object.blastn_path.with_suffix(
                                                              ".BLASTn_filtered"))
            os.remove(str(query_object.blastn_path))
            query_object.blastn_path = None
            query_object.init_target_dict(targets)
            query_object.fasta_name = get_fasta_headers(query_object.fasta_path)[0].replace(" ", "_").replace(",", "")
            query_object_list.append(query_object)

    generate_final_report(query_object_list, outdir=outdir)
    logging.info("Script Complete!")


def get_sample_name_dict(indir: Path) -> dict:
    fasta_files = list(indir.glob("*.fna"))
    fasta_files += list(indir.glob("*.fasta"))
    fasta_files += list(indir.glob("*.fa"))

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
                fasta_headers.append(header)
    return fasta_headers


def parse_blastn(blastn_file: Path, header: list = None) -> pd.DataFrame:
    """
    Parses *.BLASTn file generated by call_blastn(), then returns the df
    :param blastn_file: file path to *.BLASTn file
    :param header: List of values of expected headers in blastn file
    :return: DataFrame contents of *.BLASTn file
    """
    if header is None:
        header = ["qseqid", "stitle", "sseqid", "slen", "length", "qstart", "qend",
                  "pident", "score", "sstrand", "qseq", "bitscore"]
    df = pd.read_csv(blastn_file, delimiter="\t", names=header, index_col=None)
    df['lratio'] = df['length'] / df['slen']
    df['qseq_strand_aware'] = df.apply(get_reverse_complement_row, axis=1)
    # TODO: Ask Nick about values to hardcode here
    df = df.query("lratio >= 0.95 & pident >= 95.0 & lratio <=1.05")
    return df


def get_highest_bitscore_hit(df: pd.DataFrame) -> pd.DataFrame:
    df = df.sort_values(by=['bitscore']).reset_index()
    return df.head(1)


def call_blastn(infile: Path, database: Path, outdir: Path):
    # Need to perform first BLAST phase, if something is found go onto the next
    # logging.info(f"Running blastn on {infile}")
    out_blast = outdir / infile.with_suffix(".BLASTn").name
    blast_cmd = f"blastn -db {database} -query {infile} -word_size 15 -outfmt " \
                f"'6 qseqid stitle sseqid slen length qstart qend pident score sstrand qseq bitscore' " \
                f"> {out_blast}"
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


def generate_final_report(query_object_list: [QueryObject], outdir: Path):
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
    xlsx_path = outdir / f'TargetOutput.xlsx'
    writer = pd.ExcelWriter(str(xlsx_path), engine='xlsxwriter')
    df_pivot.to_excel(writer, sheet_name='TargetOutput', index=None)
    worksheet = writer.sheets['TargetOutput']
    worksheet.conditional_format('B2:AMJ4000', {'type': '2_color_scale'})
    writer.save()

    df_pivot.to_csv(csv_path, sep='\t', index=None)

    logging.info(f"Created .csv of count data at {csv_path}")
    logging.info(f"Created .xlsx of count data at {xlsx_path}")


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
    :param db_file: Path to database file
    """
    db_name = db_file.with_suffix(".blastDB")
    if db_file.suffix == ".gz":
        cmd = f"gunzip -c {db_file} | makeblastdb -in - -parse_seqids -dbtype nucl -out {db_name} -title {db_name}"
        run_subprocess(cmd, get_stdout=True)
    elif db_file.suffix == ".fasta":
        cmd = f"makeblastdb -in {db_file} -parse_seqids -dbtype nucl -out {db_name} -title {db_name}"
        run_subprocess(cmd, get_stdout=True)
    elif db_file.suffix == "":
        os.rename(str(db_file), str(db_file.with_suffix(".fasta")))
        cmd = f"makeblastdb -in {db_file} -parse_seqids -dbtype nucl -out {db_name} -title {db_name}"
        run_subprocess(cmd, get_stdout=True)
    else:
        logging.debug("Invalid file format provided to call_makeblastdb()")
    return db_name


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
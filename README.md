# seqPresenceAbsence

### Requirements
- Python 3.6
- [ncbi-blast+](https://blast.ncbi.nlm.nih.gov/Blast.cgi?PAGE_TYPE=BlastDocs&DOC_TYPE=Download) (makeblastdb and blastn must be in your $PATH)

### Installation
```
pip install seqPresenceAbsence
```

### Usage
```
Usage: seqPresenceAbsence [OPTIONS]

  seqPresenceAbsence is a simple script for querying an input nucleotide
  FASTA file against a database of sequences. Will return an .xlsx and .csv
  report of presence/absence of the sequences. Version: 0.1.3.

Options:
  -i, --indir PATH           Path to directory containing FASTA files you want
                             to query  [required]
  -t, --targets PATH         Path to multi-FASTA containing targets of
                             interest  [required]
  -o, --outdir PATH          Root directory to store all output files
                             [required]
  -p, --perc_identity FLOAT  Equivalent to the -perc_identity argument in
                             blastn. Defaults to 95.00.
  -k, --keep_db_seqs         Set this flag to keep the target sequence in
                             addition to the query sequence from BLAST.
  -v, --verbose              Set this flag to enable more verbose logging.
  --version                  Specify this flag to print the version and exit.
  --help                     Show this message and exit.
```
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

  seqPresenceAbsence is a simple script for querying an input FASTA file
  against a database of sequences. Will return an .xlsx and .csv report
  of presence/absence of the sequences.

Options:
  -i, --indir PATH    Path to directory containing FASTA files you want to
                      query  [required]
  -t, --targets PATH  Path to multi-FASTA containing targets of interest
                      [required]
  -o, --outdir PATH   Root directory to store all output files  [required]
  -v, --verbose       Set this flag to enable more verbose logging.
  --help              Show this message and exit.
```
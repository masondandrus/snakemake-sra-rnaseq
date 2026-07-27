"""Build a transcript-to-gene map from a GENCODE-style transcript FASTA.

GENCODE headers are pipe-delimited:
    >ENST00000456328.2|ENSG00000290825.1|-|OTTHUMT...|DDX11L2-202|DDX11L2|1657|lncRNA|
     [0]                [1]              ...                      [5]

Fields 0, 1 and 5 give transcript ID, gene ID and gene symbol.
Version suffixes are stripped so IDs join cleanly against Ensembl/biomaRt.
"""

import gzip
import sys

fasta = snakemake.input[0]
out = snakemake.output[0]
log = open(snakemake.log[0], "w")

opener = gzip.open if str(fasta).endswith(".gz") else open

n = 0
with opener(fasta, "rt") as fh, open(out, "w") as w:
    w.write("tx_id\tgene_id\tgene_symbol\n")
    for line in fh:
        if not line.startswith(">"):
            continue
        fields = line[1:].strip().split("|")
        if len(fields) < 6:
            print(f"skipping unparseable header: {line[:80]}", file=log)
            continue
        tx, gene, symbol = fields[0], fields[1], fields[5]
        w.write(f"{tx.split('.')[0]}\t{gene.split('.')[0]}\t{symbol}\n")
        n += 1

print(f"wrote {n} transcript-gene pairs to {out}", file=log)
log.close()

if n == 0:
    sys.exit("No transcripts parsed — is this a GENCODE-format FASTA?")

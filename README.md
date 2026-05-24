# Ensembl Gene Isoform & Exon Fetcher

A minimal Python script to fetch exon and isoform-level information for Ensembl human genes via the REST API.

## Features

- Retrieves exon coordinates, nucleotide sequences, and CDS-relative positions
- Reports amino acids encoded by complete codons per exon
- Handles split-codon junctions at exon boundaries (`junctional_Cterm`)
- Batch downloads with retry logic and rate-limiting
- Optional protein-sequence validation against Ensembl

## Installation

```bash
pip install requests pandas
```

## Usage

### Single gene
```bash
python scripts/fetch_isoforms.py \
  --genes ENSG00000141510 \
  --output-dir data/human_gene_exons
```

### Multiple genes
```bash
python scripts/fetch_isoforms.py \
  --genes ENSG00000141510,ENSG00000139618 \
  --output-dir data/human_gene_exons
```

### From file
```bash
python scripts/fetch_isoforms.py \
  --input-file genes.txt \
  --output-dir data/human_gene_exons
```

## Output

One TSV per gene: `data/human_gene_exons/<ENSG>_exons.tsv`

| Column | Description |
|--------|-------------|
| `gene_id` | Ensembl gene ID |
| `transcript_id` | Ensembl transcript ID |
| `exon_id` | Ensembl exon ID |
| `genomic_start` / `genomic_end` | Chromosomal coordinates |
| `exon_sequence` | Nucleotide sequence |
| `exon_protein_sequence` | Amino acids from complete codons inside the exon |
| `junctional_Cterm` | C-terminal amino acid from a split codon at the exon boundary |

## Options

| Flag | Description |
|------|-------------|
| `--overwrite` | Replace existing output files |
| `--reverse-complement-negative-strand` | Reverse-complement negative-strand exon sequences |
| `--no-protein-validation` | Skip Ensembl protein sequence validation |
| `-v, --verbose` | Print preview of each output dataframe |

## Requirements

- Python 3.8+
- `requests`
- `pandas`

#!/usr/bin/env python3
"""Fetch two test genes using only API-retrieved data."""

import sys
from pathlib import Path
from typing import Any, Dict, List

import pandas as pd

from exon_protein_resolver.api import EnsemblClient


OUTPUT_COLUMNS = [
    "gene_id",
    "gene_name",
    "strand",
    "transcript_id",
    "transcript_name",
    "transcript_biotype",
    "is_coding",
    "protein_id",
    "exon_id",
    "total_exons",
    "exon_number",
    "chromosome",
    "genomic_start",
    "genomic_end",
    "cds_genomic_start",
    "cds_genomic_end",
    "exon_length",
    "exon_sequence",
]


def fetch_gene_dataframe(gene_id: str) -> pd.DataFrame:
    client = EnsemblClient()

    print(f"\n=== Fetching {gene_id} ===", file=sys.stderr)
    gene_data = client.get_gene_info(gene_id)

    transcripts = gene_data.get("Transcript", [])
    all_exon_ids: List[str] = []

    for transcript in transcripts:
        for exon in transcript.get("Exon", []):
            eid = exon.get("id")
            if eid:
                all_exon_ids.append(eid)

    exon_sequences = client.get_sequences_batch(sorted(set(all_exon_ids)), "genomic")

    rows: List[Dict[str, Any]] = []

    for transcript in transcripts:
        transcript_id = transcript["id"]
        transcript_name = transcript.get("display_name", "")
        transcript_biotype = transcript.get("biotype", "")
        transcript_strand = int(transcript.get("strand", 1))
        translation = transcript.get("Translation")
        protein_id = translation.get("id") if translation else None
        cds_start = translation.get("start") if translation else None
        cds_end = translation.get("end") if translation else None

        exons = transcript.get("Exon", [])
        if transcript_strand == 1:
            exons = sorted(exons, key=lambda e: int(e["start"]))
        else:
            exons = sorted(exons, key=lambda e: -int(e["end"]))

        for exon_idx, exon in enumerate(exons, start=1):
            exon_id = exon["id"]
            start = int(exon.get("start", 0))
            end = int(exon.get("end", 0))
            strand = int(exon.get("strand", transcript_strand))
            chrom = exon.get("seq_region_name", "")

            rows.append({
                "gene_id": gene_id,
                "gene_name": gene_data.get("display_name", ""),
                "strand": strand,
                "transcript_id": transcript_id,
                "transcript_name": transcript_name,
                "transcript_biotype": transcript_biotype,
                "is_coding": bool(translation),
                "protein_id": protein_id or "",
                "exon_id": exon_id,
                "total_exons": len(exons),
                "exon_number": exon_idx,
                "chromosome": chrom,
                "genomic_start": start,
                "genomic_end": end,
                "cds_genomic_start": cds_start if cds_start is not None else "",
                "cds_genomic_end": cds_end if cds_end is not None else "",
                "exon_length": end - start + 1 if start and end else 0,
                "exon_sequence": exon_sequences.get(exon_id, ""),
            })

    return pd.DataFrame(rows)[OUTPUT_COLUMNS]


if __name__ == "__main__":
    out_dir = Path("test_output")
    out_dir.mkdir(parents=True, exist_ok=True)

    # Positive strand: TNF (ENSG00000232810)
    df_pos = fetch_gene_dataframe("ENSG00000232810")
    pos_path = out_dir / "ENSG00000232810_exons.tsv"
    df_pos.to_csv(pos_path, sep="\t", index=False)
    print(f"\nSaved + strand gene ({len(df_pos)} rows) -> {pos_path}")

    # Negative strand: TP53 (ENSG00000141510)
    df_neg = fetch_gene_dataframe("ENSG00000141510")
    neg_path = out_dir / "ENSG00000141510_exons.tsv"
    df_neg.to_csv(neg_path, sep="\t", index=False)
    print(f"Saved - strand gene ({len(df_neg)} rows) -> {neg_path}")
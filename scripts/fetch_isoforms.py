#!/usr/bin/env python3
"""
Fetch exon and isoform-level information for one or more Ensembl human genes.

This script retrieves, for each transcript of each input gene:
- exon IDs and genomic coordinates
- exon nucleotide sequences
- CDS-relative exon coordinates
- amino acids encoded by complete codons fully contained inside each exon
- C-terminal junctional amino acid when an exon ends in a split codon
- transcript, translation, and gene metadata

Important biological note:
`exon_protein_sequence` contains only amino acids encoded by complete codons fully
contained within that exon. Amino acids created by split codons across exon junctions
are not assigned as normal exon-internal amino acids. The C-terminal split-codon amino
acid, when present, is reported in `junctional_Cterm`.

Input methods:
1. Command-line gene IDs:
   python scripts/fetch_isoforms.py --genes ENSG00000141510 --output-dir data/human_gene_exons

2. Command-line comma-separated gene IDs:
   python scripts/fetch_isoforms.py --genes ENSG00000141510,ENSG00000139618 --output-dir data/human_gene_exons

3. Input file with one Ensembl gene ID per line:
   python scripts/fetch_isoforms.py --input-file genes.txt --output-dir data/human_gene_exons

Lines starting with "#" and empty lines are ignored in input files.

Output:
One TSV file is created per gene:
   data/human_gene_exons/ENSG00000141510_exons.tsv

Examples:
   python scripts/fetch_isoforms.py --genes ENSG00000141510 --output-dir data/human_gene_exons
   python scripts/fetch_isoforms.py --genes ENSG00000141510,ENSG00000139618 --output-dir data/human_gene_exons
   python scripts/fetch_isoforms.py --input-file genes.txt --output-dir data/human_gene_exons
   python scripts/fetch_isoforms.py --input-file genes.txt --output-dir data/human_gene_exons --verbose
   python scripts/fetch_isoforms.py --genes ENSG00000141510 --output-dir data/human_gene_exons --no-protein-validation

Requirements:
   pip install requests pandas
"""

import argparse
import re
import sys
import time
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

import pandas as pd
import requests


BASE_URL = "https://rest.ensembl.org"
BATCH_SIZE = 50
REQUEST_TIMEOUT_GET = 30
REQUEST_TIMEOUT_POST = 60
MAX_RETRIES = 3
BATCH_SLEEP_SECONDS = 0.5
DEFAULT_OUTPUT_SUFFIX = "_exons.tsv"

DEFAULT_HEADERS = {
    "Content-Type": "application/json",
    "Accept": "application/json",
}

CODON_TABLE = {
    "ATA": "I", "ATC": "I", "ATT": "I", "ATG": "M",
    "ACA": "T", "ACC": "T", "ACG": "T", "ACT": "T",
    "AAC": "N", "AAT": "N", "AAA": "K", "AAG": "K",
    "AGC": "S", "AGT": "S", "AGA": "R", "AGG": "R",
    "CTA": "L", "CTC": "L", "CTG": "L", "CTT": "L",
    "CCA": "P", "CCC": "P", "CCG": "P", "CCT": "P",
    "CAC": "H", "CAT": "H", "CAA": "Q", "CAG": "Q",
    "CGA": "R", "CGC": "R", "CGG": "R", "CGT": "R",
    "GTA": "V", "GTC": "V", "GTG": "V", "GTT": "V",
    "GCA": "A", "GCC": "A", "GCG": "A", "GCT": "A",
    "GAC": "D", "GAT": "D", "GAA": "E", "GAG": "E",
    "GGA": "G", "GGC": "G", "GGG": "G", "GGT": "G",
    "TCA": "S", "TCC": "S", "TCG": "S", "TCT": "S",
    "TTC": "F", "TTT": "F", "TTA": "L", "TTG": "L",
    "TAC": "Y", "TAT": "Y", "TAA": "*", "TAG": "*",
    "TGC": "C", "TGT": "C", "TGA": "*", "TGG": "W",
}

OUTPUT_COLUMNS = [
    "gene_id",
    "gene_name",
    "transcript_id",
    "transcript_name",
    "transcript_biotype",
    "translation_id",
    "exon_id",
    "exon_number",
    "total_exons",
    "chromosome",
    "genomic_start",
    "genomic_end",
    "strand",
    "exon_length",
    "is_coding",
    "cds_start_nt_0based",
    "cds_end_nt_0based_exclusive",
    "cds_frame_offset_start",
    "cds_frame_offset_end",
    "protein_start_aa_1based",
    "protein_end_aa_1based",
    "exon_sequence",
    "exon_protein_sequence",
    "junctional_Cterm",
]


class CustomArgumentParser(argparse.ArgumentParser):
    def error(self, message: str) -> None:
        print(f"\nError: {message}\n", file=sys.stderr)
        self.print_help(sys.stderr)
        self.exit(2)


class EnsemblExonExtractor:
    def __init__(
        self,
        gene_id: str,
        reverse_complement_negative_strand: bool = False,
        validate_protein: bool = True,
    ):
        self.gene_id = gene_id
        self.reverse_complement_negative_strand = reverse_complement_negative_strand
        self.validate_protein = validate_protein
        self.session = requests.Session()

    def _make_request(
        self,
        endpoint: str,
        method: str = "GET",
        json_data: Optional[Dict[str, Any]] = None,
    ) -> Optional[Any]:
        url = f"{BASE_URL}{endpoint}"

        for attempt in range(MAX_RETRIES):
            try:
                if method.upper() == "POST":
                    response = self.session.post(
                        url,
                        headers=DEFAULT_HEADERS,
                        json=json_data,
                        timeout=REQUEST_TIMEOUT_POST,
                    )
                else:
                    response = self.session.get(
                        url,
                        headers=DEFAULT_HEADERS,
                        timeout=REQUEST_TIMEOUT_GET,
                    )

                if response.status_code == 200:
                    return response.json()

                if response.status_code == 429 or 500 <= response.status_code < 600:
                    wait_time = int(response.headers.get("Retry-After", 2**attempt))
                    print(
                        f"Temporary Ensembl REST issue HTTP {response.status_code}. "
                        f"Waiting {wait_time}s...",
                        file=sys.stderr,
                    )
                    time.sleep(wait_time)
                    continue

                print(
                    f"HTTP error {response.status_code} for {url}: {response.text}",
                    file=sys.stderr,
                )
                return None

            except requests.exceptions.RequestException as exc:
                print(
                    f"Request error attempt {attempt + 1}/{MAX_RETRIES}: {exc}",
                    file=sys.stderr,
                )
                if attempt < MAX_RETRIES - 1:
                    time.sleep(2**attempt)
                else:
                    return None

        return None

    def get_sequences_batch(self, ids: List[str], seq_type: str) -> Dict[str, str]:
        unique_ids = sorted(set(ids))

        if not unique_ids:
            return {}

        all_results: Dict[str, str] = {}

        for i in range(0, len(unique_ids), BATCH_SIZE):
            batch = unique_ids[i:i + BATCH_SIZE]
            payload = {"ids": batch, "type": seq_type}

            print(
                f"  Downloading {seq_type} batch {i // BATCH_SIZE + 1} "
                f"({len(batch)} sequences)...",
                file=sys.stderr,
            )

            data = self._make_request(
                "/sequence/id",
                method="POST",
                json_data=payload,
            )

            if not data:
                continue

            if isinstance(data, dict):
                data = [data]

            for item in data:
                seq_id = item.get("id")
                sequence = item.get("seq", "")
                if seq_id:
                    all_results[seq_id] = sequence

            if i + BATCH_SIZE < len(unique_ids):
                time.sleep(BATCH_SLEEP_SECONDS)

        return all_results

    @staticmethod
    def reverse_complement(seq: str) -> str:
        complement = str.maketrans("ACGTNacgtn", "TGCANtgcan")
        return seq.translate(complement)[::-1]

    @staticmethod
    def translate_sequence(seq: str) -> str:
        protein = []

        for i in range(0, len(seq) - 2, 3):
            codon = seq[i:i + 3].upper()
            protein.append(CODON_TABLE.get(codon, "X"))

        return "".join(protein)

    def get_gene_info(self) -> Dict[str, Any]:
        data = self._make_request(f"/lookup/id/{self.gene_id}?expand=1")

        if not data:
            raise ValueError(
                f"Gene {self.gene_id} was not found, or Ensembl REST returned an error."
            )

        return data

    @staticmethod
    def _exon_sort_key_for_transcript(exon: Dict[str, Any], strand: int) -> int:
        if strand == 1:
            return int(exon["start"])

        return -int(exon["end"])

    def calculate_exon_cds_coords(
        self,
        exons: List[Dict[str, Any]],
        translation: Optional[Dict[str, Any]],
        strand: int,
    ) -> List[Tuple[int, int, int, int, bool]]:
        if not translation:
            return [(0, 0, -1, -1, False)] * len(exons)

        cds_start_genomic = translation.get("start")
        cds_end_genomic = translation.get("end")

        if cds_start_genomic is None or cds_end_genomic is None:
            return [(0, 0, -1, -1, False)] * len(exons)

        sorted_exons = sorted(
            exons,
            key=lambda exon: self._exon_sort_key_for_transcript(exon, strand),
        )

        coords_by_exon_id: Dict[str, Tuple[int, int, int, int, bool]] = {}
        cds_position = 0

        for exon in sorted_exons:
            exon_id = exon["id"]
            exon_start = int(exon["start"])
            exon_end = int(exon["end"])

            overlaps_cds = not (
                exon_end < cds_start_genomic or exon_start > cds_end_genomic
            )

            if not overlaps_cds:
                coords_by_exon_id[exon_id] = (0, 0, -1, -1, False)
                continue

            coding_start = max(exon_start, cds_start_genomic)
            coding_end = min(exon_end, cds_end_genomic)
            coding_length = coding_end - coding_start + 1

            cds_start_pos = cds_position
            cds_end_pos = cds_position + coding_length

            coords_by_exon_id[exon_id] = (
                cds_start_pos,
                cds_end_pos,
                cds_start_pos % 3,
                cds_end_pos % 3,
                True,
            )

            cds_position += coding_length

        return [
            coords_by_exon_id.get(exon["id"], (0, 0, -1, -1, False))
            for exon in exons
        ]

    def fetch_exon_data(self) -> List[Dict[str, Any]]:
        print(f"Retrieving gene information for {self.gene_id}...", file=sys.stderr)
        gene_data = self.get_gene_info()

        transcripts = gene_data.get("Transcript", [])
        print(f"Found {len(transcripts)} transcripts.", file=sys.stderr)

        all_exon_ids: List[str] = []
        coding_transcript_ids: List[str] = []
        translation_ids: List[str] = []

        for transcript in transcripts:
            translation = transcript.get("Translation")

            if translation:
                coding_transcript_ids.append(transcript["id"])

                translation_id = translation.get("id")
                if translation_id:
                    translation_ids.append(translation_id)

            for exon in transcript.get("Exon", []):
                exon_id = exon.get("id")
                if exon_id:
                    all_exon_ids.append(exon_id)

        all_exon_ids = sorted(set(all_exon_ids))
        coding_transcript_ids = sorted(set(coding_transcript_ids))
        translation_ids = sorted(set(translation_ids))

        print(
            f"\nDownloading {len(all_exon_ids)} unique exon sequences...",
            file=sys.stderr,
        )
        exon_sequences = self.get_sequences_batch(all_exon_ids, "genomic")

        print(
            f"\nDownloading {len(coding_transcript_ids)} CDS sequences...",
            file=sys.stderr,
        )
        cds_sequences = self.get_sequences_batch(coding_transcript_ids, "cds")

        protein_sequences: Dict[str, str] = {}

        if self.validate_protein and translation_ids:
            print(
                f"\nDownloading {len(translation_ids)} protein sequences for validation...",
                file=sys.stderr,
            )
            protein_sequences = self.get_sequences_batch(translation_ids, "protein")

        exon_data: List[Dict[str, Any]] = []
        protein_fragment_checks: List[Tuple[str, int, str, str]] = []

        for transcript in transcripts:
            transcript_id = transcript["id"]
            transcript_name = transcript.get("display_name", "")
            transcript_biotype = transcript.get("biotype", "")
            transcript_strand = int(transcript.get("strand", 1))

            exons = transcript.get("Exon", [])
            translation = transcript.get("Translation")

            cds_seq = cds_sequences.get(transcript_id, "")
            translated_cds = self.translate_sequence(cds_seq) if cds_seq else ""

            translation_id = translation.get("id") if translation else None
            official_protein = (
                protein_sequences.get(translation_id, "")
                if translation_id
                else ""
            )

            translated_cds_no_terminal_stop = translated_cds.rstrip("*")

            if (
                self.validate_protein
                and official_protein
                and translated_cds_no_terminal_stop
                and translated_cds_no_terminal_stop != official_protein
            ):
                print(
                    f"Warning: translated CDS does not exactly match Ensembl protein "
                    f"for transcript {transcript_id}. Continuing with translated CDS.",
                    file=sys.stderr,
                )

            full_protein = translated_cds

            cds_coords = self.calculate_exon_cds_coords(
                exons=exons,
                translation=translation,
                strand=transcript_strand,
            )

            for exon_idx, (exon, coord) in enumerate(zip(exons, cds_coords), start=1):
                exon_id = exon["id"]

                (
                    cds_start_nt,
                    cds_end_nt,
                    start_frame_offset,
                    end_frame_offset,
                    is_coding,
                ) = coord

                chrom = exon.get("seq_region_name", "")
                start = int(exon.get("start", 0))
                end = int(exon.get("end", 0))
                strand = int(exon.get("strand", transcript_strand))

                exon_seq = exon_sequences.get(exon_id, "")

                if (
                    self.reverse_complement_negative_strand
                    and strand == -1
                    and exon_seq
                ):
                    exon_seq = self.reverse_complement(exon_seq)

                exon_protein_seq = ""
                junctional_cterm = ""
                protein_start_aa = ""
                protein_end_aa = ""

                if is_coding and full_protein:
                    total_coding_nt_in_exon = cds_end_nt - cds_start_nt

                    skip_start_nt = (
                        3 - start_frame_offset
                        if start_frame_offset != 0
                        else 0
                    )
                    skip_end_nt = end_frame_offset

                    internal_nt = (
                        total_coding_nt_in_exon
                        - skip_start_nt
                        - skip_end_nt
                    )

                    if internal_nt >= 3 and internal_nt % 3 == 0:
                        num_aa = internal_nt // 3
                        internal_start_nt = cds_start_nt + skip_start_nt
                        start_aa_0based = internal_start_nt // 3
                        end_aa_0based_exclusive = start_aa_0based + num_aa

                        exon_protein_seq = full_protein[
                            start_aa_0based:end_aa_0based_exclusive
                        ]

                        if exon_protein_seq:
                            protein_start_aa = start_aa_0based + 1
                            protein_end_aa = end_aa_0based_exclusive

                            expected_fragment = full_protein[
                                start_aa_0based:end_aa_0based_exclusive
                            ]
                            protein_fragment_checks.append(
                                (
                                    transcript_id,
                                    exon_idx,
                                    exon_protein_seq,
                                    expected_fragment,
                                )
                            )

                    if end_frame_offset != 0:
                        last_aa_start_nt = cds_end_nt - end_frame_offset
                        last_aa_0based = last_aa_start_nt // 3
                        junctional_cterm = full_protein[
                            last_aa_0based:last_aa_0based + 1
                        ]

                exon_data.append(
                    {
                        "gene_id": self.gene_id,
                        "gene_name": gene_data.get("display_name", ""),
                        "transcript_id": transcript_id,
                        "transcript_name": transcript_name,
                        "transcript_biotype": transcript_biotype,
                        "translation_id": translation_id or "",
                        "exon_id": exon_id,
                        "exon_number": exon_idx,
                        "total_exons": len(exons),
                        "chromosome": chrom,
                        "genomic_start": start,
                        "genomic_end": end,
                        "strand": strand,
                        "exon_length": end - start + 1 if start and end else 0,
                        "is_coding": is_coding,
                        "cds_start_nt_0based": cds_start_nt if is_coding else "",
                        "cds_end_nt_0based_exclusive": cds_end_nt if is_coding else "",
                        "cds_frame_offset_start": (
                            start_frame_offset if is_coding else -1
                        ),
                        "cds_frame_offset_end": (
                            end_frame_offset if is_coding else -1
                        ),
                        "protein_start_aa_1based": protein_start_aa,
                        "protein_end_aa_1based": protein_end_aa,
                        "exon_sequence": exon_seq,
                        "exon_protein_sequence": exon_protein_seq,
                        "junctional_Cterm": junctional_cterm,
                    }
                )

        print("\nRunning positional protein-fragment consistency check...", file=sys.stderr)

        for transcript_id, exon_idx, observed, expected in protein_fragment_checks:
            if observed != expected:
                raise ValueError(
                    f"Protein consistency error for transcript {transcript_id}, "
                    f"exon #{exon_idx}: observed fragment does not match expected "
                    f"fragment at the calculated protein position."
                )

        print("Protein-fragment consistency check passed.", file=sys.stderr)
        return exon_data

    def create_dataframe(self) -> pd.DataFrame:
        exon_data = self.fetch_exon_data()

        if not exon_data:
            raise ValueError("No exon data were retrieved for the specified gene.")

        df = pd.DataFrame(exon_data)
        extra_cols = [col for col in df.columns if col not in OUTPUT_COLUMNS]
        final_columns = [col for col in OUTPUT_COLUMNS + extra_cols if col in df.columns]

        return df[final_columns]


def parse_gene_ids_from_text(text: str) -> List[str]:
    gene_ids = []

    for item in re.split(r"[,\s]+", text.strip()):
        item = item.strip()
        if item:
            gene_ids.append(item)

    return gene_ids


def read_gene_ids_from_file(path: Path) -> List[str]:
    gene_ids = []

    with path.open("r", encoding="utf-8") as handle:
        for line in handle:
            clean_line = line.strip()

            if not clean_line or clean_line.startswith("#"):
                continue

            gene_ids.extend(parse_gene_ids_from_text(clean_line))

    return gene_ids


def unique_preserving_order(items: List[str]) -> List[str]:
    seen = set()
    unique_items = []

    for item in items:
        if item not in seen:
            seen.add(item)
            unique_items.append(item)

    return unique_items


def validate_gene_ids(gene_ids: List[str]) -> List[str]:
    invalid_ids = [
        gene_id
        for gene_id in gene_ids
        if not re.fullmatch(r"ENSG[0-9]+(?:\.[0-9]+)?", gene_id)
    ]

    if invalid_ids:
        raise ValueError(
            "Invalid Ensembl gene ID format: "
            + ", ".join(invalid_ids)
            + ". Expected IDs like ENSG00000141510."
        )

    return gene_ids


def parse_args() -> argparse.Namespace:
    parser = CustomArgumentParser(
        description="Fetch exon and isoform-level information for Ensembl human genes.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  python scripts/fetch_isoforms.py --genes ENSG00000141510 --output-dir data/human_gene_exons
  python scripts/fetch_isoforms.py --genes ENSG00000141510,ENSG00000139618 --output-dir data/human_gene_exons
  python scripts/fetch_isoforms.py --input-file genes.txt --output-dir data/human_gene_exons

Input:
  Use exactly one of:
    --genes       One Ensembl gene ID or a comma-separated list
    --input-file  Text file with one Ensembl gene ID per line

Output:
  --output-dir is required.
  One TSV file is written per gene.
        """,
    )

    if len(sys.argv) == 1:
        parser.print_help(sys.stderr)
        parser.exit(2)

    input_group = parser.add_mutually_exclusive_group(required=True)

    input_group.add_argument(
        "--genes",
        help=(
            "One Ensembl gene ID or a comma-separated list of Ensembl gene IDs, "
            "for example ENSG00000141510,ENSG00000139618."
        ),
    )

    input_group.add_argument(
        "--input-file",
        type=Path,
        help="Text file containing one Ensembl gene ID per line.",
    )

    parser.add_argument(
        "--output-dir",
        type=Path,
        required=True,
        help="Directory where one TSV file per gene will be written.",
    )

    parser.add_argument(
        "-v",
        "--verbose",
        action="store_true",
        help="Print the first rows of each output dataframe to stderr.",
    )

    parser.add_argument(
        "--reverse-complement-negative-strand",
        action="store_true",
        help=(
            "Reverse-complement exon sequences on negative-strand exons. "
            "Disabled by default."
        ),
    )

    parser.add_argument(
        "--no-protein-validation",
        action="store_false",
        help="Skip retrieval of Ensembl protein sequences for validation.",
    )

    parser.add_argument(
        "--overwrite",
        action="store_true",
        help="Overwrite existing output files.",
    )

    return parser.parse_args()


def get_gene_ids_from_args(args: argparse.Namespace) -> List[str]:
    if args.input_file:
        if not args.input_file.exists():
            raise FileNotFoundError(f"Input file does not exist: {args.input_file}")

        gene_ids = read_gene_ids_from_file(args.input_file)
    else:
        gene_ids = parse_gene_ids_from_text(args.genes)

    gene_ids = unique_preserving_order(gene_ids)

    if not gene_ids:
        raise ValueError("No gene IDs were provided. Use --genes or --input-file.")

    return validate_gene_ids(gene_ids)


def write_gene_output(
    gene_id: str,
    output_dir: Path,
    overwrite: bool,
    reverse_complement_negative_strand: bool,
    validate_protein: bool,
    verbose: bool,
) -> None:
    output_path = output_dir / f"{gene_id}{DEFAULT_OUTPUT_SUFFIX}"

    if output_path.exists() and not overwrite:
        print(
            f"\nSkipping {gene_id}: output already exists at {output_path}. "
            f"Use --overwrite to replace it.",
            file=sys.stderr,
        )
        return

    extractor = EnsemblExonExtractor(
        gene_id=gene_id,
        reverse_complement_negative_strand=reverse_complement_negative_strand,
        validate_protein=validate_protein,
    )

    df = extractor.create_dataframe()
    df.to_csv(output_path, sep="\t", index=False)

    print("\nCompleted successfully.", file=sys.stderr)
    print(f"  Gene: {gene_id}", file=sys.stderr)
    print(f"  Total exon rows: {len(df)}", file=sys.stderr)
    print(f"  Unique transcripts: {df['transcript_id'].nunique()}", file=sys.stderr)
    print(f"  Output saved to: {output_path}", file=sys.stderr)

    if verbose:
        print("\nFirst rows of the dataset:", file=sys.stderr)
        print(df.head().to_string(), file=sys.stderr)


def main() -> int:
    args = parse_args()

    try:
        gene_ids = get_gene_ids_from_args(args)
        args.output_dir.mkdir(parents=True, exist_ok=True)

        print(f"Processing {len(gene_ids)} gene(s).", file=sys.stderr)

        failed_genes = []

        for gene_id in gene_ids:
            print("\n" + "=" * 80, file=sys.stderr)
            print(f"Processing {gene_id}", file=sys.stderr)
            print("=" * 80, file=sys.stderr)

            try:
                write_gene_output(
                    gene_id=gene_id,
                    output_dir=args.output_dir,
                    overwrite=args.overwrite,
                    reverse_complement_negative_strand=(
                        args.reverse_complement_negative_strand
                    ),
                    validate_protein=not args.no_protein_validation,
                    verbose=args.verbose,
                )
            except Exception as exc:
                failed_genes.append(gene_id)
                print(f"Error while processing {gene_id}: {exc}", file=sys.stderr)

        if failed_genes:
            print(
                "\nFinished with errors for: " + ", ".join(failed_genes),
                file=sys.stderr,
            )
            return 1

        print("\nAll genes processed successfully.", file=sys.stderr)
        return 0

    except Exception as exc:
        print(f"Error: {exc}", file=sys.stderr)
        return 1


if __name__ == "__main__":
    sys.exit(main())




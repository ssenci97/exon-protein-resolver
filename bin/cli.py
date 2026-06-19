"""Command-line argument parsing."""

import argparse
import re
import sys
from pathlib import Path
from typing import List


class CustomArgumentParser(argparse.ArgumentParser):
    def error(self, message: str) -> None:
        print(f"\nError: {message}\n", file=sys.stderr)
        self.print_help(sys.stderr)
        self.exit(2)


def parse_gene_ids_from_text(text: str) -> List[str]:
    return [item.strip() for item in re.split(r"[,\s]+", text.strip()) if item.strip()]


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
    return [item for item in items if not (item in seen or seen.add(item))]


def validate_gene_ids(gene_ids: List[str]) -> List[str]:
    invalid = [gid for gid in gene_ids if not re.fullmatch(r"ENSG[0-9]+(?:\.[0-9]+)?", gid)]
    if invalid:
        raise ValueError(
            f"Invalid Ensembl gene ID format: {', '.join(invalid)}. "
            "Expected IDs like ENSG00000141510."
        )
    return gene_ids


def parse_args() -> argparse.Namespace:
    parser = CustomArgumentParser(
        description="Fetch exon and isoform-level information for Ensembl human genes.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Input:
  Use exactly one of:
    --genes       One Ensembl gene ID or a comma-separated list
    --input-file  Text file with one Ensembl gene ID per line

Output:
  --output-dir is required. One TSV file is written per gene.
        """,
    )

    if len(sys.argv) == 1:
        parser.print_help(sys.stderr)
        parser.exit(2)

    input_group = parser.add_mutually_exclusive_group(required=True)
    input_group.add_argument(
        "--genes",
        help="Comma-separated Ensembl gene IDs, e.g. ENSG00000141510,ENSG00000139618.",
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
        help="Reverse-complement exon sequences on negative-strand exons.",
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


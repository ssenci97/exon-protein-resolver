"""CDS coordinate and exon ordering calculations."""

from typing import Any, Dict, List, Optional, Tuple


def exon_sort_key_for_transcript(exon: Dict[str, Any], strand: int) -> int:
    if strand == 1:
        return int(exon["start"])
    return -int(exon["end"])


def sort_exons_by_transcript_order(
    exons: List[Dict[str, Any]], strand: int
) -> List[Dict[str, Any]]:
    """Return exons sorted in transcript order (5'->3') for the given strand."""
    return sorted(exons, key=lambda e: exon_sort_key_for_transcript(e, strand))


def calculate_exon_cds_coords(
    exons: List[Dict[str, Any]],
    translation: Optional[Dict[str, Any]],
    strand: int,
) -> List[Tuple[int, int, int, int, bool]]:
    """Calculate CDS-relative coordinates for each exon.

    Exons must already be sorted in transcript order (5'->3') on entry.
    The returned list is parallel to the input exon list.

    Returns per-exon tuples of:
      (cds_start_nt_0based, cds_end_nt_0based_exclusive,
       frame_offset_start, frame_offset_end, is_coding)

    frame_offset_start: bases into the first codon at the 5' edge (0-2).
    frame_offset_end:   bases of the last split codon extending into the
                        next exon at the 3' edge (0-2). Always 0 for the
                        final coding exon.
    """
    if not translation:
        return [(0, 0, -1, -1, False)] * len(exons)

    cds_start_genomic = translation.get("start")
    cds_end_genomic = translation.get("end")
    if cds_start_genomic is None or cds_end_genomic is None:
        return [(0, 0, -1, -1, False)] * len(exons)

    # --- find last coding exon in transcript order ---
    last_coding_exon_id: Optional[str] = None
    for exon in exons:
        exon_start = int(exon["start"])
        exon_end = int(exon["end"])
        overlaps = not (exon_end < cds_start_genomic or exon_start > cds_end_genomic)
        if overlaps:
            last_coding_exon_id = exon["id"]

    coords_by_exon_id: Dict[str, Tuple[int, int, int, int, bool]] = {}
    cds_position = 0

    for exon in exons:
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

        frame_offset_end = cds_end_pos % 3 if exon_id != last_coding_exon_id else 0

        coords_by_exon_id[exon_id] = (
            cds_start_pos,
            cds_end_pos,
            cds_start_pos % 3,
            frame_offset_end,
            True,
        )
        cds_position += coding_length

    return [
        coords_by_exon_id.get(exon["id"], (0, 0, -1, -1, False))
        for exon in exons
    ]


#!/usr/bin/env python3

"""Library for processing FASTA files.

Functions
---------
rename_header: Rename a header.
prefix_to_sequence_ids: Prefix to sequence ids.
slice_records_by_exact_ids: Slice records by exact match of sequence ids.
slice_records_by_partial_ids: Slice records by partial match of sequence ids.

"""

import re
from Bio import SeqIO

def rename_header(input_filename: str, output_filename: str, output_id: str, output_name: str, output_description: str) -> None:
    """Rename a header.

    Args
    ----
    input_filename : str
        Input filename.
    output_filename : str
        Output filename.
    output_id : str
        Sequence id.
    output_name : str
        Sequence name.
    output_description : str
        Sequence description.
    
    """
    with open(input_filename, mode="r") as input_handle:
        record = SeqIO.read(input_handle, "fasta")
        record.id = output_id
        record.name = output_name
        record.description = output_description
    with open(output_filename, mode="w") as output_handle:
        SeqIO.write(record, output_handle, "fasta")


def prefix_to_sequence_ids(input_filename: str, output_filename: str, prefix: str) -> None:
    """
    Prefix to sequence ids.

    Args
    ----
    input_filename : str
        Input filename.
    output_filename : str
        Output filename.
    prefix : str
        String to prefix.

    """
    with open(input_filename, "r") as input_handle, open(output_filename, "w") as output_handle:
        for record in SeqIO.parse(input_handle, "fasta"):
            record.id = f"{prefix}_{record.id}"
            record.name = ""
            record.description = ""
            SeqIO.write(record, output_handle, "fasta")


def slice_records_by_exact_ids(input_filename: str, output_filename: str, *input_ids: str) -> None:
    """Slice records by exact match of sequence ids.

    Args
    ----
    input_filename : str
        Input filename.
    output_filename : str
        Output filename.
    input_ids : tuple
        Sequence ids to slice records.

    """
    pattern = "|".join(r"^" + re.escape(input_id) + r"$" for input_id in input_ids)

    with open(input_filename, "r") as input_handle, open(output_filename, "w") as output_handle:
        for record in SeqIO.parse(input_handle, "fasta"):
            if not re.fullmatch(pattern, record.id):
                continue
            SeqIO.write(record, output_handle, "fasta")


def slice_records_by_partial_ids(input_filename: str, output_filename: str, *input_ids: str) -> None:
    """Slice records by partial match of sequence ids.

    Args
    ----
    input_filename : str
        Input filename.
    output_filename : str
        Output filename.
    input_ids : tuple
        Sequence ids to slice records.

    """
    pattern = "|".join(r"\b" + re.escape(input_id) + r"\b" for input_id in input_ids)

    with open(input_filename, "r") as input_handle, open(output_filename, "w") as output_handle:
        for record in SeqIO.parse(input_handle, "fasta"):
            if not re.search(pattern, record.id):
                continue
            SeqIO.write(record, output_handle, "fasta")
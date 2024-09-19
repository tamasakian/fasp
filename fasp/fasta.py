#!/usr/bin/env python3

"""Library for processing FASTA files.

Functions
---------
rename_header: Rename a header.
prefix_to_sequence_ids: Prefix to sequence ids.
split_multi_to_single: Split a multi FASTA file into individual single sequence FASTA files.
merge_msa_by_ids: Merge MSAs by sequence ids.
slice_records_by_exact_ids: Slice records by exact match of sequence ids.
slice_records_by_partial_ids: Slice records by partial match of sequence ids.
measure_lengths: Measure sequence lengths.

"""

import re
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

def rename_header(input_filename: str, output_filename: str, output_id: str, output_description: str) -> None:
    """Rename a header.

    Args
    ----
    input_filename : str
        Input filename.
    output_filename : str
        Output filename.
    output_id : str
        Sequence id.
    output_description : str
        Sequence description.
    
    """
    with open(input_filename, mode="r") as input_handle:
        record = SeqIO.read(input_handle, "fasta")
        record.id = output_id
        record.description = output_description
    with open(output_filename, mode="w") as output_handle:
        SeqIO.write(record, output_handle, "fasta")


def prefix_to_sequence_ids(input_filename: str, output_filename: str, prefix: str) -> None:
    """Prefix to sequence ids.

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


def split_multi_to_single(input_filename: str, output_dirname: str) -> None:
    """Split a multi FASTA file into individual single sequence FASTA files.

    This function reads a multi-FASTA file and creates a separate FASTA file 
    for each sequence in the input file. 
    
    The output files are named based on the sequence IDs and 
    saved in the specified output directory.

    Args
    ----
    input_filename : str
        Input filename.
    output_dirname : str
        Output directory where the individual FASTA files will be saved.

    """
    with open(input_filename, "r") as input_handle:
        for record in SeqIO.parse(input_handle, "fasta"):
            output_filename = f"{output_dirname}/{record.id}.fasta"
            SeqIO.write(record, output_filename, "fasta")


def merge_msa_by_ids(input_filename: str, output_filename: str) -> None:
    """Merge MSAs by sequence ids.

    This function reads a multi-FASTA file and concatenates sequences 
    with the same sequence id. 
    
    The merged sequences are then saved to a new FASTA file.

    Args
    ----
    input_filename : str
        Input filename.
    output_filename : str
       Output filename where the merged sequences will be saved.

    """
    msa = {}

    with open(input_filename, "r") as input_handle:
        for record in SeqIO.parse(input_handle, "fasta"):
            if record.id not in msa:
                msa[record.id] = ""
            msa[record.id] += str(record.seq)

    with open(output_filename, "w") as output_handle:
        records = []
        for index, sequence in msa.items():
            record = SeqRecord(Seq(sequence), id=index, description="")
            records.append(record)
        SeqIO.write(records, output_handle, "fasta")


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


def measure_lengths(input_filename: str, output_filename: str) -> None:
    """Measure sequence lengths.

    Args
    ----
    input_filename : str
        Input filename.
    output_filename : str
        Output filename (TSV format).

    """
    with open(input_filename, "r") as input_handle, open(output_filename, "w") as output_handle:
        output_handle.write("name\tlength\n")
        for record in SeqIO.parse(input_handle, "fasta"):
            name = record.id
            length = len(record.seq)
            output_handle.write(f"{name}\t{length}\n")

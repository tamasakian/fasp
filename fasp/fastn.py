#!/usr/bin/env python3

"""Library for processing nucleotide FASTA files.

Functions
---------
rename_headers_feature: Rename headers feature values.
slice_records_by_seqids: Slice records by match of seqids.
slice_sequence_by_flanking_region: Slice a sequence with the flanking region and display start and stop codons.
slice_sequence_by_upstream_region: Slice a sequence with the upstream region and display start codon.

"""

import re
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

def rename_headers_feature(input_filename: str, output_filename: str, feature: str) -> None:
    """Rename headers feature values.

    Args
    ----
    input_filename : str
        Input filename.
    output_filename : str
        Output filename.
    feature : str
        Feature to rename headers. 
        e.g., locus_tag, protein_id.

    """
    pattern = re.compile(r"\[" + re.escape(feature) + r"=([^\]]+)]")

    def extract_feature_value(text: str, default: str = "") -> str:
        """Extract feature value from text if it matches the pattern.

        Args
        ----
        text : str
            Text to extract feature value.
        default : str
            Default value which is returned if the feature is not found. 
        
        Returns
        -------
        match.group(1) : str
            Feature value if the feature is found.
            If this is not found, default is returned instead.

        """
        match = pattern.search(text)
        return match.group(1) if match else default

    with open(input_filename, "r") as input_handle, open(output_filename, "w") as output_handle:
        for record in SeqIO.parse(input_handle, "fasta"):
            feature_value = (extract_feature_value(record.description) or
                             extract_feature_value(record.name) or
                             extract_feature_value(record.id, default=record.id))
            record.id = feature_value
            record.name = ""
            record.description = ""  
            SeqIO.write(record, output_handle, "fasta")


def slice_records_by_seqids(input_filename: str, output_filename: str, *input_seqids: str) -> None:
    """Slice records by match of seqids.

    Args
    ----
    input_filename : str
        Input filename.
    output_filename : str
        Output filename.
    input_seqids : tuple
        Seqids to slice records.

    Example
    -------
    If record.id is 'lcl|ABC123_cds_...':
    The function extracts 'ABC123' and checks if it is in input_seqids.

    """
    with open(input_filename, "r") as input_handle, open(output_filename, "w") as output_handle:
        for record in SeqIO.parse(input_handle, "fasta"):
            seqid = record.id.split("lcl|")[1].split("_cds_")[0]
            if seqid not in input_seqids:
                continue
            SeqIO.write(record, output_handle, "fasta")


def slice_sequence_by_flanking_region(
        input_filename: str, 
        output_filename: str, 
        output_id: str, 
        output_name: str, 
        output_description: str, 
        strand: str, 
        start: int, 
        end: int, 
        bp: int) -> None:
    """Slice a sequence with the flanking region and display start and stop codons.

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
    strand : str
        Sequence strand.
        '+' for forward, '-' for reverse.
    start : int
        Start position of the sequence to slice.
    end : int
        End position of the sequence to slice.
    bp : int
        Number of base pairs to include as the flanking regions.

    """
    def to_int(value: str) -> int:
        """Convert string to integer.
        
        Args
        ----
        value : str
            Value to convert to integer.

        Returns
        -------
        int(value) : int
            Integer value.

        Raises
        ------
        ValueError
            If value is not valid integers.

        """
        if not value.isdigit():
            raise ValueError(f"Invalid integer value: {value}")
        return int(value)

    start = to_int(start)
    end = to_int(end)
    bp = to_int(bp)

    with open(input_filename, mode="r") as input_handle:
        record = SeqIO.read(input_handle, "fasta")
    
    if strand == "+":
        sequence = record.seq[start-1-bp:end+bp]
    elif strand == "-":
        sequence = record.seq[start-1-bp:end+bp].reverse_complement()
    else:
        raise ValueError("Strand must be either '+' or '-'.")

    start_codon = sequence[bp:bp+3]
    stop_codon = sequence[-(bp+3):-bp]
    length = len(sequence)
    print(f"Start codon: {start_codon}, Stop codon: {stop_codon}, Length: {length}")

    with open(output_filename, mode="w") as output_handle:
        output_record = SeqRecord(
            seq=sequence,
            id=output_id,
            name=output_name,
            description=output_description
        )
        SeqIO.write(output_record, output_handle, "fasta")


def slice_sequence_by_upstream_region(
        input_filename: str,
        output_filename: str,
        output_id: str,
        output_name: str,
        output_description: str,
        strand: str,
        start: int,
        end: int,
        bp: int) -> None:
    """Slice a sequence with the upstream region and display start codon.

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
    strand : str
        Sequence strand.
        '+' for forward, '-' for reverse.
    start : int
        Start position of the sequence to slice.
    end : int
        End position of the sequence to slice.
    bp : int
        Number of base pairs to include as the upstream region.

    """
    def to_int(value: str) -> int:
        """Convert string to integer.
        
        Args
        ----
        value : str
            Value to convert to integer.

        Returns
        -------
        int(value) : int
            Integer value.

        Raises
        ------
        ValueError
            If value is not valid integers.

        """
        if not value.isdigit():
            raise ValueError(f"Invalid integer value: {value}")
        return int(value)

    start = to_int(start)
    end = to_int(end)
    bp = to_int(bp)

    with open(input_filename, mode="r") as input_handle:
        record = SeqIO.read(input_handle, "fasta")
    
    if strand == "+":
        sequence = record.seq[start-1-bp:start+2]
    elif strand == "-":
        sequence = record.seq[end-3:end+bp].reverse_complement()
    else:
        raise ValueError("Strand must be either '+' or '-'.")

    start_codon = sequence[bp:bp+3]
    length = len(sequence)
    print(f"Start codon: {start_codon}, Length: {length}")

    with open(output_filename, mode="w") as output_handle:
        output_record = SeqRecord(
            seq=sequence,
            id=output_id,
            name=output_name,
            description=output_description
        )
        SeqIO.write(output_record, output_handle, "fasta")


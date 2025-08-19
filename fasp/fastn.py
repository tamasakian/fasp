#!/usr/bin/env python3

"""Library for processing nucleotide FASTA files.

Functions
---------
rename_headers_feature: Rename headers feature values.
slice_records_by_seqids: Slice records by match of seqids.
slice_sequence_by_flanking_region: Slice a sequence with the flanking region and display start and stop codons.
slice_sequence_by_upstream_region: Slice a sequence with the upstream region and display start codon.
slice_sequence_by_downstream_region: Slice a sequence with the downstream region and display stop codon.
generate_introns: Generate introns for genes.
generate_upstream_regions: Generate upstream sequences for genes.
generate_downstream_regions: Generate downstream sequences for genes.

"""

import re
from Bio import SeqIO
from Bio.Seq import Seq
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
    start_position = int(start) - 1
    end_position = int(end) - 1 
    bp = int(bp)

    with open(input_filename, mode="r") as input_handle:
        record = SeqIO.read(input_handle, "fasta")
    
    if strand == "+":
        start_codon = record.seq[start_position:start_position+3]
        stop_codon = record.seq[end_position-2:end_position+1]
        flanking = record.seq[start_position-bp:end_position+1+bp]
        length = len(flanking)
        print(f"Start codon: {start_codon}, Stop codon: {stop_codon}, Length: {length}")
    elif strand == "-":
        start_codon = record.seq[end_position-2:end_position+1].reverse_complement()
        stop_codon = record.seq[start_position:start_position+3].reverse_complement()
        flanking = record.seq[start_position-bp:end_position+1+bp].reverse_complement()
        length = len(flanking)
        print(f"Start codon: {start_codon}, Stop codon: {stop_codon}, Length: {length}")
    else:
        raise ValueError("Strand must be either + or -.")

    with open(output_filename, mode="w") as output_handle:
        output_record = SeqRecord(seq=flanking, id=output_id, description=output_description)
        SeqIO.write(output_record, output_handle, "fasta")


def slice_sequence_by_upstream_region(
        input_filename: str,
        output_filename: str,
        output_id: str,
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
    start_position = int(start) - 1
    end_position = int(end) - 1
    bp = int(bp)

    with open(input_filename, mode="r") as input_handle:
        record = SeqIO.read(input_handle, "fasta")
    
    if strand == "+":
        start_codon = record.seq[start_position:start_position+3]
        upstream = record.seq[start_position-bp:start_position]
        length = len(upstream)
        print(f"Start codon: {start_codon}, Length: {length}")
    elif strand == "-":
        start_codon = record.seq[end_position-2:end_position+1].reverse_complement()
        upstream = record.seq[end_position+1:end_position+1+bp].reverse_complement()
        length = len(upstream)
        print(f"Start codon: {start_codon}, Length: {length}")
    else:
        raise ValueError("Strand must be either + or -.")

    with open(output_filename, mode="w") as output_handle:
        output_record = SeqRecord(seq=upstream, id=output_id, description=output_description)
        SeqIO.write(output_record, output_handle, "fasta")


def slice_sequence_by_downstream_region(
        input_filename: str,
        output_filename: str,
        output_id: str,
        output_description: str,
        strand: str,
        start: int,
        end: int,
        bp: int) -> None:
    """Slice a sequence with the downstream region and display stop codon.

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
    strand : str
        Sequence strand.
        '+' for forward, '-' for reverse.
    start : int
        Start position of the sequence to slice.
    end : int
        End position of the sequence to slice.
    bp : int
        Number of base pairs to include as the downstream region.

    """
    start_position = int(start) - 1
    end_position = int(end) - 1
    bp = int(bp)

    with open(input_filename, mode="r") as input_handle:
        record = SeqIO.read(input_handle, "fasta")

    if strand == "+":
        stop_codon = record.seq[end_position-2:end_position+1]
        downstream = record.seq[end_position+1:end_position+1+bp]
        length = len(downstream)
        print(f"Stop codon: {stop_codon}, Length: {length}")
    elif strand == "-":
        stop_codon = record.seq[start_position:start_position+3].reverse_complement()
        downstream = record.seq[start_position-bp:start_position].reverse_complement()
        length = len(downstream)
        print(f"Stop codon: {stop_codon}, Length: {length}")
    else:
        raise ValueError("Strand must be either + or -.")

    with open(output_filename, mode="w") as output_handle:
        output_record = SeqRecord(seq=downstream, id=output_id, description=output_description)
        SeqIO.write(output_record, output_handle, "fasta")

def generate_introns(input_filename: str, output_filename: str, gff_filename: str) -> None:
    """Generate introns.

    Args
    ----
    input_filename : str
        Input filename.
    output_filename : str
        Output filename.
    gff_filename : str
        Generated intron GFF filename.

    """
    introns = []
    with open(gff_filename, "r") as gff_handle:
        for line in gff_handle:
            if line.startswith("#"):
                continue
            li = line.strip().split("\t")
            if len(li) != 9:
                continue
            seqid, src, kind, start, end, score, strand, phase, attributes = li
            if kind == "intron":
                attr_dict = {}
                for attr in attributes.split(";"):
                    key, value = attr.split("=")
                    attr_dict[key] = value
                intron_id = attr_dict.get("intron_id")
                introns.append((seqid, int(start), int(end), strand, intron_id))
    
    genome = SeqIO.to_dict(SeqIO.parse(input_filename, "fasta"))

    with open(output_filename, "w") as output_handle:
        for intron in introns:
            seqid, start, end, strand, intron_id = intron
            if seqid in genome:
                intron_seq = genome[seqid].seq[start-1:end]
                if strand == "-":
                    intron_seq = intron_seq.reverse_complement()
                intron_seq = str(intron_seq)
                intron_seq = "\n".join([intron_seq[i:i + 80] for i in range(0, len(intron_seq), 80)])
                output_handle.write(f">{intron_id}\n{intron_seq}\n")


def generate_upstream_regions(input_filename: str, output_filename: str, gff_filename: str) -> None:
    """Generate upstream sequences for genes.

    Args
    ----
    input_filename : str
        Input FASTA filename.
    output_filename : str
        Output FASTA filename.
    gff_filename : str
        Generated upstream GFF filename.

    """
    upstreams = []
    with open(gff_filename, "r") as gff_handle:
        for line in gff_handle:
            if line.startswith("#"):
                continue
            li = line.strip().split("\t")
            if len(li) != 9:
                continue
            seqid, src, kind, start, end, score, strand, phase, attributes = li
            if kind == "upstream":
                attr_dict = {}
                for attr in attributes.split(";"):
                    key, value = attr.split("=")
                    attr_dict[key] = value
                upstream_id = attr_dict.get("upstream_id")
                upstreams.append((seqid, int(start), int(end), strand, upstream_id))
    
    genome = SeqIO.to_dict(SeqIO.parse(input_filename, "fasta"))

    with open(output_filename, "w") as output_handle:
        for upstream in upstreams:
            seqid, start, end, strand, upstream_id = upstream
            if seqid in genome:
                upstream_seq = genome[seqid].seq[start-1:end]
                if strand == "-":
                    upstream_seq = upstream_seq.reverse_complement()
                upstream_seq = str(upstream_seq)
                upstream_seq = "\n".join([upstream_seq[i:i + 80] for i in range(0, len(upstream_seq), 80)])
                output_handle.write(f">{upstream_id}\n{upstream_seq}\n")


def generate_downstream_regions(input_filename: str, output_filename: str, gff_filename: str) -> None:
    """Generate downstream sequences for genes.

    Args
    ----
    input_filename : str
        Input FASTA filename.
    output_filename : str
        Output FASTA filename.
    gff_filename : str
        Generated downstream GFF filename.

    """
    downstreams = []
    with open(gff_filename, "r") as gff_handle:
        for line in gff_handle:
            if line.startswith("#"):
                continue
            li = line.strip().split("\t")
            if len(li) != 9:
                continue
            seqid, src, kind, start, end, score, strand, phase, attributes = li
            if kind == "downstream":
                attr_dict = {}
                for attr in attributes.split(";"):
                    key, value = attr.split("=")
                    attr_dict[key] = value
                downstream_id = attr_dict.get("downstream_id")
                downstreams.append((seqid, int(start), int(end), strand, downstream_id))
    
    genome = SeqIO.to_dict(SeqIO.parse(input_filename, "fasta"))

    with open(output_filename, "w") as output_handle:
        for downstream in downstreams:
            seqid, start, end, strand, downstream_id = downstream
            if seqid in genome:
                downstream_seq = genome[seqid].seq[start-1:end]
                if strand == "-":
                    downstream_seq = downstream_seq.reverse_complement()
                downstream_seq = str(downstream_seq)
                downstream_seq = "\n".join([downstream_seq[i:i + 80] for i in range(0, len(downstream_seq), 80)])
                output_handle.write(f">{downstream_id}\n{downstream_seq}\n")
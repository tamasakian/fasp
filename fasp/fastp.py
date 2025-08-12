#!/usr/bin/env python3

"""Library for processing protein FASTA files.

Functions
---------
exclude_isoforms_by_length: Exclude isoforms based on length. 
exclude_non_nuclear_proteins: Exclude mitochondrial and chloroplast proteins. 
og_prefixer: Add ortholog group (OG) prefixes to the sequence IDs based on the OG table by SonicParanoid2. 

"""

from collections import defaultdict
from Bio import SeqIO

def exclude_isoforms_by_length(input_filename: str, output_filename: str, gff3_file: str) -> None:
    """Exclude isoforms based on length.
    
    Args
    ----
    input_filename : str
        Input protein FASTA filename.
    output_filename : str
        Output protein FASTA filename.
    gff3_file : str
        Input genome GFF3 filename.

    """
    
    def parse_gff3(gff3_file: str) -> dict:
        """Parse GFF3 file and make dict with 'protein_id', 'start', 'end' and 'length' of each gene. 

        Args
        ----
        gff3_file : str

        Returns
        -------
        genes : dict
            Dict with 'protein_id', 'start', 'end' and 'length' of each genes. 

        """
        genes = {}
        with open(gff3_file, "r") as gff3_handle:
            for line in gff3_handle:
                ## Exclude comments
                if line.startswith("#"):
                    continue

                li = line.strip().split("\t")
                if len(li) != 9:
                    continue
                seqid, src, kind, start, end, score, strand, phase, attributes = li

                ## Exclude lines other than CDS.
                if kind != "CDS":
                    continue

                ## Handle attributes.
                attr_dict = {}
                for attr in attributes.split(";"):
                    key_value = attr.split("=")
                    if len(key_value) != 2:
                        continue
                    key, value = key_value
                    attr_dict[key] = value
                
                ## Exclude CDS without protein_id.
                if "protein_id" not in attr_dict:
                    continue

                ## Exclude CDS without locus_tag and gene.
                if "locus_tag" not in attr_dict and "gene" not in attr_dict:
                    continue

                ## Read information of CDS.
                protein_id = attr_dict["protein_id"]
                if "locus_tag" in attr_dict:
                    gene = attr_dict["locus_tag"]
                elif "gene" in attr_dict:
                    gene = attr_dict["gene"]
                start, end, length = int(start), int(end), int(end) - int(start)
                if gene not in genes:
                    genes[gene] = []
                genes[gene].append({"protein_id": protein_id, "start": start, "end": end, "length": length})

        return genes

    def select_longest_protein(genes: dict) -> dict:
        """Select the longest proteins for each gene based on CDS information.

        Args
        ----
        genes : dict

        Returns
        -------
        longest_proteins : dict

        """
        longest_proteins = {}
        for gene, cds_list in genes.items():
            protein_lengths = {}
            for cds in cds_list:
                protein_id = cds["protein_id"]
                length = cds["length"]
                if protein_id not in protein_lengths:
                    protein_lengths[protein_id] = length
                else:
                    protein_lengths[protein_id] += length
            ## Select the longest protein.
            longest_proteins[gene] = max(protein_lengths, key=protein_lengths.get)

        return longest_proteins

    def slice_proteins(input_filename: str, output_filename: str, longest_proteins: dict) -> None:
        """Slice FASTA file to retain only the longest proteins for each gene.

        Args
        ----
        input_filename : str
        output_filename : str
        longest_proteins : dict

        """
        input_proteins = SeqIO.to_dict(SeqIO.parse(input_filename, "fasta"))
        selected_protein_ids = set(longest_proteins.values())

        output_proteins = []
        for selected_protein_id in selected_protein_ids:
            if selected_protein_id not in input_proteins:
                continue
            output_proteins.append(input_proteins[selected_protein_id])
        
        with open(output_filename, "w") as output_handle:
            SeqIO.write(output_proteins, output_handle, "fasta")

    genes = parse_gff3(gff3_file)
    longest_proteins = select_longest_protein(genes)
    slice_proteins(input_filename, output_filename, longest_proteins)


def exclude_non_nuclear_proteins(input_filename: str, output_filename: str) -> None:
    """Exclude mitochondrial and chloroplast proteins from the protein sequences.
    
    Args
    ----
    input_filename : str
        Input protein FASTA filename.
    output_filename : str
        Output protein FASTA filename.

    """
    with open(input_filename, mode="r") as input_handle, open(output_filename, "w") as output_handle:
        for protein in SeqIO.parse(input_handle, "fasta"):
            description = protein.description.lower()
            ## Exclude mitochondrion and chloroplast.
            if "(mitochondrion)" in description or "(chloroplast)" in description:
                continue
            SeqIO.write(protein, output_handle, "fasta")

def og_prefixer(input_filename: str, output_filename: str, og_filename: str) -> None:
    """
    Add ortholog group (OG) prefixes to the sequence IDs in a FASTA file based on the OG table by SonicPranoid2.

    Args
    ----
    input_filename : str
        Path to the input FASTA file. 
        This file should contain sequences with IDs that represent Protein IDs.

    output_filename : str
        Path to the output FASTA file where the updated sequences with OG-prefixed IDs will be saved. 

    og_filename : str
        Path to the OG table file in tab-delimited format.
        The table should include columns with gene IDs (starting from column 5) and OG names.

    """

    og_dict = defaultdict(set)

    with open(og_filename, "r") as og_handle:
        for line in og_handle:
            line = line.strip()
            if not line:
                continue
            li = line.split("\t")
            og_name = f"OG{li[0]}"

            for og_members in li[4:]:
                gene_list = og_members.split(",")
                og_dict[og_name].update(gene_list)
    
    seq_list = list(SeqIO.parse(input_filename, "fasta"))

    for seq in seq_list:
        seq_name = seq.id
        new_name = seq_name

        for og_name, genes in og_dict.items():
            if any(gene in seq_name for gene in genes):
                new_name = f"{og_name}_{seq_name}"
                break

        seq.id = new_name
        seq.description = ""

    SeqIO.write(seq_list, output_filename, "fasta")
    print(f"updated sequences written to {output_filename}.")

def extract_protein_hmmsearch(input_file: str, output_file: str, hmm_file) -> None:
    """
    Extract protein sequences from a FASTA file based on a HMM search result.

    Args
    ----
    input_file : str
        Path to the input FASTA file containing protein sequences.
    output_file : str
        Path to the output FASTA file where the extracted protein sequences will be saved.
    hmm_file : str
        Path to the HMM search result file.
        This file should contain protein names in the second column, which will be used to filter the sequences.
    """
    protein_set = set()

    with open(hmm_file, mode="r") as hmm_handle:
        for line in hmm_handle:
            if line.startswith("#") or not line.strip():
                continue
            li = line.strip().split()
            if len(li) < 14:
                continue
            protein_name = li[0]
            protein_set.add(protein_name)

    # === Sort protein names and write to output file ===
    protein_set = sorted(protein_set)
    protein_list = []
    for record in SeqIO.parse(input_file, "fasta"):
        if record.id in protein_set:
            protein_list.append(record)

    with open(output_file, "w") as output_handle:
        SeqIO.write(protein_list, output_handle, "fasta")



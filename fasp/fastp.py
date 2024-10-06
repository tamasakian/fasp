#!/usr/bin/env python3

"""Library for processing protein FASTA files.

Functions
---------
exclude_isoforms_by_length: Exclude isoforms based on length.

"""

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
            print(f"Gene: {gene}, CDS List: {cds_list}, Protein Lengths: {protein_lengths}")
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

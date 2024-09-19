#!/usr/bin/env python3

import argparse
import sys
from fasp import fasta
from fasp import fastp
from fasp import fastn

def main():
    args = parse_args()
    function = functions[args.function]
    try:
        function(*args.args)
    except TypeError as e:
        print(e)
        sys.exit(1)

def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("function", choices=functions)
    parser.add_argument("args", type=str, nargs="*")
    args = parser.parse_args()
    return args

functions = {
    "rename_header": fasta.rename_header, 
    "prefix_to_sequence_ids": fasta.prefix_to_sequence_ids, 
    "split_multi_to_single": fasta.split_multi_to_single, 
    "merge_msa_by_ids": fasta.merge_msa_by_ids, 
    "slice_records_by_exact_ids": fasta.slice_records_by_exact_ids, 
    "slice_records_by_partial_ids": fasta.slice_records_by_partial_ids, 
    "summarize_lengths": fasta.summarize_lengths, 
    "rename_headers_feature": fastn.rename_headers_feature, 
    "slice_records_by_seqids": fastn.slice_records_by_seqids, 
    "slice_sequence_by_flanking_region": fastn.slice_sequence_by_flanking_region,
    "slice_sequence_by_upstream_region": fastn.slice_sequence_by_upstream_region,
    "slice_sequence_by_downstream_region": fastn.slice_sequence_by_downstream_region,
    "generate_introns": fastn.generate_introns
}

if __name__ == "__main__":
    main()

#!python3

"""
Author: "Akhil Kumar"
Email: "akhil.kumar@alumni.iitd.ac.in"

This script analyzes the fasta file and saves the desired statistics 
into a csv file with a same name as the input fasta file just with a 
.csv extension instead of a .fasta or a .fa extension
"""

import math
import os
import pandas as pd
from pathlib import Path

from Bio import SeqIO
from Bio.Seq import Seq
from Bio import SeqUtils
from tqdm import tqdm


def fasta_file_analyzer(fasta_file_path, csv_file_path):
    """Input: Path where fasta files are located, Path where csv files should be saved;
    Output: Creates and saves dataframes to the csv_file_path location
    """
    seq_objects = SeqIO.parse(fasta_file_path, 'fasta')

    datalist = []

    seq_records = (seq_obj for seq_obj in seq_objects)

    for record in tqdm(seq_records):

        temp_seq = record.seq
        seq_len = len(record)
        if not seq_len:
            continue
        count_A = temp_seq.count('A')
        count_C = temp_seq.count('C')
        count_G = temp_seq.count('G')
        count_T = temp_seq.count('T')

        ambiguous_base_count = seq_len - (count_A + count_C + count_G + count_T)
        seq_len_adjusted = seq_len - ambiguous_base_count

        count_CG = temp_seq.count('CG')
        count_GC = temp_seq.count('GC')

        percent_A = 100 * (count_A / seq_len_adjusted)
        percent_C = 100 * (count_C / seq_len_adjusted)
        percent_G = 100 * (count_G / seq_len_adjusted)
        percent_T = 100 * (count_T / seq_len_adjusted)

        percent_CG = 100 * (count_CG / seq_len_adjusted)
        percent_GC = 100 * (count_GC / seq_len_adjusted)
        GC_content = percent_G + percent_C

        ObyE_CG = (count_CG * seq_len_adjusted) / (count_C * count_G)
        ObyE_GC = (count_GC * seq_len_adjusted) / (count_G * count_C)

        datalist.append((record.name, seq_len, seq_len_adjusted, 
                         ambiguous_base_count, count_CG, 
                         ObyE_CG, ObyE_GC, 
                         percent_A, percent_C, percent_G, percent_T, 
                         percent_CG, GC_content))

    df = pd.DataFrame(datalist, columns=['info', 'seq_len', 'seq_len_adjusted', 
                                         'ambiguous_base_count', 'count_CG', 
                                         'ObyE_CG', 'ObyE_GC', 
                                         'percent_A', 'percent_C', 'percent_G', 'percent_T', 
                                         'percent_CG', 'GC_content'])

    df.to_csv(csv_file_path, index=False)


if __name__ == '__main__':

    fasta_files_dir = Path('../fasta_files/')
    for file in fasta_files_dir.iterdir():
        if (file.suffix in ('.fasta', '.fa') and
            not file.stem.startswith('.')):
            # For users with Python v<3.6, convert the input and output file path to 
            # strings before passing them to fasta_file_analyzer()
            input_file_path = file.resolve()
            output_file_path = fasta_files_dir.parent.joinpath(f'csv_files/{file.stem}.csv')
            print(f'Preliminary analysis of {file.name} begins')
            fasta_file_analyzer(input_file_path, output_file_path)

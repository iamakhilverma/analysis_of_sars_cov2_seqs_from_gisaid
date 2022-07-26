#!python3

"""
Author: "Akhil Kumar", "Nishank Goyal"
Email: "akhil.kumar@alumni.iitd.ac.in", "nishankaligarh@gmail.com"

This python script maps ZAP-binding motifs (i.e., C(n_{m})G(n)CG,
where m = 4/5/6/7/8) and CpG sites by analyzing MSA for the
1,410,423 full-length sequences. It matches ZAP-binding motifs and
CpG sites in each sequence and extracts their locations in the
genome(mapping was done with the SARS-CoV-2 WIV04 reference sequence).
It also gives us the frequency of the ZAP-binding motifs and CpG sites
at the locations where the motifs were found within the MSA file.

MSA- multiple sequence alignment
"""

import re
import itertools

# import shutil
from collections import Counter

import pandas as pd
from tqdm import tqdm


# # Uncomment only if the MSA file isn't unpacked
# # Unpacking the MSA file downloaded from GISAID
# shutil.unpack_archive(
#   '../fasta_files/alignment/(full)msa_0902.tar.xz',
#   '../fasta_files/'
#   )


def find_motif(motif, genome):
    """Locates all occurences of the input motif in the
    input genome and returns a list of their start positions
    within that genome
    """
    return [m.start() for m in re.finditer(rf"(?={motif})", genome)]


def process(info, genome, motif, datafile):
    """Uses find_motif() and appends a list of three entries to the datafile:
    1. the sequence's (i.e, input genome's)
    header,
    2. number of input motifs found within that input genome, and
    3. the matched motifs themselves

    Input:
    - datafile is where results will be stored,
    - genome is the sequence where motifs will be searched,
    - info is header of the sequence, and
    - motif is the pattern which should be matched in the genome

    Acceptable datatypes:
    - datafile should be a list/tuple,
    - genome should be a string,
    - info should be a string, and
    - motif should be a raw string.
    """
    genome = genome.replace("\n", "")
    matches = find_motif(motif, genome)
    datafile.append([info, len(matches), matches])


def convert_to_int(loc):
    """Converts the type of the locations to int"""
    return list(map(int, loc))


def motif_analyzer(motif, outfile_name):
    """Finds the locations of the motif within the full MSA file. It
    then maps it to the reference genome (mapping was done with the
    SARS-CoV-2 WIV04 reference sequence). It also gives us the frequency
    and percentage of the ZAP-binding motifs and CpG sites at the locations
    where the motifs were found within the MSA file.

    Input:
    - motif is the pattern which should be matched in the genome. It should
        be a raw string
    - outfile_name is the name of the output file. It should have an
        .csv extension

    Output:
    - Saves a .csv file with the name outfile_name in the csv_files directory
    """
    print(f'Analysis for "{motif}" motif begins...')
    data = []

    with open("../fasta_files/msa_0902/msa_0902.fasta") as infile:
        # we can use itertools.zip_longest(infile, infile) as well
        for line1, line2 in tqdm(itertools.zip_longest(*[infile] * 2)):
            if line1 and line2 and line1.startswith(">"):
                process(info=line1[1:], genome=line2, motif=motif, datafile=data)

    data = pd.DataFrame(data, columns=["info", "number", "identified_motif_locations"])
    data["info"] = data["info"].replace(r"\n", "", regex=True)
    data["accession_ids"] = data["info"].str.split("|").str[1]

    df_accession_ids = pd.read_csv("../other_data_files/accession_ids.txt", header=None)
    df_accession_ids = df_accession_ids.rename(columns={0: "accession_ids"})
    # Keeping only the sequences with desired accession ids
    data = data[data.accession_ids.isin(df_accession_ids.accession_ids)]

    print(
        f"Total entries analyzed: {len(df_accession_ids)}\n"
        "No. of entries from desired accession ids not found"
        f" in our df with 2.8 mil seqs: {len(df_accession_ids) - len(data)}"
    )

    count = Counter()
    data.apply(
        lambda x: count.update(convert_to_int(x["identified_motif_locations"])), axis=1
    )
    count = pd.DataFrame.from_dict(count, orient="index", columns=["frequency"])
    count["percentage"] = (count["frequency"] * 100) / len(data)

    mapper = pd.read_excel("../other_data_files/mapping.xlsx")
    mapper = mapper.set_index(["Unnamed: 0"])
    count["mapping"] = mapper.loc[count.index]["PreAlignedNumber"]

    count.to_csv(f"../csv_files/{outfile_name}", index_label="python_based_indices")


def zap_n_motif_constructor(m):
    """Constructs the pattern for ZAP-binding motif given n
    Input: m is the m in C(n_{m})G(n)CG, where m = 4/5/6/7/8
    Returns the raw string for the desired ZAP-binding motif
    """
    zap_motif = r"[cC]-*"
    for _ in range(m):
        zap_motif = zap_motif + r"[aAtTcCgG]-*"
    return zap_motif + r"[gG]-*[aAtTcCgG]-*[cC]-*[gG]"


if __name__ == "__main__":
    motif_analyzer(motif=r"[cC]-*[gG]", outfile_name="counter_file_cpg.csv")
    motif_analyzer(
        motif=zap_n_motif_constructor(4), outfile_name="counter_file_zap_4_motif.csv"
    )
    motif_analyzer(
        motif=zap_n_motif_constructor(5), outfile_name="counter_file_zap_5_motif.csv"
    )
    motif_analyzer(
        motif=zap_n_motif_constructor(6), outfile_name="counter_file_zap_6_motif.csv"
    )
    motif_analyzer(
        motif=zap_n_motif_constructor(7), outfile_name="counter_file_zap_7_motif.csv"
    )
    motif_analyzer(
        motif=zap_n_motif_constructor(8), outfile_name="counter_file_zap_8_motif.csv"
    )

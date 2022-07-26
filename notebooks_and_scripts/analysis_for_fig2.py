#!python3

"""
Author: "Nishank Goyal"
Email: "nishankaligarh@gmail.com"

This script contains the analysis for Fig. 2.(a) and 2.(b) of our 
publication.
"""

import pandas as pd
import numpy as np
import re
import matplotlib.pyplot as plt
import matplotlib.axes as ax
import math

data = pd.read_csv("../csv_files/filtered_sars_cov2.csv")

data = data.groupby(["dates"]).mean()  # this step is done on 14,10,423 sequences
data = data.rename(
    columns={
        "seq_len": "with",
        "seq_len_adjusted": "without",
        "ambiguous_base_count": "N",
        "count_CG": "CG",
        "ObyE_CG": "o/eCG",
        "ObyE_GC": "o/eGC",
        "percent_CG": "%CG",
        "percent_A": "perA",
        "percent_C": "perC",
        "percent_G": "perG",
        "percent_T": "perT",
    }
)
data = data.drop(["GC_content"], axis=1)

cases = pd.read_excel(
    "../other_data_files/world_data.xlsx"
)  # OWiD dataset maintained using Johns Hopkins University CSSE COVID-19 Data

# taking mean of 52 entries till 21st January (15 days), assigning the mean value to 21st January and value of cases 0, because OWiD data starts from 22nd January
data = data.reset_index()
data.loc[511] = data[:15].mean()
data = data.drop(data[:15].index)
data["dates"][511] = "2020-01-21"

data = (
    data.set_index("dates").sort_index().reset_index()
)  # setting date as iindex, sorting the index and resetting the index to make database sorted according to date

data["cases"] = np.append(
    [0], cases["Total cases"].values
)  # adding 0 on top of cases value from OWiD dataset for 21st January 2020
data["month"] = [
    int(data["dates"][i].split("-")[1])
    + ((int(data["dates"][i].split("-")[0]) - 2020) * 12)
    for i in range(len(data))
]  # converting dates into months

data = data[["dates", "without", "CG", "%CG", "o/eCG", "o/eGC", "cases"]]

final = data.set_index(["cases"])

# function to get change percentage of each data point w.r.t series maximum and minimum
def changePercentage(change_column):
    change = []
    for i in range(len(final.index)):
        change.append(
            (
                (max(final[change_column]) - final[change_column].iloc[i])
                * 100
                / (max(final[change_column]) - min(final[change_column]))
            )
        )
    return change


final["perChange"] = changePercentage(
    "%CG"
)  # getting change percentage for each data point in CG percentage data
final["numChange"] = changePercentage(
    "CG"
)  # getting change percentage for each data point in CG number data

for i in ["CG", "%CG"]:

    name = i

    left_labels = {"CG": "Number of CpGs", "%CG": "Percentage of CpGs"}
    right_labels = {
        "CG": "Percentage reduction in no. of CpGs",
        "%CG": "Percentage reduction in CpG percentage",
    }
    save_labels = {"CG": "Number_CG", "%CG": "Percentage_CG"}

    fig = plt.figure()
    ax = fig.add_axes([0.1, 0.1, 0.8, 0.8])

    ax.plot(
        final.index,
        final[name],
        linewidth=0.75,
        color="#0173b2",
        label="Daily average values",
    )
    x = final.index
    y = final[name].rolling(10).mean()
    ci = 2.262 * final[name].rolling(10).std() / np.sqrt(10)
    ax.plot(x, y, linewidth=1.25, color="#de8f05", label="10-days Moving Average")

    ax.yaxis.grid(True, which="major")
    ax.xaxis.grid(True, which="major")
    if name in ["CG", "%CG"]:
        ax1 = ax.twinx()
        ax1.set_ylim(105, -5)
        ax1.set_yticks([100, 80, 60, 40, 20, 0])
        ax1.set_yticklabels(["100%", "80%", "60%", "40%", "20%", "0%"], fontsize=11)
        ax1.set_ylabel(right_labels[name], fontsize=12)
        ax.fill_between(
            x / 3,
            [min(final[name])] * len(x),
            [max(final[name])] * len(x),
            color="#d55e00",
            alpha=0.1,
        )
        ax.fill_between(
            (max(x) / 3) + ((2 * x) / 3),
            [min(final[name])] * len(x),
            [max(final[name])] * len(x),
            color="#029e73",
            alpha=0.1,
        )
    ax.fill_between(
        x, (y - ci), (y + ci), color="g", alpha=0.3, label="95% Confidence Interval"
    )

    ax.set_xticks(
        [
            0,
            20000000,
            40000000,
            60000000,
            80000000,
            100000000,
            120000000,
            140000000,
            160000000,
            180000000,
        ]
    )
    ax.set_xticklabels(
        ["0M", "20M", "40M", "60M", "80M", "100M", "120M", "140M", "160M", "180M"],
        fontsize=11,
    )
    ax.set_xlabel("Number of Infections", fontsize=12)
    ax.set_ylabel(left_labels[name], fontsize=12)
    ax.tick_params(axis="both", which="major", labelsize=11)
    fig.legend(loc="center left", bbox_to_anchor=(0.515, 0.845), prop={"size": 11})
    plt.show()

#!/usr/bin/env python
import random
import argparse
from os.path import basename, splitext

import HTSeq
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from scipy.stats import pearsonr

parser = argparse.ArgumentParser()
parser.add_argument("-c", "--chromsize", metavar="FILE", help="Index file of "
                    "reference genome")
parser.add_argument("-x", "--wxs", metavar="FILE", help="Output from GATK4 "
                    "CNV using whole exome sequencing (WXS) data")
parser.add_argument("-g", "--wgs", metavar="FILE", help="Output from BIC-seq2 "
                    "using whole genome sequencing (WGS) data")
parser.add_argument("-r", "--round", metavar="INT", type=int, default=10,
                    help="Sample points from input data for a given rounds "
                    "(default: %(default)s)")
parser.add_argument("-s", "--sample", metavar="INT", type=int, default=100000,
                    help="Randomly sample a given number of points from input "
                    "data (default: %(default)s)")
args = parser.parse_args()

included_chroms = list(map(str, range(1, 23))) + ["X", "Y"]
chrom_lens = {}
with open(args.chromsize) as infile:
    for line in infile:
        line = line.strip().split("\t")
        if line[0] in included_chroms:
            chrom = "chr" + line[0]
            length = int(line[1])
            chrom_lens[chrom] = length

wxs_array = HTSeq.GenomicArray(chrom_lens, stranded=False, typecode="d")
with open(args.wxs) as infile:
    for i, line in enumerate(infile):
        if i > 0:
            line = line.strip().split("\t")
            chrom = "chr" + line[1] 
            start = int(line[2]) - 1
            end = int(line[3])
            cnv = float(line[5])
            iv = HTSeq.GenomicInterval(chrom, start, end, ".")
            wxs_array[iv] += cnv

wgs_array = HTSeq.GenomicArray(chrom_lens, stranded=False, typecode="d")
with open(args.wgs) as infile:
    for i, line in enumerate(infile):
        if i > 0:
            line = line.strip().split("\t")
            chrom = line[0]
            start = int(line[1]) - 1
            end = int(line[2])
            cnv = 2**float(line[8])
            iv = HTSeq.GenomicInterval(chrom, start, end, ".")
            wgs_array[iv] += cnv

pool = set()
pearsonr_data = {"pearsonr": [], "pval": []}
prefix1 = splitext(basename(args.wxs))[0]
prefix2 = splitext(basename(args.wgs))[0]
for i in range(args.round):
    wxs_data = []
    wgs_data = []
    while len(wxs_data) < args.sample:
        chrom = random.choice(list(chrom_lens))
        pos = random.randrange(chrom_lens[chrom])
        if (chrom, pos) not in pool:
            iv = HTSeq.GenomicInterval(chrom, pos, pos+1, ".")
            wxs_cnv = list(wxs_array[iv])[0]
            wgs_cnv = list(wgs_array[iv])[0]
            if wxs_cnv > 0 and wgs_cnv > 0:
                pool.add((chrom, pos))
                wxs_data.append(wxs_cnv)
                wgs_data.append(wgs_cnv)

    sns.set(style="white", color_codes=True)
    minval = min(wxs_data + wgs_data)
    maxval = max(wxs_data + wgs_data)

    r, pval = pearsonr(wxs_data, wgs_data)
    pearsonr_data["pearsonr"].append(r)
    pearsonr_data["pval"].append(pval)
    g = sns.jointplot("WXS", "WGS", data=pd.DataFrame({"WXS": wxs_data,
                      "WGS": wgs_data}), xlim=(minval, maxval),
                      ylim=(minval, maxval), space=0, kind="reg",
                      joint_kws={"scatter_kws": {"alpha": 0.1}})
    plt.savefig("{0}_{1}.{2}.{3}.png".format(prefix1, prefix2, args.sample,
                str(i+1).zfill(len(str(args.round)))), dpi=300)
    plt.close()

pearsonr_data = pd.DataFrame(pearsonr_data)
pearsonr_data = pearsonr_data[["pearsonr", "pval"]]
pearsonr_data.to_csv("{0}_{1}.{2}.pearsonr.txt".format(prefix1, prefix2,
                     args.sample), index=False, sep="\t", float_format="%.6f")

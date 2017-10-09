#!/usr/bin/env python
import argparse
from os.path import basename, splitext

import HTSeq
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

parser = argparse.ArgumentParser()
parser.add_argument("-c", "--chromsize", help="Index file of reference genome")
parser.add_argument("-x", "--wxs", help="Output from GATK4 CNV using whole exome sequencing (WXS) data")
parser.add_argument("-g", "--wgs", help="Output from BIC-Seq2 using whole genome sequencing (WGS) data")
parser.add_argument("-w", "--winsize", type=int, default=10000, help="Moving window size. (default: %(default)s)")
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
            sample, chrom, start, end, numprobes, segmean, segcall = line.strip().split("\t")
            chrom = "chr" + chrom 
            start = int(start) - 1
            end = int(end)
            cnv = float(segmean)
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

plt.switch_backend('Agg')
sns.set(style="white", color_codes=True)
fig = plt.figure(figsize=(30, 12))
ax1 = fig.add_subplot(211)
ax2 = fig.add_subplot(212)
prelen = 0
x = []
gy = []
xy = []
for chrom in included_chroms:
    chrom = "chr" + chrom
    for start in range(0, chrom_lens[chrom] - args.winsize, args.winsize):
        end = start + args.winsize
        iv = HTSeq.GenomicInterval(chrom, start, end)
        x.append(0.5*(prelen + start + prelen + end))
        gy.append(np.mean(list(wgs_array[iv])))
        xy.append(np.mean(list(wxs_array[iv])))
    prelen += chrom_lens[chrom]

ax1.plot(x, xy, "o", ms=6, alpha=1)
prelen = 0
xticks = []
xlabels = []
for chrom in included_chroms:
    chrom = "chr" + chrom
    xlabels.append(chrom)
    if prelen > 0:
        ax1.axvline(prelen, color="grey", ls="--", alpha=0.6)
    xticks.append(0.5*(prelen + prelen + chrom_lens[chrom]))
    prelen += chrom_lens[chrom]
ax1.set_xlim(0, prelen)
ax1.set_ylim(0, 2)
ax1.set_xticks(xticks)
ax1.set_xticklabels(xlabels)
ax1.set_ylabel("WXS", size="xx-large")

ax2.plot(x, gy, "o", ms=6, alpha=1)
prelen = 0
xticks = []
xlabels = []
for chrom in included_chroms:
    chrom = "chr" + chrom
    xlabels.append(chrom)
    if prelen > 0:
        ax2.axvline(prelen, color="grey", ls="--", alpha=0.6)
    xticks.append(0.5*(prelen + prelen + chrom_lens[chrom]))
    prelen += chrom_lens[chrom]
ax2.set_xlim(0, prelen)
ax2.set_ylim(0, 2)
ax2.set_xticks(xticks)
ax2.set_xticklabels(xlabels)
ax2.set_ylabel("WGS", size="xx-large")

plt.tight_layout()
plt.savefig("{0}_{1}.genomecnv.png".format(splitext(basename(args.wxs))[0], splitext(basename(args.wgs))[0]), dpi=300)
plt.close()

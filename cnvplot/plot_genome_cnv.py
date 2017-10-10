#!/usr/bin/env python
import argparse
from os.path import basename, splitext

import HTSeq
import seaborn as sns
import matplotlib.pyplot as plt

parser = argparse.ArgumentParser()
parser.add_argument("-c", "--chromsize", metavar="FILE", help="Index file of "
                    "reference genome")
parser.add_argument("-x", "--wxs", metavar="FILE", help="Output from GATK4 "
                    "CNV using whole exome sequencing (WXS) data")
parser.add_argument("-g", "--wgs", metavar="FILE", help="Output from BIC-seq2 "
                    "using whole genome sequencing (WGS) data")
parser.add_argument("-w", "--winsize", metavar="INT", type=int, default=10000,
                    help="Moving window size. (default: %(default)s)")
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

plt.switch_backend('Agg')
sns.set(style="white", color_codes=True)
fig = plt.figure(figsize=(30, 12))
ax1 = fig.add_subplot(211)
ax2 = fig.add_subplot(212)
prelen = 0
x = []
wgs_y = []
wxs_y = []
for chrom in included_chroms:
    chrom = "chr" + chrom
    for start in range(0, chrom_lens[chrom] - args.winsize, args.winsize):
        end = start + args.winsize
        iv = HTSeq.GenomicInterval(chrom, start, end)
        x.append(prelen + 0.5*(start + end))
        wxs_window = list(wxs_array[iv])
        wxs_y.append(sum(wxs_window)/len(wxs_window))
        wgs_window = list(wgs_array[iv])
        wgs_y.append(sum(wgs_window)/len(wgs_window))
    prelen += chrom_lens[chrom]

ax1.plot(x, wxs_y, "o", ms=6, alpha=0.6)
prelen = 0
xticks = []
xticklabels = []
for chrom in included_chroms:
    xticklabels.append(chrom)
    if prelen > 0:
        ax1.axvline(prelen, color="gray", ls="--", alpha=0.6)
    length = chrom_lens["chr" + chrom]
    xticks.append(prelen + 0.5*(length))
    prelen += length
ax1.set_xlim(0, prelen)
ax1.set_ylim(0, 2)
ax1.set_xticks(xticks)
ax1.set_xticklabels(xticklabels)
ax1.tick_params(labelsize="xx-large")
ax1.set_xlabel("Chromosome", size="xx-large")
ax1.set_ylabel("WXS", size="xx-large")

ax2.plot(x, wgs_y, "o", ms=6, alpha=0.6)
prelen = 0
xticks = []
xticklabels = []
for chrom in included_chroms:
    xticklabels.append(chrom)
    if prelen > 0:
        ax2.axvline(prelen, color="gray", ls="--", alpha=0.6)
    length = chrom_lens["chr" + chrom]
    xticks.append(prelen + 0.5*(length))
    prelen += length
ax2.set_xlim(0, prelen)
ax2.set_ylim(0, 2)
ax2.set_xticks(xticks)
ax2.set_xticklabels(xticklabels)
ax2.tick_params(labelsize="xx-large")
ax2.set_xlabel("Chromosome", size="xx-large")
ax2.set_ylabel("WGS", size="xx-large")

plt.tight_layout()
prefix1 = splitext(basename(args.wxs))[0]
prefix2 = splitext(basename(args.wgs))[0]
plt.savefig("{0}_{1}.genomecnv.png".format(prefix1, prefix2, dpi=300))
plt.close()

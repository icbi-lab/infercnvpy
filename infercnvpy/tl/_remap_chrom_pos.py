"""Remap gene buckets on actual chromosome positions"""

from .._util import tqdm
import re
import pandas as pd
import numpy as np


def rescale_to_chromosome(X, offset_col, chr_col, binwidth=50000):
    """Rescale gene buckets to chromosome lengths

    The infercnv algorithm (and also SCEVAN) yield a result matrix where the
    number of entries per chromosome is not proportional to the length
    of the chromosome, but to the number of genes per chromosome.

    This can look misleading in the chromosome heatmap.

    This function remaps buckets of arbitrary length to buckets proportional
    to the chromosome length by computing the mean of all input buckets that
    fall within a fixed interval (e.g. 50kb) of genomic positions.

    TODO: this is untested, slow and needs to be improved

    Parameters
    ----------
    X
        input CNV matrix
    offset_col
        Numpy array giving the genomic position of each bucket in X
    chr_col
        numpy array giving the chromosome for each bucket in X
    binwidth
        size of the output bins

    Returns
    -------
    rescaled matrix and chr_pos dictionary
    """
    chrom_sizes = {
        chr: length
        for _, chr, length in pd.read_csv(
            "https://github.com/igvteam/igv/raw/master/genomes/sizes/hg38.chrom.sizes",
            sep="\t",
            names=["chr", "length"],
        ).itertuples()
        if re.match("chr\d+$", chr) is not None
    }

    bins = {}
    for c, length in tqdm(chrom_sizes.items()):
        tmp_bins = []
        tmp_x = X[:, chr_col == c]
        tmp_offset = offset_col[chr_col == c].tolist()
        xi = 0
        for i in range(0, length, binwidth):
            tmp_ind = []
            while len(tmp_offset) and tmp_offset[0] < i:
                tmp_offset.pop(0)
                tmp_ind.append(xi)
                xi += 1
            if len(tmp_ind):
                tmp_bins.append(np.mean(tmp_x[:, tmp_ind], axis=1)[:, np.newaxis])
            else:
                #                 tmp_bins.append(np.zeros(X.shape[0])[:, np.newaxis])
                # repeat last state until there's a new one
                try:
                    tmp_bins.append(tmp_bins[-1])
                except IndexError:
                    tmp_bins.append(np.zeros(X.shape[0])[:, np.newaxis])

        bins[c] = np.hstack(tmp_bins)

    chr_pos = {
        c: p
        for c, p in zip(
            bins.keys(), np.cumsum([0] + [x.shape[1] for x in bins.values()])
        )
    }
    return np.hstack(bins.values()), chr_pos

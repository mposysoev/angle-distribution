#!/usr/bin/python3

import numpy as np
from matplotlib import pyplot as plt
import sys


if __name__ == "__main__":
    if len(sys.argv) != 4:
        print("""USAGE:
    python3 gistagramator <type> <bins> <file_name> 
    <type>: hist or scatter
    <bins>: int number
    <file_name>: file name""")
        sys.exit(0)
        
    type_plot = sys.argv[1]
    bins_numb = int(sys.argv[2])
    file_name = sys.argv[3]

    data = np.loadtxt(file_name)

    if type_plot == 'hist':
        plt.hist(data, bins=bins_numb)
        plt.savefig(f"{file_name}-histogram.png")
    if type_plot == 'scatter':
        hist, bin_edges = np.histogram(data, bins=5)
        shift = (bin_edges[1] - bin_edges[0]) / 2.0
        bins = bin_edges - shift
        plt.scatter(bins[1:-1], hist[1:])
        plt.savefig(f"{file_name}-scatter-histogram.png")

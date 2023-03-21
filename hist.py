import numpy as np
from matplotlib import pyplot as plt
import sys


if __name__ == "__main__":
    file_name = sys.argv[1]
    data = np.loadtxt(file_name)
    hist, bin_edges = np.histogram(data, bins=50)
    shift = (bin_edges[1] - bin_edges[0]) / 2.0
    bins = bin_edges - shift
    plt.scatter(bins[1:-1], hist[1:] / np.sin(bins[1:-1]))
    plt.savefig(f"{file_name}-out.png")

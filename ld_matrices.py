import pandas as pd
import argparse
import gzip
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap


def main(args):
    ld1 = pd.read_table(args.ld1, header=None)
    ld2 = pd.read_table(args.ld2, header=None)

    ld1 = ld1.values
    ld2 = ld2.values

    i_upper = np.triu_indices(ld1.shape[0])
    ld2[i_upper] = ld1.T[i_upper]

    # Set up the matplotlib figure
    f, ax = plt.subplots(figsize=(6, 5))

    # Generate a custom diverging colormap
    cmap = ListedColormap(sns.color_palette('Blues', 256))

    # Draw the heatmap with the mask and correct aspect ratio
    sns.set(font_scale=2.5)
    sns.heatmap(ld2, cmap=cmap)
    ax.set(title=args.ld1_pop, ylabel=args.ld2_pop)
    for item in ([ax.title, ax.yaxis.label]):
        item.set_fontsize(24)
    plt.tick_params(axis='x', bottom=False, labelbottom=False)
    plt.tick_params(axis='y', left=False, labelleft=False)
    fig = plt.gcf()
    fig.savefig(args.out, dpi=300, transparent=True)


if __name__ == '__main__':
    print('Starting run')
    parser = argparse.ArgumentParser()
    parser.add_argument('--ld1')
    parser.add_argument('--ld1_pop')
    parser.add_argument('--ld2')
    parser.add_argument('--ld2_pop')
    parser.add_argument('--out')
    args = parser.parse_args()
    main(args)

import matplotlib
import matplotlib.pyplot as plt
import seaborn as sns
matplotlib.use("Agg")

from sklearn.neighbors import KernelDensity


def init_plot(figsize=(12, 4)):
    fig = plt.figure(figsize=figsize)

    sns.set_context("paper")
    sns.set_style(
        "ticks",
        {
            "axes.grid": False,
            "axes.linewidth": "0.75",
            "grid.color": "0.75",
            "grid.linestyle": u":",
            "legend.frameon": True,
        },
    )
    plt.rc("text", usetex=True)
    plt.rc("font", family="Serif")

    return fig


def save_fig(fig, file_name, print_output=True):
    fig.tight_layout(pad=0.5)

    if print_output:
        print("\t Writing figure to {}...".format(file_name))

    fig.savefig(
        file_name, bbox_inches="tight", dpi=180, facecolor="w", transparent=True
    )


"""
Kernel Density Estimation plot based on scikit-learn

Parameters
----------
x, y
"""


def kdeplot():
    kde = KernelDensity(
        bandwidth=0.1, metric="haversine", kernel="gaussian", algorithm="kd_tree"
    )

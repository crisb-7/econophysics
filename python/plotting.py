import seaborn as sns
import matplotlib.pyplot as plt

# General pretty plot
def line_plot(x, y):
    pass

# General pretty histogram
def simple_histogram(x, bins):
    pass

# Specific plots - remove soon
def plot_entropy(entropy, figsize):
    fig, ax = plt.subplots(figsize=figsize)
    ax.plot(entropy)
    plt.show()

def plot_money_hist(hist_data, bin_edges, figsize):
    fig, ax = plt.subplots(figsize=figsize)
    sns.histplot(x=hist_data, bins=bin_edges, ax=ax)
    plt.show()
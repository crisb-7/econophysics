import random
import numpy as np
from plotting import plot_entropy, plot_money_hist

# TODO: rethink design - implement well-being & tax 
# TODO: organize constant definitions/make; script?
# TODO: modularize to use delta/uniform initial distribution
# TODO: make constants much more descriptive
# TODO: revise implementation details for performance (speed + memory)
# TODO: add readme description of constants
# TODO: develop tests + test code
# TODO: develop a visualization module 

def main():
    N = 1000    # Population (N agents)
    M = 1e5     # Total sum of money in the system
    C = 50      # Number of income classes
    T = 1000    # Total time of the simulation

    # Delta amount of money
    nu = 0.05
    dm = (M/N)*nu

    # Max amount of money an agent can have
    MMAX = (M/(np.sum(np.arange(C))*(N/C))*C)*2
    CLASS = np.linspace(0, MMAX, C)
    bin_centers = 0.5*(CLASS[1:] + CLASS[0:-1])

    money_distrib = np.ones((N,))*(M/N)
    Entropy = np.zeros((T,))

    # TODO: verify correctness - code runs without errors
        # 06 Jan - Weird entropy drops T > 18000
    verbose = True
    if verbose:
        print("Number of Agents:", N)
        print("Time iterations (unit):", T)
        print("Delta Money:", f"${dm}")

    for t in range(T):

        # Choose two agents
        agent_i, agent_j = choose_agents(N)

        # Perform money exchange between them
        money_distrib = exchange_money(agent_i, agent_j, money_distrib, dm)
        
        # Histogram for money distribution
        hist, bin_edges = np.histogram(money_distrib, bins=CLASS)
        
        # Entropy function
        entropy = N*np.log(N) - np.sum(hist[hist>0]*np.log(hist[hist>0]))
        Entropy[t] = entropy

    # note: these pop up one at a time
    plot_entropy(Entropy, (10,4))
    plot_money_hist(hist, CLASS, (10,4))

def choose_agents(n_population: int):
    i = random.randint(0, n_population-1)
    j = random.randint(0, n_population-1)
    while i == j:
        j = random.randint(0, n_population-1)
    return i, j

def exchange_money(agent_i: int, agent_j: int, distribution: np.array, money: int) -> np.array:
    """"Perform money exchange between agents.

    Returns
    -------
        distribution: money distribution updated with the money exchange
    """
    # Determine which agent will win/lose money
    s = random.choice([-1, 1])

    if (distribution[agent_i] + s*money < 0) or (distribution[agent_j] - s*money < 0):
        return distribution
    
    distribution[agent_i] = distribution[agent_i] + s*money
    distribution[agent_j] = distribution[agent_j] - s*money

    return distribution

if __name__ == "__main__":
    main()
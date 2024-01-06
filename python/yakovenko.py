import argparse
import random
import numpy as np
from plotting import plot_entropy, plot_money_hist

def main(args):
    N_AGENTS = args.population     # Population (N agents) 
    TOTAL_MONEY = args.total_money    # Total sum of money in the system
    N_CLASSES = args.nclasses       # Number of income classes
    TIME = args.time           # Unit time of the simulation
    EXCH_FRAC = args.exch_frac

    EXCH_MONEY = (TOTAL_MONEY/N_AGENTS)*EXCH_FRAC   # Amount of money per exchange

    verbose = args.verbose
    if verbose:
        print("Number of Agents:", N_AGENTS)
        print("Time iterations (unit):", TIME)
        print("Total Money in System:", f"${TOTAL_MONEY}")
        print("Money Exchange Fraction:", EXCH_FRAC)
        print("Money per Exchange:", f"${EXCH_MONEY}")
        print("Number of Income Classes:", N_CLASSES)

    # Max amount of money an agent can have
    # Break this down further
    max_agent_money = (TOTAL_MONEY/(np.sum(np.arange(N_CLASSES))*(N_AGENTS/N_CLASSES))*N_CLASSES)*2
    class_edges = np.linspace(0, max_agent_money, N_CLASSES)
    class_centers = 0.5*(class_edges[1:] + class_edges[0:-1])

    # Init distribution and entropy arrays
    money_distrib = np.ones((N_AGENTS,))*(TOTAL_MONEY/N_AGENTS)
    Entropy = np.zeros((TIME,))

    for t in range(TIME):

        # Choose two agents
        agent_i, agent_j = choose_agents(N_AGENTS)

        # Perform money exchange between them
        money_distrib = exchange_money(agent_i, agent_j, money_distrib, EXCH_MONEY)
        
        # Histogram for money distribution
        hist, bin_edges = np.histogram(money_distrib, bins=class_edges)
        
        # Entropy function
        entropy = N_AGENTS*np.log(N_AGENTS) - np.sum(hist[hist>0]*np.log(hist[hist>0]))
        
        # 06 Jan - Weird entropy drops T > 18000
        Entropy[t] = entropy

    # note: these pop up one at a time
    plot_entropy(Entropy, (10,4))
    plot_money_hist(hist, class_edges, (10,4))

def choose_agents(n_population: int):
    i = random.randint(0, n_population-1)
    j = random.randint(0, n_population-1)
    while i == j:
        j = random.randint(0, n_population-1)
    return i, j

def exchange_money(agent_i: int, agent_j: int, distribution: np.array, exch_money: int) -> np.array:
    """"Perform money exchange between agents.

    Returns
    -------
        distribution: money distribution updated with the money exchange
    """
    # Determine which agent will win/lose money
    s = random.choice([-1, 1])

    if (distribution[agent_i] + s*exch_money < 0) or (distribution[agent_j] - s*exch_money < 0):
        return distribution
    
    distribution[agent_i] = distribution[agent_i] + s*exch_money
    distribution[agent_j] = distribution[agent_j] - s*exch_money

    return distribution

if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument("--population", type=int, default=1000)
    parser.add_argument("--total_money", type=int, default=10000)
    parser.add_argument("--nclasses", type=int, default=50)
    parser.add_argument("--time", type=int, default=10000)
    parser.add_argument("--exch_frac", type=float, default=0.05)
    parser.add_argument("--verbose", type=bool, default=True)
    
    args = parser.parse_args()
    
    main(args)
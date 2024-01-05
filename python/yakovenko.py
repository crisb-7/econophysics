import random
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt

# TODO: organize constant definitions/make some user-defined
# TODO: make constants much more descriptive
# TODO: add readme description of constants
# TODO: revise implementation details for performance (speed + memory)
# TODO: develop tests + test code
# TODO: modularize to use delta/uniform initial distribution
# TODO: add logic to enable/disable tax model
# TODO: rethink design to implement well-being analysis
# TODO: Translate Entropy Plots
# TODO: Well-being Boltzmann-Gibbs distribution

N = 1000    # Population (N agents)
M = 1e5     # Total sum of money in the system
C = 50      # Number of income classes
T = 18000    # Total time of the simulation

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
    # 04 Jan - Entropy looks correct, simpler design
verbose = True
if verbose:
    print("Number of Agents:", N)
    print("Time iterations (unit):", T)
    print("Delta Money:", f"${dm}")

for t in range(T):

    # Choose two agents
    agent_i = random.randint(0, N-1)
    agent_j = random.randint(0, N-1)
    while agent_i == agent_j:
        agent_j = random.randint(0, N-1)

    s = random.choice([-1, 1])

    if (money_distrib[agent_i] + s*dm < 0) or (money_distrib[agent_j] - s*dm < 0):
        continue

    money_distrib[agent_i] = money_distrib[agent_i] + s*dm
    money_distrib[agent_j] = money_distrib[agent_j] - s*dm
    
    # Histogram for money distribution
    hist, bin_edges = np.histogram(money_distrib, bins=CLASS)
    
    # Entropy function
    Entropy[t] = N*np.log(N) - np.sum(hist*np.log(hist, where=(hist>0)))

# TODO: develop a visualization module 
# This is just a rough visualization
fig, ax = plt.subplots(1, 2, figsize=(10,4))
ax[0].plot(Entropy)

print(hist)
sns.histplot(x=hist, bins=bin_edges, ax=ax[1])
plt.show()
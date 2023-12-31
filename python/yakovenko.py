import numpy as np

# TODO: organize constant definitions/make some user-defined
# TODO: make constants much more descriptive
# TODO: add readme description of constants
# TODO: revise implementation details for performance (speed + memory)
# TODO: develop tests + test code
# TODO: develop a visualization module 
#       (this can be used for Wright's model too)
# TODO: modularize to use delta/uniform initial distribution
# TODO: add logic to enable/disable tax model
# TODO: rethink design to implement well-being analysis
# TODO: Translate Entropy Plots
# TODO: Well-being Boltzmann-Gibbs distribution

N = 1000    # Population (N agents)
M = 1e5     # Total sum of money in the system
C = 50      # Number of income classes

# Maximum amount of money an agent can have
MMAX = (M/(np.sum(np.arange(C))*(N/C))*C)*2

CLASS = np.linspace(0, MMAX, C)
T = 1000        
Ci = len(CLASS)
dm = M*4e-5
lambda_ = 0
tau = 30
a = 1/300

bin_centers = 0.5*(CLASS[1:] + CLASS[0:-1])
Mk_vec = np.linspace(0, MMAX, N)

iva = 0
delta_distrib = np.ones((N,))*(M/N)
Sdelta = np.zeros((T,))

OcaD = [None]*T
OcbD = [None]*T
B = [None]*T

# TODO: verify correctness - code runs without errors

for t in range(T):
    
    # histogram visualization used to go here
    if t % tau == 0:
        delta_distrib = delta_distrib + iva/N
        iva = 0

    for j in range(N):
        l = np.random.randint(N)
        s = np.sign(2*np.random.rand()-1)
        
        if (delta_distrib[j] + dm*s < 0) or delta_distrib[l] - dm*s < 0:
            continue
        elif s > 0:
            delta_distrib[j] = delta_distrib[j] + dm*(1-lambda_)
            delta_distrib[l] = delta_distrib[l] - dm
            iva = iva + dm*lambda_
        else:
            delta_distrib[j] = delta_distrib[j] - dm
            delta_distrib[l] = delta_distrib[l] + dm*(1-lambda_)
            iva = iva + dm*lambda_
    
    # Create histogram for money distribution
    hist, bin_edges = np.histogram(delta_distrib, bins=CLASS)
    print(id(Sdelta))
    Sdelta[t] = N*np.log(N) - np.sum(hist*np.log(hist, where=(hist>0)))

    # Well being functions
    OcaD[t] = sum(a*bin_centers*hist)
    OcbD[t] = sum(hist*(1-np.exp(-a*bin_centers)))
    
    VEm = 0
    for u in range(C-1):
        VEm = VEm + bin_centers[u]*hist[u]/N
    B[t] = VEm
    VEm = 0
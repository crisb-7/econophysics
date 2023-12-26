import numpy as np


# TODO: organize constant definitions/make some user-defined
# TODO: make constants much more descriptive
# TODO: add readme description of constants

N = 1000
M = 1e5
C = 50
MMAX = (M/(np.sum(np.arange(1, C))*(N/C))*C)*2

CLASS = np.linspace(0, MMAX, C+1)
T = 12000
Ci = len(CLASS)
dm = M*4e-5
lambda_ = 0
tau = 30
a = 1/300

myCenters = 0.5*( CLASS[1:C] + CLASS[2:C+1] )
Mk_vec = np.linspace(0, MMAX, N)

# TODO: revise implementation details for performance (speed + memory)

iva = 0
delta_distrib = np.ones(1,N)*(M/N)
Sdelta = np.zeros(1,T)

# TODO: develop tests + test code
# TODO: develop a visualization module 
#       (this can be used for Wright's model too)
# TODO: modularize to use delta/uniform initial distribution
# TODO: add logic to enable/disable tax model

# TODO: rethink design to implement well-being analysis
OcaD = [None]*T
OcbD = [None]*T
B = [None]*T

# IDK IF THIS RUNS PROPERLY OR NOT, 
# THIS IS JUST A ROUGH TRANSLATION OF MATLAB CODE;
# DO NOT RUN YET!!1!1!!!

for t in range(1, T+1):
    # histogram visualization used to go here
    if t % tau == 0:
        mD = mD + iva/N
        iva = 0

    for j in range(1, N+1):
        l = np.random.randint(1, N)
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

    Sdelta[t] = N*np.log(N) - np.sum(hist*np.log(hist))

    # Well being functions
    OcaD[t] = sum(a*myCenters*hist)
    OcbD[t] = sum(hist*(1-np.exp(-a*myCenters)))
    
    VEm = 0
    for u in range(1, C+1):
        VEm = VEm + myCenters(u)*hist[u]/N
    B[t] = VEm
    VEm = 0

# TODO: Translate Entropy Plots
    
# TODO: Well-being Boltzmann-Gibbs distribution
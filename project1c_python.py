import matplotlib.pyplot as plt
import numpy as np
import os

# Initialte list with different n-values to be evaluated
n_list = [10, 100, 1000, 10000, 1e5]

epsilon_max = np.zeros(len(n_list))
h_log10 = np.zeros(len(n_list))

k = 0
# Loop over n-values
for n in n_list:
    # Read corresponding txt file
    filename = "data/solution_n%d" % n
    f = open(filename+".txt", "r")

    # Generate arrays
    x = np.zeros(n+2)
    v = np.zeros(n+2)

    # Read data from file to arrays
    i = 0
    for line in f:
        column = line.split()
        x[i] = column[0]
        v[i] = column[1]
        i = i+1


    # Generate analytical solution
    u = 1-(1-np.exp(-10))*x - np.exp(-10*x)

    #Slice arrays to get rid of zero element:
    v = v[1:-1]
    u = u[1:-1]

    # Compute relative error:
    epsilon = np.log10(np.abs((v - u)/u));
    epsilon_max[k] = np.max(epsilon)

    # Compute h and h_log10:
    h = 1.0/(n+1)
    h_log10[k] = np.log10(h);
    k += 1

# Plot and save .eps file:
plt.plot(h_log10, epsilon_max)
plt.xlabel('log(h)')
plt.ylabel('max(epsilon)')
plt.grid('on')
plt.savefig("relative_error_python.eps")

# Move .eps file to data filder:
os.system("mv relative_error_python.eps data")

import matplotlib.pyplot as plt
import numpy as np
import os

# Initialte list with different n-values to be evaluated
n_list = [10, 100, 1000]

# Loop over n-values
for n in n_list:
    # Call cpp code to compute solution and generate txt file with answer
    call = './../Project1-build/Project1 %d' % n
    os.system(call)

    # Read corresponding txt file
    filename = "solution_n%d" % n
    f = open(filename+".txt", "r")

    # Generate arrays
    x = np.zeros(n)
    v = np.zeros(n)

    # Read data from file to arrays
    i = 0
    for line in f:
        column = line.split()
        x[i] = column[0]
        v[i] = column[1]
        i = i+1

    # Generate analytical solution
    u = 1-(1-np.exp(-10))*x - np.exp(-10*x)

    # Plot command:
    plt.plot(x,v, x,u)
    plt.xlabel('x')
    plt.ylabel('u(x)')
    plt.title('n = %d' % n)
    plt.legend(['Approximation v(x)', 'Analytical u(x)'])
    plt.grid('on')
    plt.savefig(filename+".eps")
    plt.hold(False)

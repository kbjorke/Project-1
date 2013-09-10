import matplotlib.pyplot as plt
import numpy as np
import os

# Initialte list with different n-values to be evaluated
n_list = [10, 100, 1000, 10000, 1e5]

# Loop over n-values
for n in n_list:
    n = int(n)
    call = './../Project1-build/Project1 -s %d' % n
    os.system(call)

    # Read corresponding txt file
    filename = "solution_n%d" % n
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
    
    # CLose and move file to data folder:
    f.close()
    os.system("mv %s.txt data" % filename)

    # Generate analytical solution
    u = 1-(1-np.exp(-10))*x - np.exp(-10*x)
    
    # Plot and save .eps file:
    plt.plot(x,v, x,u)
    plt.xlabel('x')
    plt.ylabel('u(x)')
    plt.title('n = %d' % n)
    plt.legend(['Approximation v(x)', 'Analytical u(x)'])
    plt.grid('on')
    plt.savefig(filename+".eps")
    plt.hold(False)

    # Move .eps file to data filder:
    os.system("mv %s.eps data" % filename)

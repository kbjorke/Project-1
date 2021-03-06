"""
Python-script to run the cpp program Project1 to time solving
the problem with the LU decompostion method and the optimized
method. Also uses the cpp program Project1 to find relative error 
for solutions of LU decomposition method by a list of n-values and 
plots the results, in a similar way as the python-script 
'project1c.py'

Usage:  $~ python project1d.py 
        $~ python project1d.py 1 2
        $~ python project1d.py 1 3 40


Made by: Kristian Bjørke
"""
import matplotlib.pyplot as plt
import numpy as np
import os
import sys

n_list_time = [10,100,1000]

for n in n_list_time:
    os.system("./../Project1-build/Project1 -t %d" %n)

# Initialte list with different n-values to be evaluated
option = len(sys.argv) 
if option == 1:
    # Default list when no commandline argument given
    n_list = [10, 100, 1000]
elif option == 3:
    start = float(sys.argv[1])
    end = float(sys.argv[2])
    n_log_list = range(start, end)
    n_list = [10**power for power in n_log_list]
elif option == 4:
    start = float(sys.argv[1])
    end = float(sys.argv[2])
    resolution = int(sys.argv[3])
    n_list = np.logspace(start, end, resolution)


# Loop over n-values
for n in n_list:
    if int(n) > 1e3:
        print "n = %e to high value." % int(n)
        break
    call = "./../Project1-build/Project1 -eLU %d" %int(n)
    os.system(call)


# Read file containing relative errors
f = open("relative_error.txt", "r")

# Generate lists
h_log10 = []
epsilon_max = []

# Read data from file to lists 
for line in f:
    if line.strip():
        column = line.split()
        h_log10.append(column[1])
        epsilon_max.append(column[2])

# Close and move file
f.close()        
os.system("mv relative_error.txt data")

# Plot and save .eps file:
plt.plot(h_log10, epsilon_max)
plt.xlabel('log(h)')
plt.ylabel('max(epsilon)')
plt.grid('on')
plt.savefig("relative_error_LU.eps")
plt.show()

# Move .eps file to data filder:
os.system("mv relative_error_LU.eps data")

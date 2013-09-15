"""
Python-script to run the cpp program Project1 to get the relative
error for solutions made by the optimized method for a set of
n-values. This set is decided by user through commadline arguments.
The resulting set of relative error is read from file and ploted.
File containing the list of errors and a .eps file og the plot 
are moved to folder /data.

Usage:  $~ python project1c.py         (n_list = [10, 100, 1000, 10000, 1e5])
        $~ python project1c.py 1 3     (n_list = [10, 100, 1000])
        $~ python project1c.py 1 7 100 (n_list = np.logspace(1, 7, 100))

Made by: Kristian Bjørke
"""

import matplotlib.pyplot as plt
import numpy as np
import os
import sys

# Initialte list with different n-values to be evaluated
option = len(sys.argv) 
if option == 1:
    # Default list for no commandline argument
    n_list = [10, 100, 1000, 10000, 1e5]
elif option == 3:
    start = int(sys.argv[1])
    end = int(sys.argv[2])
    n_log_list = range(start, end)
    n_list = [10**power for power in n_log_list]
elif option == 4:
    start = float(sys.argv[1])
    end = float(sys.argv[2])
    resolution = int(sys.argv[3])
    n_list = np.logspace(start, end, resolution)


# Loop over n-values
for n in n_list:
    if int(n) > 1e8:
        print "n = %e to high value." % int(n)
        break
    call = "./../Project1-build/Project1 -e %d" %int(n)
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
plt.savefig("relative_error.eps")
plt.show()

# Move .eps file to data filder:
os.system("mv relative_error.eps data")

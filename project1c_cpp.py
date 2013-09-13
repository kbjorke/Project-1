import matplotlib.pyplot as plt
import numpy as np
import os
import sys

# Initialte list with different n-values to be evaluated
option = len(sys.argv) 
if option == 1:
    n_list = [10, 100, 1000, 10000, 1e5]
elif option == 3:
    n_log_list = range(int(sys.argv[1]), int(sys.argv[2]))
    n_list = [10**power for power in n_log_list]
elif option == 4:
    n_list = np.logspace(int(sys.argv[1]), int(sys.argv[2]),
                                                       int(sys.argv[3]))


# Loop over n-values

for n in n_list:
    if int(n) > 1e7:
        print "n = %e to high value." % int(n)
        break
    call = "./../Project1-build/Project1 -e %d" %int(n)
    os.system(call)



# Read file containing relative errors
filename = "data/relative_error"
f = open(filename+".txt", "r")

# Generate lists
h_log10 = []
epsilon_max = []

# Read data from file to lists 
for line in f:
    if line.strip():
        column = line.split()
        h_log10.append(column[1])
        epsilon_max.append(column[2])

f.close()        

# Plot and save .eps file:
plt.plot(h_log10, epsilon_max)
plt.xlabel('log(h)')
plt.ylabel('max(epsilon)')
plt.grid('on')
plt.savefig("relative_error_cpp.eps")
plt.show()

# Move .eps file to data filder:
os.system("mv relative_error_cpp.eps data")

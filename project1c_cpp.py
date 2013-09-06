import matplotlib.pyplot as plt
import numpy as np
import os

# Initialte list with different n-values to be evaluated
n_list = [10, 100, 1000, 10000, 1e5]

# Loop over n-values

for n in n_list:
    call = "./../Project1-build/Project1 %d 2" %n
    os.system(call)

os.system("mv relative_error.txt data")


# Read file containing relative errors
filename = "data/relative_error"
f = open(filename+".txt", "r")

# Generate arrays
h_log10 = []
epsilon_max = []

# Read data from file to lists 
for line in f:
    if line.strip():
        column = line.split()
        h_log10.append(column[1])
        epsilon_max.append(column[2])

# Plot and save .eps file:
plt.plot(h_log10, epsilon_max)
plt.xlabel('log(h)')
plt.ylabel('max(epsilon)')
plt.grid('on')
plt.savefig("relative_error_cpp.eps")

# Move .eps file to data filder:
os.system("mv relative_error_cpp.eps data")

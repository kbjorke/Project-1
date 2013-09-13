import os

n_list = [10, 100, 1000]

for n in n_list:
    os.system("./../Project1-build/Project1 -t %d" %n)

from scitools.std import *

n = 100

f = open("../Project1-build/solution.txt", "r")

x = zeros(100)
u = zeros(100)

i = 0
for line in f:
    column = line.split()
    x[i] = column[0]
    u[i] = column[1]
    i = i+1

plot(x,u)

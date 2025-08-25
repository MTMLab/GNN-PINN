#!/usr/bin/env python
import math

# Open the file and read the number of solvent molecules
with open('number_solvent.txt', 'r') as file:
    num_solvent = int(file.readline())

# Start with the cubic root of the number of solvent molecules
start = round(num_solvent ** (1/3))

# Try to find three factors
for i in range(start, 0, -1):
    if num_solvent % i == 0:
        temp = num_solvent // i
        for j in range(start, 0, -1):
            if temp % j == 0:
                k = temp // j
                factors = [i, j, k]
                break
        if 'factors' in locals():
            break

# Sort the factors in descending order
factors.sort(reverse=True)

# Open a new file and write the factors
with open('factors.txt', 'w') as file:
    file.write(' '.join(map(str, factors)) + '\n')

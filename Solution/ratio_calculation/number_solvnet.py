#!/usr/bin/env python
# Open the files and read the masses
with open('extracted_mass.txt', 'r') as file:
    mass_A = float(file.readline())
with open('extracted_solvent_mass.txt', 'r') as file:
    mass_B = float(file.readline())

# Calculate the number of solvent molecules and round to the nearest integer
num_solvent = round((90*mass_A) / (10*mass_B))

# Print the result
print(num_solvent)

with open('number_solvent.txt', 'w') as file:
    file.write(str(num_solvent) + '\n')



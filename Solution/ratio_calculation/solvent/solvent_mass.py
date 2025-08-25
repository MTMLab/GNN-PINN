#!/usr/bin/env python
# Open the file
with open('solvent_mass.txt', 'r') as file:
    lines = file.readlines()

# Get the last line
last_line = lines[-1]

# Split the line into a list of words
words = last_line.split()

# Get the second word (the mass value)
mass = float(words[1])

# Open a new file and write the mass value
with open('extracted_solvent_mass.txt', 'w') as file:
    file.write(str(mass) + '\n')



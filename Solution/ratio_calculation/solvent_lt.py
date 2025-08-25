#!/usr/bin/env python
# Open the file and read the factors
with open('factors.txt', 'r') as file:
    factors = list(map(int, file.readline().split()))

# Define the moltemplate code with the read factors
moltemplate_code = f"""import solvent.lt
write_once("Data Boundary") {{
0.0    {15*factors[0]} xlo xhi
0.0    200 ylo yhi
0.0    200 zlo zhi
}}
solvents = new solvent[{factors[0]}].move(10, 0, 0)
              [{factors[1]}].move(0, {200/factors[1]}, 0)
              [{factors[2]}].move(0, 0, {200/factors[2]})
solvents[*][*][*].move(4.0, 4.0, 4.0)"""

# Open a new file and write the updated moltemplate code
with open('system_solvent.lt', 'w') as file:
    file.write(moltemplate_code)


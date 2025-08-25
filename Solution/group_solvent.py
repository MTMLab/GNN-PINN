#!/usr/bin/env python
with open('system_solvent.data', 'r') as f:
    lines = f.readlines()

start = False
with open('system_solvent_group.data', 'w') as f:
    for line in lines:
        parts = line.split()
        if 'Atoms  # full' in line:
            start = True
            f.write(line)
        elif start and len(parts) == 7:  # Assuming you have atom-ID, molecule-ID, atom-type, charge, x, y, and z
            parts[1] = '1'  # replace the molecule-ID with your new tag
            f.write(' '.join(parts) + '\n')
        else:
            f.write(line)


import sys


parm_bond = dict()
parm_angle = dict()
parm_dihedral = dict()
parm_improper = dict()
f = open("gaff.lt")
gaff = f.readlines()
for i in range(269, 1076):
    parm_bond[gaff[i][21:].split()[0]] = gaff[i][21 + len(gaff[i][21:].split()[0]):-1]
for i in range(1889, 6137):
    parm_angle[gaff[i][23:].split()[0]] = gaff[i][23 + len(gaff[i][23:].split()[0]):-1]
for i in range(10391, 11032):
    parm_dihedral[gaff[i][29:].split()[0]] = gaff[i][29 + len(gaff[i][29:].split()[0]):-1]
for i in range(11754, 11792):
    parm_improper[gaff[i][20:].split()[0]] = gaff[i][20:-1]
f.close()

arg = sys.argv[1]
a = open(arg, 'r').readlines()
index = [0, 0, 0, 0, 0]
for i in range(len(a)):
    if 'BOND' in a[i]:
        index[0] = i
    if 'ANGLE' in a[i]:
        index[1] = i
    if 'DIHE' in a[i]:
        index[2] = i
    if 'IMPROPER' in a[i]:
        index[3] = i
    if 'NONBON' in a[i]:
        index[4] = i
bond = a[index[0] + 1:index[1] - 1]
bond = [["".join(i[:5].split()), "".join(i[36:41].split())] for i in bond]
bo = [i[0].split("-") for i in bond]
angle = a[index[1] + 1:index[2] - 1]
angle = [["".join(i[:8].split()), "".join(i[40:48].split())] for i in angle]
an = [i[0].split("-") for i in angle]
dihedral = a[index[2] + 1:index[3] - 1]
dihedral = [["".join(i[:11].split()), "".join(i[68:79].split())] for i in dihedral]
di = [i[0].split("-") for i in dihedral]
improper = a[index[3] + 1:index[4] - 1]
improper = ["".join(i[:11].split()) for i in improper]
im = [i.split("-") for i in improper]
a = open(sys.argv[2], 'w')
print('write_once("In Settings") {', file=a)
for i in bond:
    try:
        print(f"    bond_coeff @bond:{i[0] + parm_bond[i[1]]}", file=a)
    except:
        pass
for i in angle:
    try:
        print(f"    angle_coeff @angle:{i[0] + parm_angle[i[1]]}", file=a)
    except:
        pass
for i in dihedral:
    try:
        print(f"    dihedral_coeff @dihedral:{i[0] + parm_dihedral[i[1]]}", file=a)
    except:
        pass
for i in improper:
    try:
        print(f"    improper_coeff @improper:{i} cvff 0.3667 -1 2", file=a)
    except:
        pass
print("}", file=a)

print('write_once("Data Bonds By Type") {', file=a)
for i in range(len(bo)):
    print(f"    @bond:{bond[i][0]} @atom:{bo[i][0]} @atom:{bo[i][1]}", file=a)
print("}", file=a)

print('write_once("Data Angles By Type") {', file=a)
for i in range(len(an)):
    print(f"    @angle:{angle[i][0]} @atom:{an[i][0]} @atom:{an[i][1]} @atom:{an[i][2]}", file=a)
print("}", file=a)

print('write_once("Data Dihedrals By Type") {', file=a)
for i in range(len(di)):
    print(f"    @dihedral:{dihedral[i][0]} @atom:{di[i][0]} @atom:{di[i][1]} @atom:{di[i][2]} @atom:{di[i][3]}", file=a)
print("}", file=a)

print('write_once("Data Impropers By Type") {', file=a)  # (gaff_imp.py) are not used here
for i in range(len(im)):
    print(f"    @improper:{improper[i]} @atom:{im[i][2]} @atom:{im[i][1]} @atom:{im[i][0]} @atom:{im[i][3]}", file=a)
print("}", file=a)

a.close()

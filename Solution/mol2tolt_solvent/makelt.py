import sys

a1 = open(sys.argv[1], 'r').readlines()
index = [0, 0, 0]
for i in range(len(a1)):
    if '@<TRIPOS>ATOM' in a1[i]:
        index[0] = i
    if '@<TRIPOS>BOND' in a1[i]:
        index[1] = i
    if '@<TRIPOS>SUBSTRUCTURE' in a1[i]:
        index[2] = i
a = a1[index[0] + 1:index[1]]
cation = [
    [a[i].split()[5], float(a[i].split()[-1]), float(a[i].split()[2]), float(a[i].split()[3]), float(a[i].split()[4])]
    for i in range(len(a))]
a = a1[index[1] + 1:index[2]]
cation_bond = [[i.split()[1], i.split()[2]] for i in a]

a = open(sys.argv[2], 'w')
print('import "gaff.lt"', file=a)
print("solvent inherits GAFF{", file=a)
# print(f"{sys.argv[2]}""aaaa" " inherits GAFF{", file=a)   cation이 아니라 분자이름으로 바뀜
print(f"#include '{sys.argv[3]}'", file=a)
print('write("Data Atoms"){', file=a)
for i in range(len(cation)):
    print("    $atom:atom" + str(i+1), "$mol", "@atom:" + cation[i][0], cation[i][1], cation[i][2], cation[i][3],
          cation[i][4], file=a)
    # print("    $atom:atom"+str(i+1),"$mol","@atom:"+cation[i][0],0.0,cation[i][2],cation[i][3],cation[i][4],file=a)
print("}", file=a)
print("write('Data Bond List') {", file=a)
for i in range(len(cation_bond)):
    print("    $bond:bond" + str(i + 1), "$atom:atom" + cation_bond[i][0], "$atom:atom" + cation_bond[i][1], file=a)
print("}", file=a)
print("}", file=a)
a.close()

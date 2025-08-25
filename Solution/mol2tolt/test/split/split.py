import sys
f=open(sys.argv[1],"r").readlines()
ids=[id  for id,i in enumerate(f) if "@<TRIPOS>MOLECULE" in i]
ids_=ids+[len(f)]
for molid_id,mol_id in enumerate(ids):
    with open(f"{sys.argv[1].split('.')[0]}_{molid_id}.mol2","w") as ff:
        print("".join(f[mol_id:ids_[molid_id+1]]+["@<TRIPOS>SUBSTRUCTURE"]),file=ff)

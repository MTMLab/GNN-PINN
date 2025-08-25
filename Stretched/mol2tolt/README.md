# Install `Ambertools`

conda install in current **ENV**

```shell
conda install -c conda-forge ambertools
```

or create a new **ENV** and enter

```shell
conda create -n ambertools -c conda-forge ambertools
conda activate ambertools
```

# convert mol2 to lt

```shell
# example
cd test
sh ../mol2tolt.sh C12H16N4_0.mol2
```

and it will generate 3 lt files

```shell
.
├── C12H16N4_0_add.lt # this file contains recommended parameters (ambertools parmchk2) for some topologies that gaff2.1 not include
├── C12H16N4_0.lt  # this file is the molecule template
├── C12H16N4_0.mol2 
└── gaff2.1.lt     # gaff2 parameters (some bugs are fixed, don't make any change)
```

## SDF.mol2 file

for **sdf.mol2** file, it should firstly **split sdf** into mol2 files.

```
# example
cd test/split
python split.py C12H16N4.sdf.mol2
```

 


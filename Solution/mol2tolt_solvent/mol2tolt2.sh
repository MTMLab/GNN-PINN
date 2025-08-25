#!/bin/sh

SCRIPTDIR=`dirname $0`

fullname=$1
name=${fullname%.*}

cp $SCRIPTDIR/gaff.lt ./

#antechamber -i tmp.mol2 -fi mol2 -fo mol2 -o monomer.mol2 -pf y -at gaff

parmchk2 -i solvent.mol2 -f mol2 -o ${name}.frcmod -s 2
python $SCRIPTDIR/addp.py ${name}.frcmod ${name}_add.lt
python $SCRIPTDIR/makelt.py solvent.mol2 $name.lt ${name}_add.lt
#rm ${name}.frcmod monomer.mol2

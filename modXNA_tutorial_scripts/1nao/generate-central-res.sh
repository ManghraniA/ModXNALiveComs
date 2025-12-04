#!/bin/bash

WORKDIR=`pwd`

modxna="../../modxna.sh"

declare -A ResNames

ResNames["5MG"]="5PO OME DGG"
ResNames["3MC"]="DPO OME DCC"
ResNames["OMC"]="DPO OME DCC"
ResNames["OMU"]="DPO OME RUU"

for res in OMC OMU ; do

mkdir $res/
cd $res/

echo "${ResNames[$res]}" > $res.modxna.in

$modxna -i $res.modxna.in -m $res

cd $WORKDIR

done

#5' 5CG cap
for res5pr in 5MG ; do

mkdir $res5pr/
cd $res5pr/

echo "${ResNames[$res5pr]}" > $res5pr.modxna.in

$modxna -i $res5pr.modxna.in -m $res5pr --5cap

cd $WORKDIR

done

#3' 3EC cap
for res3pr in 3MC ; do

mkdir $res3pr/
cd $res3pr/

echo "${ResNames[$res3pr]}" > $res3pr.modxna.in

$modxna -i $res3pr.modxna.in -m $res3pr --3cap

cd $WORKDIR

done

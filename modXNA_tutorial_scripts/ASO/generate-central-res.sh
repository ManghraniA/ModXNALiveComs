#!/bin/bash

WORKDIR=`pwd`

modxna="/home/bergonzoc/GitHub/modXNA.dan/modxna.sh"

declare -A ResNames

ResNames["5CG"]="5PO CET DGG"
ResNames["CEG"]="RPO CET RGG"
ResNames["MEC"]="RPO CET M5C"
ResNames["MOG"]="RPO MOE RGG"
ResNames["MOA"]="RPO MOE RAA"
ResNames["POG"]="PS1 MOE RGG"
ResNames["POA"]="RPO MOE RAA"
ResNames["MSC"]="PS1 CET M5C"
ResNames["OMU"]="PS1 OME RUU"
ResNames["3EC"]="PS1 CET M5C"


for res in CEG MEC MOG MOA POG POA MSC OMU ; do

mkdir $res/
cd $res/

echo "${ResNames[$res]}" > $res.modxna.in

$modxna -i $res.modxna.in -m $res

cd $WORKDIR

done

#5' 5CG cap
for res5pr in 5CG ; do

mkdir $res5pr/
cd $res5pr/

echo "${ResNames[$res5pr]}" > $res5pr.modxna.in

$modxna -i $res5pr.modxna.in -m $res5pr --5cap

cd $WORKDIR

done

#3' 3EC cap
for res3pr in 3EC ; do

mkdir $res3pr/
cd $res3pr/

echo "${ResNames[$res3pr]}" > $res3pr.modxna.in

$modxna -i $res3pr.modxna.in -m $res3pr --3cap

cd $WORKDIR

done

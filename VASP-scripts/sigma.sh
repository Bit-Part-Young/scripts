#!/bin/bash

: '
检查 entropy T*S 是否小于 1 meV/atom

reference: https://github.com/tamaswells/VASP_script/blob/master/sigma.sh
'

entropy=$(grep 'entropy T\*S' OUTCAR | tail -1 | awk '{print $5}')
natoms=$(grep 'NIONS' OUTCAR | tail -1 | awk '{print $12}')

entropy_average=$(echo "scale=5; ${entropy}/${natoms}" | bc)
value_abs=${entropy_average#-}

if grep -q 'SIGMA' INCAR; then
  sigma=$(grep SIGMA INCAR | grep -Eo '([0-9.]+)' | head -1)
else
  sigma=0.2
fi

flag=$(echo "$value_abs < 0.001" | bc)

if [ "${flag}" -eq 1 ]; then
  echo "Your SIGMA $sigma is okay; for average entropy T*S $value_abs is < 0.001 eV/atom."
else
  echo "Your SIGMA $sigma is BAD; for average entropy T*S $value_abs is > 0.001 eV/atom. Try to decrease SIGMA!"
fi

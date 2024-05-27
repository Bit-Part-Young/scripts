#!/bin/bash

# To determinr whether entropy T*S is less than 1 mev
# reference: https://github.com/tamaswells/VASP_script/blob/master/sigma.sh

entropy_ts=$(grep 'entropy T\*S' OUTCAR | tail -n 1 | awk '{print $5}')
natoms=$(sed -n 7p POSCAR | tr -d '\r' | awk '{for(i=1; i<=NF; i++) sum+=$i; print sum}')

energy_ts_average=$(echo "scale=5; ${entropy_ts}/${natoms}" | bc)
value_abs=${energy_ts_average#-}

sigma=$(echo "SIGMA = 0.1 (for metal: 0.2)" | grep -Eo '([0-9.]+)' | head -1)

flag=$(echo "$value_abs < 0.001" | bc)

if [ "${flag}" -eq 1 ]; then
	echo "Your SIGMA $sigma is okay; for average T*S $value_abs is < 0.001 eV/atom."
else
	echo "Your SIGMA $sigma is BAD; for average T*S $value_abs is > 0.001 eV/atom. Try to decrease SIGMA!"
fi

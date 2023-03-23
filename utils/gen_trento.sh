#!/bin/sh

for i in $@; do
  dx=$(echo 20 / $i | bc -l)
  for b in 0 3 7; do
    for sig in 4.23 6.4 7.32; do
      mkdir -p sig$sig/b$b/e$i
      trento Pb Pb --random-seed 1 --b-min $b --b-max $b --cross-section $sig --grid-max 10 --grid-step $dx 100 -o sig$sig/b$b/e$i
    done
  done
done

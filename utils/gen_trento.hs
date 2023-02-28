#!/bin/sh

for i in $@; do
  dx=$(echo 20 / $i | bc -l)
  trento Pb Pb --random-seed 1 --b-min 7 --b-max 7 --grid-max 10 --grid-step $dx 100 -o e$i
done

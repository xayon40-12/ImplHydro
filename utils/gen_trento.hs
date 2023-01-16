#!/bin/sh

trento Pb Pb --random-seed 1 --grid-max 10 --grid-step 0.2 100 -o e100
trento Pb Pb --random-seed 1 --grid-max 10 --grid-step 0.1 100 -o e200
trento Pb Pb --random-seed 1 --grid-max 10 --grid-step 0.0666666666666 100 -o e300

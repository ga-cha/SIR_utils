#!/bin/bash

for i in {11..999}
do
    python voxelwise_brainsmash.py $i | tee -a output.log
done
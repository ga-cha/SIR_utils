#!/bin/bash

out=slurm_vox_bs.log
py_sir=/fs03/kg98/gchan/conda_envs/sir/bin/python

for i in {1..999}
do
    sbatch --job-name="$i brainsmash" --output="$out" <<EOL
#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=12
#SBATCH --time=24:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=gabriella.chan@monash.edu

$py_sir voxelwise_brainsmash.py $i
EOL
done
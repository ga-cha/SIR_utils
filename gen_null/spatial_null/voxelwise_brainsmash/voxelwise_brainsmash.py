# 10/02/25 
# Gabriella Chan
# gabriella.chan@moansh.edu

# a py workflow for voxelwise_brainsmash.ipynb
# run the ipynb prior to this to save masked gmv and distance matrix and select parameters

import nibabel as nib
import pandas as pd
import numpy as np
from brainsmash.workbench import volume
from brainsmash.mapgen.sampled import Sampled
import sys

if len(sys.argv) < 2:
    print("Usage: python voxelwise_brainsmash.py map_id")
i = sys.argv[1]

# Load atlas and set mask
s132_img = nib.load('/fs03/kg98/gchan/Atlases/Tian/Schaefer_Tian/reordered/Schaefer2018_100Parcels_' +
    '7Networks_order_Tian_Subcortex_S2_MNI152NLin6Asym_1.5mm_reordered.nii.gz')
atlas = s132_img.get_fdata()
roi_i = 1
roi_j = 66
mask = (atlas >= roi_i) & (atlas <= roi_j)

data_dir = './data/lh/'
out_dir = './results/lh/'
brain_map = data_dir + 'lme_betas.txt'
filenames = {"D": data_dir + "distmat.npy",
            "index": data_dir + "index.npy"
            }

kwargs = {'ns': 500,
          'knn': 1500,
          'pv': 60,
          'kernel': 'gaussian',
          'n_jobs': 12
          }
          
gen = Sampled(x=brain_map, D=filenames['D'], index=filenames['index'], **kwargs)
surrogate_maps = gen(n=1)

# fill in the surrogate map values
null_map = np.zeros_like(atlas, dtype=float)
null_map[mask] = surrogate_maps

# Combine output with s132_img header and save
null_img = nib.Nifti1Image(null_map, s132_img.affine, s132_img.header)
nib.save(null_img, out_dir + "null_map_" + i + ".nii.gz")
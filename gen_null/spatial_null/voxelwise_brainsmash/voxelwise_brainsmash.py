# 10/02/25 
# Gabriella Chan
# gabriella.chan@moansh.edu

# a py workflow for voxelwise_brainsmash.ipynb
# run the ipynb prior to this to save masked gmv and distance matrix and select parameters

import nibabel as nib
import numpy as np
from brainsmash.mapgen.sampled import Sampled
import sys
import os
import pickle
from sklearn.utils import check_random_state

import time

def get_map_generator(structure, kwargs):
    data_dir = structure["data_dir"]
    brain_map = data_dir + 'lme_betas.txt'
    filenames = {
        "D": data_dir + "distmat.npy",
        "index": data_dir + "index.npy"
    }

    try:
        with open(data_dir + 'gen.pkl', 'rb') as f:
            gen = pickle.load(f)
            # reset random state in generator
            gen._rs = check_random_state(None)
    except FileNotFoundError:
        print("Map generator not found. Creating map generator for " + structure["name"])
        start_time = time.time()

        gen = Sampled(x=brain_map, D=filenames['D'], index=filenames['index'], **kwargs)
        with open(data_dir + 'gen.pkl', 'wb') as f:
            pickle.dump(gen, f)
        print(f"Map generator created in {time.time() - start_time:.2f} seconds")


    return gen

def save_null_map(surrogate_maps, structure, atlas_img, i):
    mask = structure["mask"]
    if not os.path.exists(structure["out_dir"]):
        os.makedirs(structure["out_dir"])

    null_map = np.zeros(atlas_img.shape)
    null_map[mask] = surrogate_maps

    null_img = nib.Nifti1Image(null_map, atlas_img.affine, atlas_img.header)
    nib.save(null_img, structure["out_dir"] + "null_map_" + i + ".nii.gz")

def combine_null_maps(i):
    null_cx_img = nib.load('./results/l_cx/null_map_' + i + '.nii.gz')
    null_subcx_img = nib.load('./results/l_subcx/null_map_' + i + '.nii.gz')
    null_cx = null_cx_img.get_fdata()
    null_subcx = null_subcx_img.get_fdata()

    null_combined = null_cx + null_subcx
    nib.save(nib.Nifti1Image(null_combined, null_cx_img.affine, null_cx_img.header), './results/lh/null_map_' + i + '.nii.gz')

if __name__ == "__main__":

    if len(sys.argv) < 2:
        print("Usage: python voxelwise_brainsmash.py map_id")
        sys.exit(1)
    i = sys.argv[1]

    # Load atlas
    atlas_img = nib.load('/fs03/kg98/gchan/Atlases/Tian/Schaefer_Tian/reordered/Schaefer2018_100Parcels_' +
        '7Networks_order_Tian_Subcortex_S2_MNI152NLin6Asym_1.5mm_reordered.nii.gz')
    atlas = atlas_img.get_fdata()

    structures = [
        {
            "name": "l_cx",
            "roi_i": 1,
            "roi_j": 50,
            "data_dir": './data/l_cx/',
            "out_dir": './results/l_cx/'
        },
        {
            "name": "l_subcx",
            "roi_i": 51,
            "roi_j": 66,
            "data_dir": './data/l_subcx/',
            "out_dir": './results/l_subcx/'
        }
    ]

    kwargs = {
        'ns': 500,
        'knn': 1500,
        'pv': 60,
        'kernel': 'gaussian',
        'n_jobs': 13
    }

    for structure in structures:
        structure["mask"] = (atlas >= structure["roi_i"]) & (atlas <= structure["roi_j"])

        gen = get_map_generator(structure, kwargs)
        
        start_time = time.time()
        surrogate_maps = gen(n=1)
        print(f"Null maps for {structure['name']} generated in {time.time() - start_time:.2f} seconds")

        save_null_map(surrogate_maps, structure, atlas_img, i)

    combine_null_maps(i)
    print(f"Null map {i} complete")
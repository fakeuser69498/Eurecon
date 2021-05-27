# Eurecon

Equidistant and Uniform Data Augmentation for 3D Objects.

Python library for 3D data augmentation, based on Eurecon algorithm. 

- The library works with a variety of 3D file formats, including point cloud formats (.XYZ, .PTS, .PCD), polygon mesh formats (.STL, .OFF, .OBJ, .PLY, .GLTF) and biomolecular formats (.MOL2, .PDB, .XYZ).
- Eurecon is computationally efficient taking ~0.1 seconds to generate 1,000 samples  of an object of 1,000 3D points. 
- Based on numpy, open3d and OpenBabel.
- Simple, flexible API that allows the library to be used in any machine learning pipeline.
- Large, diverse set of transformations based on the amount of input tessellation axes.
- Easy to extend the library to wrap around other libraries.
- Easy to extend to other tasks.



## Requirements:
- Python 3.8+
- click==7.1.2+
- click-plugins==1.1.1+
- numpy==1.17.2+
- open3d==0.10.0.0+
- openbabel==3.1.1.1+
- tqdm==4.46.1+


## How to use:

Easiest way to run the augmentation via Eurecon is running eurerun.py with the following command line:

python3.8 examples/eurerun.py -r <desired_RMSD> -p <desired_partition_parameter> -a <path/to/tessellation/axes/file> -rr <enabling_or_disabling_relmin> -rrr <desired_relmin_ratio> -i <path/to/input/files> -o <path/to/the/output/directory>

## Important notes:

- You have to create the output directory before running the algorithm in case it was not created beforehand
- Using augmentation on biomolecular file formats requires having OpenBabel (https://github.com/openbabel/openbabel/releases/tag/openbabel-3-1-1) installed; also please consider the fact that running Eurecon would effectively delete most of the metadata, leaving only an array of 3D coordinates and atom types

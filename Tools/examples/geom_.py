import trimesh
import sys
sys.path.insert(1, r'..\Tools') ## Replace with proper packaging one day
print(sys.path)
import Geometry.from_stl_to_voxel as stl2vox
import numpy as np

# Step 1: Create Basic Geometric Shapes
cube = trimesh.creation.box(extents=(10, 10, 10))  # Create a cube of size 1x1x1
sphere = trimesh.creation.icosphere(radius=6)  # Create a sphere of radius 0.6

# Step 2: Perform a Boolean Operation (Union, Intersection, Difference)
# Example: Difference between cube and sphere
difference_mesh = cube.difference(sphere)
intersection_mesh = cube.intersection(sphere)
# Step 3: Export the resulting mesh as an STL file
output_file = "difference_cube_sphere.stl"  # Output STL file path

# Exporting the mesh to STL format
difference_mesh.export(output_file)
print(f"Exported the resulting mesh to {output_file}")
output_file = "inner_core.stl"  # Output STL file path

# Exporting the mesh to STL format
intersection_mesh.export(output_file)
print(f"Exported the resulting mesh to {output_file}")

# 4 set the volume resolution
resolution=[0.5,0.5,1]

# 5 import the largest encompassing surface and store the translation data
stl_mesh, transp = stl2vox.load_and_translate_stl_to_positive_quadrant('difference_cube_sphere.stl')
stl2vox.display_stl_with_matplotlib(stl_mesh)
#grid=stl2vox.VoxelGrid(resolution,np.array([0,0,0]),np.array([15,15,15]))

# 6 voxelize this volume
voxel_grid = stl2vox.create_voxel_grid(stl_mesh, resolution)

# 7 create grid for the second surface from the larger one
min_bound, max_bound = stl2vox.get_bounds(stl_mesh)
voxel_grid2 =stl2vox.VoxelGrid(resolution,min_bound,max_bound)

#8 load inner surface
stl_mesh2, _ = stl2vox.load_and_translate_stl_to_positive_quadrant('inner_core.stl',transp)

# voxelize second surface
voxel_grid2 = stl2vox.create_voxel_grid(stl_mesh2, resolution,voxel_grid2)

# 10 add together assign a second material to the inner surface
output=voxel_grid.data.astype(np.int32)
output[voxel_grid2.data]=2

# Save the voxel grid as a raw file
raw_filename = "cube_sphere_grid.raw"
mhd_filename = "cube_sphere_grid.mhd"
#save_voxel_grid_as_raw(voxel_grid, raw_filename)
stl2vox.save_voxel_data_as_raw_mhd(output, resolution,raw_filename,mhd_filename)
print(f"Voxelization complete. Voxel grid saved as {raw_filename}.")

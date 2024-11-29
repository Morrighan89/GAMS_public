import trimesh
import numpy as np
import sys
sys.path.insert(1, r'..\Tools') ## Replace with proper packaging one day
print(sys.path)
import Geometry.from_stl_to_voxel as stl2vox

# Step 1: Create Two Concentric Disks (Cylinders)
# Disk 1 (radius 1, height 0.1)
disk_1 = trimesh.creation.cylinder(radius=10.0, height=1)

# Disk 2 (radius 0.5, height 0.1)
disk_2 = trimesh.creation.cylinder(radius=5.0, height=1)

# Exporting the mesh to STL format
disk_1.export('disk_out.stl')
disk_2.export('disk_in.stl')

resolution=[0.1,0.1,1]
stl_mesh, transp = stl2vox.load_and_translate_stl_to_positive_quadrant('disk_out.stl')
stl2vox.display_stl_with_matplotlib(stl_mesh)
#grid=stl2vox.VoxelGrid(resolution,np.array([0,0,0]),np.array([15,15,15]))
voxel_grid = stl2vox.create_voxel_grid(stl_mesh, resolution)

min_bound, max_bound = stl2vox.get_bounds(stl_mesh)
voxel_grid2 =stl2vox.VoxelGrid(resolution,min_bound,max_bound)
stl_mesh2, _ = stl2vox.load_and_translate_stl_to_positive_quadrant('disk_in.stl',transp)
voxel_grid2 = stl2vox.create_voxel_grid(stl_mesh2, resolution,voxel_grid2)

output=voxel_grid.data.astype(np.int32)
output[voxel_grid2.data]=2
# Save the voxel grid as a raw file
raw_filename = "2disk_grid.raw"
mhd_filename = "2disk_grid.mhd"
#save_voxel_grid_as_raw(voxel_grid, raw_filename)
stl2vox.save_voxel_data_as_raw_mhd(output, resolution,raw_filename,mhd_filename)
print(f"Voxelization complete. Voxel grid saved as {raw_filename}.")
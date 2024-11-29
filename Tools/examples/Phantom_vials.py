
import numpy as np
import os
import sys
sys.path.insert(1, r'..\Tools') ## Replace with proper packaging one day
import Geometry.from_stl_to_voxel as stl2vox
path= r'..\..\Test_Object'

resolution=[0.4375, 1, 0.4375]

stl_mesh, transp = stl2vox.load_and_translate_stl_to_positive_quadrant(os.path.join(path,'Shell000.stl'))
stl2vox.display_stl_with_matplotlib(stl_mesh)
#grid=stl2vox.VoxelGrid(resolution,np.array([0,0,0]),np.array([15,15,15]))
voxel_grid = stl2vox.create_voxel_grid(stl_mesh, resolution)
output=voxel_grid.data.astype(np.int32)
min_bound, max_bound = stl2vox.get_bounds(stl_mesh)
for i in range(1,40):
    voxel_grid2 =stl2vox.VoxelGrid(resolution,min_bound,max_bound)
    stl_mesh2, _ = stl2vox.load_and_translate_stl_to_positive_quadrant(os.path.join(path,f'Shell{i:03d}.stl'),transp)
    voxel_grid2 = stl2vox.create_voxel_grid(stl_mesh2, resolution,voxel_grid2)
    output[voxel_grid2.data]=i+1
# Save the voxel grid as a raw file
raw_filename = "tub2_grid.raw"
mhd_filename = "tub2_grid.mhd"
#save_voxel_grid_as_raw(voxel_grid, raw_filename)
stl2vox.save_voxel_data_as_raw_mhd(output,resolution,raw_filename,mhd_filename)
print(f"Voxelization complete. Voxel grid saved as {raw_filename}.")
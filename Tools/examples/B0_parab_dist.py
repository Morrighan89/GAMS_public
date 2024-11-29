import B0_Distortion.B0_dist_calculator as b0d


def main():
    # Read grid information from the file
    grid_info = b0d.read_grid_info("2d_sphere_HIres.msh")
    
    # Define grid ranges and number of points
    x_range = (0.0+grid_info["ox"], grid_info["SizeX"]+grid_info["ox"])  # range for x-axis
    y_range = (0.0+grid_info["oy"], grid_info["SizeY"]+grid_info["oy"])  # range for y-axis
    z_range = (0.0+grid_info["oz"], grid_info["SizeZ"]+grid_info["oz"])  # range for z-axis
    nx, ny, nz = grid_info["nx"], grid_info["ny"], grid_info["nz"]  # number of grid points along each axis

    # Call the function to compute and save the data
    
    data, dx, dy, dz = b0d.compute_on_grid(lambda x, y, z: example_function_2(x, y, z, a), x_range, y_range, z_range, nx, ny, nz)
    
    
    #plot_xy_slice_surface(data, y_range, x_range, nz, z_target=0)
    b0d.plot_xy_slice(data,'B0_dist_noise5e-7', y_range, x_range, nz, z_target=0)
    
    b0d.save_to_file_as_2d_slices("output_2d_slices.dat", data, nx,ny,nz)
    
    print("Grid Information:")
    for key, value in grid_info.items():
        print(f"{key}: {value}")

    
    
if name == "__main__":
    main()
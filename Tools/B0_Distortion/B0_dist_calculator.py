import numpy as np
import matplotlib
from matplotlib.backends.backend_pgf import FigureCanvasPgf
matplotlib.backend_bases.register_backend('pdf', FigureCanvasPgf)

import matplotlib.pyplot as plt

pgf_with_custom_preamble = {
    "pgf.texsystem": "lualatex",
    "font.family": "sans-serif",  # use serif/main font for text elements
    "text.usetex": True,     # use inline math for ticks
    "pgf.rcfonts": False,    # don't setup fonts from rc parameters
    "axes.labelsize": 22,
    "axes.formatter.useoffset": True,
    "font.size": 22,
    "legend.fontsize": 20,
    "axes.titlesize": 20,           # Title size when one figure
    "xtick.labelsize": 20,
    "ytick.labelsize": 20,
    "figure.titlesize": 22,         # Overall figure title
    "pgf.preamble": r"\usepackage{fontspec} \usepackage{units} \usepackage{metalogo} \usepackage{unicode-math} \usepackage{amsmath} \setmathfont{Fira Math} \setmonofont{Libertinus Mono} "
}
matplotlib.rcParams.update(pgf_with_custom_preamble)


def compute_on_grid(f, x_range, y_range, z_range, nx, ny, nz, filename="output.dat"):
    """
    Computes the function f(x, y, z) on a grid and saves the results to an ASCII file.
    
    Parameters:
        f (function): A function of three variables, f(x, y, z).
        x_range (tuple): A tuple specifying the range of x values (min_x, max_x).
        y_range (tuple): A tuple specifying the range of y values (min_y, max_y).
        z_range (tuple): A tuple specifying the range of z values (min_z, max_z).
        nx (int): Number of grid points along the x-axis.
        ny (int): Number of grid points along the y-axis.
        nz (int): Number of grid points along the z-axis.
        filename (str): The name of the file to save the output.
    """
    # Calculate grid spacing #### CHECK THIS!!!!
    dx = (x_range[1] - x_range[0]) / (nx)
    dy = (y_range[1] - y_range[0]) / (ny)
    dz = (z_range[1] - z_range[0]) / (nz)

    data = np.zeros((nx, ny, nz))

    # Generate linearly spaced values for x, y, z
    x_vals = np.linspace(x_range[0], x_range[1], nx)
    y_vals = np.linspace(y_range[0], y_range[1], ny)
    z_vals = np.linspace(z_range[0], z_range[1], nz)

    # Compute the function value for each point on the grid
    for i, xi in enumerate(x_vals):
        for j, yi in enumerate(y_vals):
            for k, zi in enumerate(z_vals):
                data[i, j, k] = f(xi, yi, zi)  # Store computed value in the 3D array

    return data, dx, dy, dz

def save_to_file(filename, data, dx, dy, dz, x_range, y_range, z_range, nx, ny, nz):
    """
    Saves the grid metadata and data values to an ASCII file.
    
    Parameters:
        filename (str): The name of the file to save.
        data (np.ndarray): 3D array of computed function values.
        dx, dy, dz (float): Grid spacing in each dimension.
        x_range, y_range, z_range (tuple): Range for x, y, z values.
        nx, ny, nz (int): Number of grid points in each dimension.
    """
    with open(filename, "w") as file:
        # Write header information
        file.write(f"{x_range[1] - x_range[0]} {y_range[1] - y_range[0]} {z_range[1] - z_range[0]}\n")
        file.write(f"{dx} {dy} {dz}\n")
        file.write(f"{nx} {ny} {nz}\n")
        
        # Write the grid values in x, y, z, value format
        for i in range(nx):
            for j in range(ny):
                for k in range(nz):
                    # Compute x, y, z positions based on grid indices
                    xi = x_range[0] + i * dx
                    yi = y_range[0] + j * dy
                    zi = z_range[0] + k * dz
                    # Retrieve the computed value from the 3D array
                    value = data[i, j, k]
                    file.write(f"{xi} {yi} {zi} {value}\n")

    print(f"Data saved to {filename}")

# Function to save the computed values ordered by z-slices to an ASCII file
def save_to_file_by_z_slices(filename, data, dx, dy, dz, x_range, y_range, z_range, nx, ny, nz):
    """
    Saves the grid data ordered by slices along the z-axis to an ASCII file.
    
    Parameters:
        filename (str): The name of the file to save.
        data (np.ndarray): 3D array of computed function values.
        dx, dy, dz (float): Grid spacing in each dimension.
        x_range, y_range, z_range (tuple): Range for x, y, z values.
        nx, ny, nz (int): Number of grid points in each dimension.
    """
    with open(filename, "w") as file:
        # Write header information
        file.write(f"{x_range[1] - x_range[0]} {y_range[1] - y_range[0]} {z_range[1] - z_range[0]}\n")
        file.write(f"{dx} {dy} {dz}\n")
        file.write(f"{nx} {ny} {nz}\n")
        
        # Iterate over each z-slice
        for k in range(nz):
            zi = z_range[0] + k * dz  # Fixed z coordinate for this slice
            file.write(f"\n# Slice at z = {zi}\n")
            # Iterate over x and y for the current z-slice
            for i in range(nx):
                for j in range(ny):
                    xi = x_range[0] + i * dx
                    yi = y_range[0] + j * dy
                    value = data[i, j, k]
                    file.write(f"{xi} {yi} {zi} {value}\n")

    print(f"Data saved to {filename} ordered by z-slices.")

def save_to_file_as_2d_slices(filename, data,nx, ny, nz):
    """
    Saves the computed data values as a series of 2D matrices in the xy plane for each z-slice.
    
    Parameters:
        filename (str): The name of the file to save.
        data (np.ndarray): 3D array of computed function values.
        nz (int): Number of grid points along the z-axis.
    """
    with open(filename, "w") as file:
        # Iterate over each z-slice and save each as a 2D matrix
        np.savetxt(file, [[nx, ny, nz]], fmt='% 5d')
        for k in range(nz):
            #file.write(f"# Slice at z-index {k}\n")  # Header indicating the z-slice
            # Write the 2D slice (matrix in the xy plane) for the current z
            np.savetxt(file, data[:, :, k], fmt="%.12f")
            #file.write("\n")  # Add a blank line between slices for readability

    print(f"Data saved to {filename} as 2D slices in the xy plane.")

# Function to read the ASCII file and extract parameters
def read_grid_info(filename):
    """
    Reads the grid size and spacing information from an ASCII file.
    
    Parameters:
        filename (str): The name of the file to read.
    
    Returns:
        A dictionary containing SizeX, SizeY, SizeZ, dx, dy, dz, nx, ny, and nz.
    """
    with open(filename, "r") as file:
        # Read first three lines of the file
        SizeX, SizeY, SizeZ = map(float, file.readline().split())
        dx, dy, dz = map(float, file.readline().split())
        nx, ny, nz = map(int, file.readline().split())
        ox, oy, oz = map(float, file.readline().split())

    # Return the extracted values as a dictionary
    return {
        "SizeX": SizeX, "SizeY": SizeY, "SizeZ": SizeZ,
        "dx": dx, "dy": dy, "dz": dz,
        "nx": nx, "ny": ny, "nz": nz,
        "ox": ox, "oy": oy, "oz": oz
    }

def plot_xy_slice_surface(data, x_range, y_range, nz, z_target=0):
    """
    Plots a 3D surface plot of the xy-plane slice of f(x, y, z) at z = z_target.
    
    Parameters:
        data (np.ndarray): 3D array of computed function values.
        x_range (tuple): The range of x values (min_x, max_x).
        y_range (tuple): The range of y values (min_y, max_y).
        nz (int): Number of grid points along the z-axis.
        z_target (float): The value of z at which to take the slice. Default is 0.
    """
    

    # Find the z index closest to z_target
    z_vals = np.linspace(z_range[0], z_range[1], nz)
    z_index = (np.abs(z_vals - z_target)).argmin()

    # Create the meshgrid for x and y
    x_vals = np.linspace(x_range[0], x_range[1], data.shape[0])
    y_vals = np.linspace(y_range[0], y_range[1], data.shape[1])
    X, Y = np.meshgrid(x_vals, y_vals)

    # Get the data slice at the closest z_index
    Z = data[:, :, z_index]*1.e6

    # Plot the 3D surface
    fig = plt.figure(figsize=(10, 7))
    ax = fig.add_subplot(111, projection="3d")
    surf = ax.plot_surface(X, Y, Z, cmap="viridis", edgecolor="k", linewidth=0.5, antialiased=True)
    ax.set_xlabel("x")
    ax.set_ylabel("y",labelpad=10)
    
    #ax.set_zlabel("f(x, y, z=0)")
    #ax.set_title(f"z = {z_target}")
    fig.colorbar(surf, ax=ax, aspect=15,shrink=0.75, label="$\\mu$T")
    plt.tight_layout()
    plt.savefig(f'B0_dist.svg', bbox_inches='tight',dpi=1000)
    plt.savefig(f'B0_dist_gradx.pdf', bbox_inches='tight',dpi=1000,backend='pgf')
    plt.savefig(f'B0_dist_gradx.png', bbox_inches='tight',dpi=600,backend='pgf')
    plt.show()

def plot_xy_slice(data,filename, y_range, x_range, nz, z_target=0):
    """
    Plots the xy-plane slice of the function f(x, y, z) at z = z_target.
    
    Parameters:
        data (np.ndarray): 3D array of computed function values.
        y_range (tuple): The range of y values (min_y, max_y).
        x_range (tuple): The range of x values (min_x, max_x).
        nz (int): Number of grid points along the z-axis.
        z_target (float): The value of z at which to take the slice. Default is 0.
    """
    # Find the z index closest to z_target
    z_vals = np.linspace(z_range[0], z_range[1], nz)
    z_index = (np.abs(z_vals - z_target)).argmin()

    plt.figure(figsize=(8, 6))
    plt.imshow(data[:, :, z_index]*1.e6, extent=(x_range[0], x_range[1], y_range[0], y_range[1]), origin='lower', aspect='equal',vmin=0,vmax=1)
    plt.colorbar(aspect=15,shrink=0.75, label="$\\mu$T")
    plt.xlabel("x")
    plt.ylabel("y")
    
    plt.savefig(f'{filename}.pdf', bbox_inches='tight',dpi=600)
    plt.savefig(f'{filename}.pdf', bbox_inches='tight',dpi=1000,backend='pgf')
    plt.show()

# Example function to compute
def example_function(x, y, z):
    return 160*np.pi*1.e-7*x #-3.e-6

def example_function_2(x, y, z, a=1.0):
    return -a * (x**2 + (y-0.05)**2) - z + 2.e-6

def example_function_noise(x, y, z, A=1.0):
    return np.random.uniform(0, A)



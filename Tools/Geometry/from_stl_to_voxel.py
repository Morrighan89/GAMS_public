import numpy as np
import os
from stl import mesh
from concurrent.futures import ThreadPoolExecutor, as_completed
from tqdm import tqdm
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
import warnings

"""
Steps for Voxelization Using Ray Tracing

    Parse the STL File
    Read the STL file to extract the triangles (facets) that form the closed surface. Each triangle is defined by three vertices.

    Define the Voxel Grid
        Establish the bounding box of the STL geometry. This is the minimum and maximum extent of the geometry in xx, yy, and zz directions.
        Choose a voxel size and create a 3D grid that spans the bounding box. Each voxel is defined by its center point.


    Ray Tracing for Internal Voxel Detection
        Cast rays parallel to one of the axes (e.g., xx-axis) through the voxel grid.
        For each ray:
            Record the intersections with the surface triangles.
            Sort the intersection points along the ray.
            Alternate between "outside" and "inside" states as you pass through intersections:
                The region between the first and second intersection is inside, the region between the third and fourth is outside, and so on.
            Mark voxels between intersection pairs as "internal voxels."


    Output the Voxel Data
        Store the voxel data in a 3D array or export it in a format of your choice.

"""

from stl import mesh
import numpy as np

class BVHNode:
    def __init__(self, bbox, triangles=None, left=None, right=None):
        self.bbox = bbox  # Axis-aligned bounding box (min, max)
        self.triangles = triangles  # List of triangle indices
        self.left = left  # Left child
        self.right = right  # Right child

def compute_aabb(triangle_indices, stl_mesh):
    """
    Compute the AABB for a set of triangles from the STL mesh.
    """
    all_vertices = stl_mesh.vectors[triangle_indices].reshape(-1, 3)
    mins = np.min(all_vertices, axis=0)
    maxs = np.max(all_vertices, axis=0)
    return mins, maxs

def build_bvh(stl_mesh, triangle_indices=None, depth=0, max_triangles_per_leaf=10):
    """
    Build a BVH for the triangles in the STL mesh.

    Args:
        stl_mesh: The `mesh.Mesh` object containing the STL triangles.
        triangle_indices: Indices of the triangles to consider (None for all).
        depth: Current depth of the BVH (used for recursion).
        max_triangles_per_leaf: Maximum number of triangles in a leaf node.

    Returns:
        Root node of the BVH.
    """
    if triangle_indices is None:
        triangle_indices = np.arange(len(stl_mesh.vectors))

    if len(triangle_indices) <= max_triangles_per_leaf:
        return BVHNode(bbox=compute_aabb(triangle_indices, stl_mesh), triangles=triangle_indices)

    # Compute AABB for current set of triangles
    bbox = compute_aabb(triangle_indices, stl_mesh)
    longest_axis = np.argmax(bbox[1] - bbox[0])

    # Sort triangle indices along the longest axis
    centroids = np.mean(stl_mesh.vectors[triangle_indices], axis=1)
    sorted_indices = triangle_indices[np.argsort(centroids[:, longest_axis])]

    # Split triangles into two groups
    mid = len(sorted_indices) // 2
    left = build_bvh(stl_mesh, sorted_indices[:mid], depth + 1, max_triangles_per_leaf)
    right = build_bvh(stl_mesh, sorted_indices[mid:], depth + 1, max_triangles_per_leaf)

    return BVHNode(bbox=bbox, left=left, right=right)


def ray_aabb_intersection(ray_origin, ray_direction, bbox):
    """
    Check if a ray intersects an axis-aligned bounding box (AABB).

    Args:
        ray_origin (np.array): Origin of the ray (3D).
        ray_direction (np.array): Direction of the ray (3D).
        bbox (tuple): Bounding box as (min_corner, max_corner), where
                      min_corner and max_corner are np.array of shape (3,).

    Returns:
        bool: True if the ray intersects the AABB, False otherwise.
    """
    min_corner, max_corner = bbox
    tmin = (min_corner - ray_origin) / ray_direction
    tmax = (max_corner - ray_origin) / ray_direction

    # Handle division by zero: if ray_direction is zero, set tmin and tmax to +inf or -inf
    tmin = np.where(ray_direction != 0, tmin, -np.inf)
    tmax = np.where(ray_direction != 0, tmax, np.inf)

    # Ensure tmin < tmax for each axis
    t1 = np.minimum(tmin, tmax)
    t2 = np.maximum(tmin, tmax)

    # Find the entry and exit points
    t_entry = np.max(t1)
    t_exit = np.min(t2)

    # Ray intersects if entry <= exit and exit > 0
    return t_entry <= t_exit and t_exit > 0




def traverse_bvh(ray_origin, ray_direction, node, stl_mesh, intersections):
    """
    Recursively traverse the BVH and find intersections with triangles in the STL mesh.
    """
    if node is None:
        return

    # Check ray-AABB intersection
    if not ray_aabb_intersection(ray_origin, ray_direction, node.bbox):
        return

    if node.triangles is not None:  # Leaf node
        for triangle_index in node.triangles:
            intersection = ray_triangle_intersection(ray_origin, ray_direction, stl_mesh, triangle_index)
            if intersection is not None:
                intersections.append(intersection)
    else:  # Internal node
        traverse_bvh(ray_origin, ray_direction, node.left, stl_mesh, intersections)
        traverse_bvh(ray_origin, ray_direction, node.right, stl_mesh, intersections)


def parallel_cast_rays_bvh(bvh, stl_mesh, grid_dims, voxel_size, axis=0, num_threads=8):
    """
    Parallelized raycasting using BVH and STL mesh.

    Args:
        bvh: Root node of the BVH.
        stl_mesh: The `mesh.Mesh` object containing the STL triangles.
        grid_dims: Dimensions of the voxel grid (nx, ny, nz).
        voxel_size: Size of each voxel.
        axis: Axis along which to cast rays (0 for x, 1 for y, 2 for z).
        num_threads: Number of threads for parallel processing.

    Returns:
        A dictionary where keys are (slice_idx, ray_idx) and values are sorted intersection points.
    """
    ray_direction = np.zeros(3)
    ray_direction[axis] = 1  # Direction of rays (+x, +y, or +z)
    orthogonal_axes = [i for i in range(3) if i != axis]

    def process_ray(slice_idx, ray_idx):
        ray_origin = np.zeros(3)
        ray_origin[orthogonal_axes[0]] = (slice_idx+0.5) * voxel_size[orthogonal_axes[0]]
        ray_origin[orthogonal_axes[1]] = (ray_idx +0.5)* voxel_size[orthogonal_axes[1]]
        ray_origin[axis]=-0.5*voxel_size[axis]
        intersections = []
        traverse_bvh(ray_origin, ray_direction, bvh, stl_mesh, intersections)
        # Sort intersections by distance from ray origin
        unique_intersections = deduplicate_intersections(intersections, ray_origin)
        return (slice_idx, ray_idx), unique_intersections

    tasks = []
    min_bound,max_bound = get_bounds(stl_mesh)
    start_voxel = (min_bound// voxel_size).astype(int)
    end_voxel =  (max_bound// voxel_size).astype(int)
    for slice_idx in range(start_voxel[orthogonal_axes[0]],end_voxel[orthogonal_axes[0]]+1):
        for ray_idx in range(start_voxel[orthogonal_axes[1]],end_voxel[orthogonal_axes[1]]+1): #for ray_idx in range(grid_dims[orthogonal_axes[1]]):
            tasks.append((slice_idx, ray_idx))

    # Parallel processing
    results = {}

    # Create a ThreadPoolExecutor for parallel processing
    with ThreadPoolExecutor(max_workers=num_threads) as executor:
        # Submit tasks to the executor and wrap them with tqdm for progress monitoring
        future_to_task = {executor.submit(process_ray, *task): task for task in tasks}

        # Use tqdm to monitor the completion of tasks
        for future in tqdm(as_completed(future_to_task), total=len(tasks), desc="Raytracing Progress", unit="ray"):
            task = future_to_task[future]
            try:
                (slice_idx, ray_idx), intersection_points = future.result()
                results[(slice_idx, ray_idx)] = intersection_points
            except Exception as e:
                print(f"Error processing ray {task}: {e}")

    return results


#def sort_intersections(intersections, ray_origin):
#    """
#    Sort intersections by distance from the ray origin.
#    """
#    if len(intersections) == 0:
#        return []
#
#    # Compute distances for each intersection point
#    distances = [np.linalg.norm(intersection - ray_origin) for intersection in intersections]
#
#    # Sort intersections by distance
#    sorted_intersections = [point for _, point in sorted(zip(distances, intersections))]
#    return sorted_intersections

def deduplicate_intersections(intersections, ray_origin, epsilon=1e-6):
    """
    Deduplicate intersections by their distance from the ray origin.

    Args:
        intersections (list of np.array): List of intersection points (3D NumPy arrays).
        ray_origin (np.array): The origin of the ray.
        ray_direction (np.array): The normalized direction of the ray.
        epsilon (float): Tolerance for considering two distances as identical.

    Returns:
        list: Deduplicated and sorted list of intersection points.
    """
    if len(intersections) == 0:
        return []

    # Calculate t values (distance along the ray direction)
    t_values =  [np.linalg.norm(intersection - ray_origin) for intersection in intersections]

    # Sort intersections by t-values and convert points to tuples for comparability
    sorted_intersections = sorted(zip(t_values, [tuple(point) for point in intersections]))

    # Deduplicate intersections based on t-values
    unique_intersections = []
    last_t = -np.inf

    for t, point_tuple in sorted_intersections:
        if abs(t - last_t) > epsilon:  # If t is different from the last unique t, add it
            unique_intersections.append(np.array(point_tuple))
            last_t = t

    return unique_intersections

def load_stl(filename):
    """_summary_

    Args:
        filename (_type_): _description_

    Returns:
        _type_: _description_
    """
    return mesh.Mesh.from_file(filename)

def load_and_translate_stl_to_positive_quadrant(stl_file_path,translation_vector=None):
    """
    Load an STL file and translate its mesh to the positive quadrant.

    Args:
        stl_file_path (str): Path to the STL file.

    Returns:
        stl.Mesh: Translated STL mesh object.
    """
    # Load the STL file
    stl_mesh = mesh.Mesh.from_file(stl_file_path)


    # Compute the bounding box (min and max coordinates of the mesh)
    min_coords = np.min(stl_mesh.vectors.reshape(-1, 3), axis=0)
    max_coords = np.max(stl_mesh.vectors.reshape(-1, 3), axis=0)

    print(f"Bounding Box Min: {min_coords}, Max: {max_coords}")

    # Calculate translation vector to shift mesh to the positive quadrant
    if translation_vector is None:
        translation_vector = -min_coords  # Shift everything by the negative of the minimum point

    print(f"Translation Vector: {translation_vector}")

    # Apply the translation to all vertices
    stl_mesh.vectors += translation_vector

    # Recalculate normals (optional but recommended after translation)
    stl_mesh.update_normals()

    return stl_mesh, translation_vector


def get_bounds(stl_mesh):
    """_summary_

    Args:
        stl_mesh (_type_): _description_

    Returns:
        _type_: _description_
    """
    min_bound = np.min(stl_mesh.points.reshape(-1, 3), axis=0)
    max_bound = np.max(stl_mesh.points.reshape(-1, 3), axis=0)
    return min_bound, max_bound

class VoxelGrid:
    def __init__(self, resolution, min_bound, max_bound):
        """
        Initialize the VoxelGrid object.

        Args:
            data (numpy.ndarray): 3D numpy array containing the voxel data.
            resolution (tuple of float): Size of each voxel in (x, y, z) directions.
        """
       
        self.dims = ((max_bound - min_bound) / resolution).astype(int)+1 # Dimensions of the grid (z, y, x)
        self.data = np.zeros(self.dims, dtype=bool)
        self.resolution = resolution  # Voxel size in physical units (e.g., mm)


def create_voxel_grid(stl_mesh, resolution,voxel_grid=None ):
    """_summary_

    Args:
        stl_mesh (_type_): _description_
        resolution (_type_): _description_

    Returns:
        voxel_grid _type_: _description_
    """

    if voxel_grid is None:
        min_bound, max_bound = get_bounds(stl_mesh)
        #grid_dims = ((max_bound - min_bound) / resolution).astype(int)+2
        #print(grid_dims)
        ## Create an empty voxel grid
        #voxel_grid = np.zeros(grid_dims, dtype=bool)
        voxel_grid = VoxelGrid( resolution, min_bound, max_bound)

    
    # Compute the voxel grid dimensions
    
    total_iterations=stl_mesh.vectors.shape[0]

    bvh = build_bvh(stl_mesh)

    # Perform parallel raycasting
    intersections = parallel_cast_rays_bvh(bvh,stl_mesh, voxel_grid.dims, resolution, axis=0, num_threads=8)

    # Mark internal voxels based on intersections
    mark_internal_voxels(voxel_grid.data, intersections, axis=0, voxel_size=resolution)



    ## Iterate through each triangle in the STL mesh and fill the surface voxels
    #with tqdm(total=total_iterations, desc="Checking surface intersection") as pbar:
    #    for triangle in stl_mesh.vectors:
    #        fill_triangle(voxel_grid, triangle, min_bound, resolution)
    #        pbar.update(1)



    return voxel_grid


def fill_triangle(voxels, triangle, min_bound, resolution):
    """_summary_

    Args:
        voxels (_type_): _description_
        triangle (_type_): _description_
        min_bound (_type_): _description_
        resolution (_type_): _description_
    """
    # Convert triangle vertices to voxel coordinates
    voxel_indices = ((triangle - min_bound) / resolution).astype(int)

    voxel_coords = ((triangle - min_bound) / resolution)
    # Get the bounding box of the triangle in voxel coordinates
    min_voxel = np.min(voxel_indices, axis=0)
    max_voxel = np.max(voxel_indices, axis=0)

    # Iterate through the bounding box
    for x in range(min_voxel[0], max_voxel[0] + 1):
        for y in range(min_voxel[1], max_voxel[1] + 1):
            for z in range(min_voxel[2], max_voxel[2] + 1):
                if voxel_intersects_triangle([x, y, z], voxel_coords):
                    voxels[x, y, z] = True

def check_ray_intersections(voxel_center, directions, triangle):
    intersections = 0
    for ray_direction in directions:
        if ray_intersects_triangle_woop(voxel_center, ray_direction, triangle):
            intersections += 1
    return intersections


def voxel_intersects_triangle_par(voxel, triangle):
    voxel_center = np.array(voxel) + 0.5  # Center of the voxel
    intersections = 0

    directions = [
        np.array([1, 0, 0]), np.array([0, 1, 0]), np.array([0, 0, 1]),
        np.array([-1, 0, 0]), np.array([0, -1, 0]), np.array([0, 0, -1]),
        np.array([1, 1, 1]), np.array([1, 1, -1]), np.array([1, -1, 1]),
        np.array([1, -1, -1]), np.array([-1, 1, 1]), np.array([-1, 1, -1]),
        np.array([-1, -1, 1]), np.array([-1, -1, -1])
    ]

    with ThreadPoolExecutor() as executor:
        future_intersections = [executor.submit(check_ray_intersections, voxel_center, [d], triangle) for d in directions]
        for future in future_intersections:
            intersections += future.result()

    # A voxel is considered on the surface if it intersects the triangle
    return intersections % 2 == 1

def voxel_intersects_triangle(voxel, triangle):
    """ The voxel_intersects_triangle function casts rays from the voxel center along the x, y, and z axes to detect intersections with the triangle. If the number of intersections is odd, the voxel is considered on the surface.

    Args:
        voxel (_type_): _description_
        triangle (_type_): _description_

    Returns:
        _type_: _description_
    """
    voxel_center = np.array(voxel) + 0.5  # Center of the voxel
    intersections = 0

    # Cast rays along the x, y, and z axes
    for axis in range(3):
        ray_direction = np.zeros(3)
        ray_direction[axis] = 1
        if ray_intersects_triangle_woop(voxel_center, ray_direction, triangle):
            intersections += 1
    # Cast rays along the -x, -y, and -z axes
    for axis in range(3):
        ray_direction = np.zeros(3)
        ray_direction[axis] = -1
        if ray_intersects_triangle_woop(voxel_center, ray_direction, triangle):
            intersections += 1
    
    # Cast rays along the main diagonals
    diagonals = [
        [1, 1, 1],
        [1, 1, -1],
        [1, -1, 1],
        [1, -1, -1],
        [-1, 1, 1],
        [-1, 1, -1],
        [-1, -1, 1],
        [-1, -1, -1]
    ]

    for diagonal in diagonals:
        ray_direction = np.array(diagonal)
        if ray_intersects_triangle_woop(voxel_center, ray_direction, triangle):
            intersections += 1
    

    # A voxel is considered on the surface if it intersects the triangle
    return intersections %2 == 1    

def ray_triangle_intersection(ray_origin, ray_direction, stl_mesh, triangle_index):
    """The ray_intersects_triangle function uses the Möller–Trumbore algorithm to determine if a ray intersects a triangle.

    Args:
        ray_origin (_type_): _description_
        ray_direction (_type_): _description_
        stl_mesh (_type_):
        triangle_index (int or int array): index of tringles

    Returns:
        intersection_point: intersection coordinate
    """
    # Möller–Trumbore intersection algorithm
    epsilon = 1e-12
    triangle = stl_mesh.vectors[triangle_index]
    v0, v1, v2 = triangle  # Unpack the vertices
    edge1 = v1 - v0
    edge2 = v2 - v0
    h = np.cross(ray_direction, edge2)
    a = np.dot(edge1, h)
    if -epsilon < a < epsilon:
        return None  # Ray is parallel to the triangle
    f = 1.0 / a
    s = ray_origin - v0
    u = f * np.dot(s, h)
    if u < 0.0 or u > 1.0:
        return None
    q = np.cross(s, edge1)
    v = f * np.dot(ray_direction, q)
    if v < 0.0 or u + v > 1.0:
        return None
    t = f * np.dot(edge2, q)
    if t > epsilon:  # Intersection exists
        intersection_point = ray_origin + t * ray_direction
        return intersection_point
    return None  # No valid intersection


def ray_intersects_triangle_woop(ray_origin, ray_direction, triangle):
    """Ray-Triangle Intersection Woop algorithm

    Args:
        ray_origin (_type_): _description_
        ray_direction (_type_): _description_
        triangle (_type_): _description_

    Returns:
        _type_: _description_
    """    
    # Transform the ray into the triangle's coordinate system

    v0, v1, v2 = triangle[0], triangle[1], triangle[2]
    edge1 = v1 - v0
    edge2 = v2 - v0

    # Calculate normal and the distance from the plane of the triangle
    normal = np.cross(edge1, edge2)
    d = np.dot(normal, v0)
    
    # Transform the ray origin and direction into the triangle's coordinate system
    ray_origin_transformed = ray_origin - v0
    ray_direction_transformed = np.dot(ray_direction, normal)
    
    if np.abs(ray_direction_transformed) < 1e-8:
        return False  # Ray is parallel to the triangle plane

    t = (d - np.dot(normal, ray_origin)) / ray_direction_transformed

    if t < 0:
        return False  # Triangle is behind the ray

    # Calculate the intersection point
    P = ray_origin + t * ray_direction

    # Check if the intersection point is inside the triangle using barycentric coordinates
    u = np.cross(P - v0, edge2)
    v = np.cross(edge1, P - v0)

    if (np.dot(normal, u) < 0) or (np.dot(normal, v) < 0):
        return False

    # Calculate barycentric coordinates
    denom = np.dot(normal, normal)
    b1 = np.dot(normal, u) / denom
    b2 = np.dot(normal, v) / denom
    b0 = 1 - b1 - b2

    return (b0 >= 0) and (b1 >= 0) and (b2 >= 0)

def cast_rays(triangles, grid_dims, voxel_size, axis=0):
    """
    Cast rays along the given axis and collect intersection points with triangles.

    Args:
        triangles: List of triangles (each triangle is a list of 3 vertices).
        grid_dims: Dimensions of the voxel grid (nx, ny, nz).
        voxel_size: Size of each voxel.
        axis: Axis along which to cast rays (0 for x, 1 for y, 2 for z).

    Returns:
        A dictionary where keys are (slice_idx, ray_idx) and values are sorted intersection points.
    """
    intersections = {}
    ray_direction = np.zeros(3)
    ray_direction[axis] = 1  # Direction of the rays (e.g., +x)
    orthogonal_axes = [i for i in range(3) if i != axis]
    
    for slice_idx in range(grid_dims[orthogonal_axes[0]]):  # For each slice
        for ray_idx in range(grid_dims[orthogonal_axes[1]]):  # For each ray in the slice
            # Define the ray origin based on the current slice and ray index
            ray_origin = np.zeros(3)
            ray_origin[orthogonal_axes[0]] = slice_idx * voxel_size
            ray_origin[orthogonal_axes[1]] = ray_idx * voxel_size

            # Collect intersections
            ray_intersections = []
            for triangle in triangles:
                intersection = ray_triangle_intersection(ray_origin, ray_direction, triangle)
                if intersection is not None:
                    ray_intersections.append(intersection[axis])  # Store the coordinate along the ray's axis
            
            # Sort intersections along the ray
            ray_intersections.sort()
            intersections[(slice_idx, ray_idx)] = ray_intersections
    return intersections

def mark_internal_voxels(voxel_grid, intersections, axis, voxel_size):
    """
    Mark internal voxels based on ray-triangle intersections.

    Args:
        voxel_grid: 3D numpy array representing the voxel grid.
        intersections: Dictionary of ray intersections (from cast_rays).
        axis: Axis along which the rays were cast.
        voxel_size: Size of each voxel.
    """
    orthogonal_axes = [i for i in range(3) if i != axis]

    for (slice_idx, ray_idx), points in intersections.items():
        inside = True
        for i in range(1, len(points)):
            # Toggle inside state on every pair of intersections
            start_voxel = (points[i - 1] // voxel_size).astype(int)
            if np.any(start_voxel < 0):
                # Raise a warning
                warnings.warn("The array contains negative values, which will be set to 0.", UserWarning)
                # Set negative values to 0
                start_voxel[start_voxel < 0] = 0
            end_voxel = (points[i] // voxel_size).astype(int)
            if inside:
                # Mark voxels between start and end as internal
                for voxel in range(start_voxel[axis], end_voxel[axis] + 1):
                    idx = [0, 0, 0]
                    idx[axis] = voxel
                    idx[orthogonal_axes[0]] = slice_idx
                    idx[orthogonal_axes[1]] = ray_idx
                    #print(idx)
                    voxel_grid[tuple(idx)] = True
                    
            inside = not inside

def ray_traced_voxelization(stl_triangles, grid_dims, voxel_size, axis=0):
    """
    Perform ray-traced voxelization to identify internal voxels.

    Args:
        stl_triangles: List of triangles from the STL surface.
        grid_dims: Dimensions of the voxel grid (nx, ny, nz).
        voxel_size: Size of each voxel.
        axis: Axis along which rays will be cast.

    Returns:
        3D numpy array of voxels (True for internal voxels).
    """
    voxel_grid = np.zeros(grid_dims, dtype=bool)

    # Step 1: Cast rays and collect intersections
    intersections = cast_rays(stl_triangles, grid_dims, voxel_size, axis)

    # Step 2: Mark internal voxels based on intersections
    mark_internal_voxels(voxel_grid, intersections, axis, voxel_size)

    return voxel_grid



def save_voxel_grid_as_raw(voxel_grid, filename):
    """Save the voxel grid as a raw binary file

    Args:
        voxel_grid numpy.array() :The numpy array voxel grid
        filename String path: output filename
    """
    with open(filename, 'wb') as f:
        voxel_grid.astype(np.uint8).tofile(f)


def save_voxel_data_as_raw_mhd(voxel_data, voxel_size, raw_filename, mhd_filename):
    """
    Save voxel data to a RAW file with accompanying MHD metadata file.

    Args:
        voxel_data (numpy.ndarray): 3D array containing the voxel data (e.g., 1 for inside, 0 for outside).
        voxel_size (tuple of float): Size of each voxel in (x, y, z) directions, typically in millimeters.
        raw_filename (str): Path to save the RAW file.
        mhd_filename (str): Path to save the MHD file.
    """
    # Ensure the directory exists
    #os.makedirs(os.path.dirname(raw_filename), exist_ok=True)

    # Write the voxel data to the RAW file
    voxel_data.astype(np.int32).tofile(raw_filename)

    # Prepare metadata for the MHD file
    mhd_content = f"""ObjectType = Image
                      NDims = 3
                      DimSize = {voxel_data.shape[2]} {voxel_data.shape[1]} {voxel_data.shape[0]}
                      ElementType = MET_INT
                      ElementSpacing = {voxel_size[2]} {voxel_size[1]} {voxel_size[0]}
                      Offset = 0 0 0
                      ElementDataFile = {os.path.basename(raw_filename)}
                      BinaryData = True
                      BinaryDataByteOrderMSB = False
                      CompressedData = False
    """

    # Write the metadata to the MHD file
    with open(mhd_filename, 'w') as mhd_file:
        mhd_file.write(mhd_content)

    print(f"Voxel data saved to {raw_filename} and {mhd_filename}")

def display_stl_with_matplotlib(stl_mesh):
    """
    Display an STL file using Matplotlib.

    Args:
        stl_file (str): Path to the STL file.
    """


    # Extract the vertices and faces
    vertices = stl_mesh.vectors.reshape(-1, 3)
    faces = np.arange(vertices.shape[0]).reshape(-1, 3)

    # Plot the mesh using Matplotlib
    fig = plt.figure(figsize=(10, 10))
    ax = fig.add_subplot(111, projection='3d')

    # Plot the surface triangles
    ax.plot_trisurf(vertices[:, 0], vertices[:, 1], faces, vertices[:, 2], linewidth=0.1, cmap='gray', edgecolor='k')

    # Set equal aspect ratio for all axes
    ax.set_box_aspect([1, 1, 1])  # Equal aspect ratio for x, y, z

    # Set axis limits to ensure equal scaling
    max_range = (vertices.max(axis=0) - vertices.min(axis=0)).max() / 2
    mid_x = (vertices[:, 0].max() + vertices[:, 0].min()) / 2
    mid_y = (vertices[:, 1].max() + vertices[:, 1].min()) / 2
    mid_z = (vertices[:, 2].max() + vertices[:, 2].min()) / 2

    ax.set_xlim(mid_x - max_range, mid_x + max_range)
    ax.set_ylim(mid_y - max_range, mid_y + max_range)
    ax.set_zlim(mid_z - max_range, mid_z + max_range)

    # Set the axes limits and labels
    ax.set_xlabel("X")
    ax.set_ylabel("Y")
    ax.set_zlabel("Z")
    ax.set_title("STL Mesh Visualization")

    # Show the plot
    plt.show()



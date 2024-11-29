# GAMS: A GPU-accelerated MRI solver for synthetic image generation

A pipeline for the generation of synthetic MRI images, including the design and acquisition of the MRI pulse sequence, the definition of the 3D computational domain (digital phantom), its discretization in voxels, the time integration of the Bloch equations in such domain, and the numerical signal processing. 
The MRI simulator includes three components: 

1) the pre-processor, which enables us to input the geometrical, physical, excitation and acquisition parameters;
2) the solver which is used to solve the Bloch equation and enables us to produce the raw-data signal by means of in-house developed CUDA-Fortran routines;
3) the post-processor, which enables us to process the generated signal and convert it into a DICOM image, with proper metadata.


## Tools
Contains some scripts to write some MRI test sequences, Generate B0 distortion field components, create voxelized geometries.

## TEST
Contains a bunch of generated sequences and the corresponding DICOM files.

## SLZ_to_hdf5_MRI

C++ program to convert the optional output solution in binary file to a xdmf + h5 file format that can e opened with paraview.
Never managed to compile with the FORTRAN H5 libraries.

Maybe one day this will be python too.


### To do list.

Clean up of missing routines.
Upload of missing routines.
Streamline of the pipeline
Friendly GUI
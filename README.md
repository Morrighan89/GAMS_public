# GAMS: A GPU-accelerated MRI solver for synthetic image generation

A pipeline for the generation of synthetic MRI images, including the design and acquisition of the MRI pulse sequence, the definition of the 3D computational domain (digital phantom), its discretization in voxels, the time integration of the Bloch equations in such domain, and the numerical signal processing. 
The MRI simulator includes three components: 

1) the pre-processor, which enables us to input the geometrical, physical, excitation and acquisition parameters;
2) the solver which is used to solve the Bloch equation and enables us to produce the raw-data signal by means of in-house developed CUDA-Fortran routines; (to be released)
3) the post-processor, which enables us to process the generated signal and convert it into a DICOM image, with proper metadata. 


## Tools

Contains some scripts to write some MRI test sequences, Generate B0 distortion field components, create voxelized geometries.

## Test Object

Digital phantom for MRI simulation. Test object with vials.

## TEST

Contains a bunch of generated sequences and the corresponding DICOM files.

## SLZ_to_hdf5_MRI

C++ program to convert the optional output solution in a binary file to an xdmf + h5 file format that can be opened with Paraview.
Never managed to compile with the FORTRAN H5 libraries.

Future revisions will be in python too.

## Acknowledgment

This research was developed in the framework of Project 20NRM05 IMET-MRI, which has received funding from the European Metrology Programme for Innovation and Research (EMPIR), co-financed by the participating states, and from the European Union’s Horizon 2020 Programme.

### To do list.

Clean up of missing routines.
Upload of missing routines.
Streamline of the pipeline
Friendly GUI

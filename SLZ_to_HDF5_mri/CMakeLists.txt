# CMakeList.txt : CMake project for SLZ_to_HDF5_mri, include source and define
# project specific logic here.
#
cmake_minimum_required (VERSION 3.10.1)


# Add source to this project's executable.
add_executable (SLZ_to_HDF5_mri "SLZ_to_HDF5_mri.cpp" "SLZ_to_HDF5_mri.h" "Mesh.hpp" "slz.hpp" "xdmf_writer.hpp" "h5_writer.hpp" "precision.hpp")

include_directories(${HDF5_INCLUDE_DIR})
link_directories( ${HDF5_LIBRARY_DIRS} )

target_link_libraries(${PROJECT_NAME}
	
	"C:/Program Files/HDF_Group/HDF5/1.12.0/lib/libzlib.lib"
	"C:/Program Files/HDF_Group/HDF5/1.12.0/lib/libszip.lib"
	${HDF5_C_LIBRARIES}
)
# TODO: Add tests and install targets if needed.

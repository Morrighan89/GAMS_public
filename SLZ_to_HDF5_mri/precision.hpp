#ifndef precision_H
#define precision_H
#include <hdf5.h>
//typedef float real_t;
//hid_t H5_real_type = H5Tcopy(H5T_NATIVE_FLOAT);
typedef double real_t;
hid_t H5_real_type = H5Tcopy(H5T_NATIVE_DOUBLE);
/*Before compiling uncomment the desired precision*/ 

#endif

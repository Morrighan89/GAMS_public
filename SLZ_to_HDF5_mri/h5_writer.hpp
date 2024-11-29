#ifndef writeH5_HEADER
#define writeH5_HEADER

#include <iostream>
#include <fstream>
#include <math.h>
#include <string>
#include <hdf5.h>
#include "Mesh.hpp"
#include "slz.hpp"
#include "precision.hpp"
//this is a class to write the H5 file complementing the hdf5 file obtained from our simulation
namespace slzToHdf5
{
	class writeH5{
	public:
		void postProcess(real_t& time);
		writeH5(std::string& Name, MeshData* ptrMesh, SlzData* ptrSlz);


		MeshData* ptrMesh;
		SlzData* ptrSlz;
		MatData* ptrMat;
		void writeInitH5();
		void writeToH5();



		void writeTopology(hid_t file_id);
		void writeGeometry();
		void writeSolution();
		void writeHext();
		void writeMean();
		void writeTimeStamp(const void* ptrData, int numTime);
		void writeMatProp(MatData* ptrMat);
		void writeScalarDatastructure(hid_t file_id, int num, std::string datasetname, const void* ptrData, int var_type=0);
		void writeVectorDatastructure(hid_t file_id, int num, std::string datasetname, const void* ptrData);
		void writeVariableDatastructure(hid_t file_id, int num1, int num2, std::string Name, const void* ptrData);
		void writeStructuredVectorDatastructure(hid_t file_id, int nx, int ny, int nz, std::string Name, const void* ptrData);
		void writeStructuredScalarDatastructure(hid_t file_id, int nx, int ny, int nz, std::string datasetName, const void* ptrData, int var_type=0);
		void writeVariable();
		void writeScalar();
		void writeVector();



		std::ofstream     M_H5;
		hid_t     file_id;
		const std::string M_closingLines;
		std::streampos    M_closingLinesPosition;
		std::string       M_outputFileName;
		std::string       M_CaseName;

	};
	//Constructors
	writeH5::writeH5(std::string& Name, MeshData* Mesh, SlzData* Slz) :

		M_CaseName(Name),
		ptrMesh(Mesh),
		ptrSlz(Slz)
		
		{
		M_outputFileName = this->M_CaseName + ".H5";
	}

	//Methods
	void writeH5::postProcess(real_t& time)
	{

		M_outputFileName = this->M_CaseName + ".H5";


		// write empty xdmf file
		writeInitH5();
	}


	
	//Protected Methods

	void writeH5::writeInitH5()
	{

		const char *fileName = M_outputFileName.c_str();
		file_id = H5Fcreate(fileName, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
		int TimeSteps[1];
		TimeSteps[0] = ptrSlz->M_numSavedInstant;
		writeScalarDatastructure(file_id, 1, "TimestepsNumber",TimeSteps, 1);
		real_t volElem[1];
		volElem[0] = ptrMesh->M_volElem;
		writeScalarDatastructure(file_id, 1, "volElem", volElem, 0);
	}


	void writeH5::writeToH5()
	{
		const char *fileName = M_outputFileName.c_str();
		file_id = H5Fopen(fileName, H5F_ACC_RDWR, H5P_DEFAULT);
		if (file_id == -1){ std::cout << "error opening h5"; }

	}

	


	void writeH5::writeTopology(hid_t file_id)
	{
		hid_t     dataset_id, dataspace_id, group_id;
		hsize_t   dims[2];
		herr_t    status;
		int arr[3];
		dims[0] = 1;
		dims[1] = 3;

		arr[0] = ptrMesh->meshSizeX()+1; //xmf Paraview takes as input the number of nodes not cells per side
		arr[1] = ptrMesh->meshSizeY()+1;
		arr[2] = ptrMesh->meshSizeZ()+1;

		std::cout << "creazione dataspace contenete edges hdf5" << std::endl;
		dataspace_id = H5Screate_simple(2, dims, NULL);

		std::cout << "Inizio scrittura meshsize in hdf5" << std::endl;
		dataset_id = H5Dcreate2(file_id, "NXNYNZ", H5T_NATIVE_INT, dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

		status = H5Dwrite(dataset_id, H5T_NATIVE_INT, H5S_ALL, dataspace_id, H5P_DEFAULT, arr);
		std::cout << "scritto" << std::endl;
		status = H5Dclose(dataset_id);
		std::cout << "chiuso dataset" << std::endl;
		status = H5Sclose(dataspace_id);
		std::cout << "chiuso dataspace" << std::endl;

		
	}

	void writeH5::writeGeometry()
	{
		// Write the data file.
		hid_t     dataset_id, dataspace_id, group_id;
		hsize_t   dims[2];
		herr_t    status;
		real_t arr[3];
		arr[0] = ptrMesh->meshDX();
		arr[1] = ptrMesh->meshDY();
		arr[2] = ptrMesh->meshDZ();
		writeVectorDatastructure(file_id, 1, "/DXDYDZ", &arr);
		writeStructuredScalarDatastructure(file_id, (*ptrMesh).meshSizeX(), (*ptrMesh).meshSizeY(), (*ptrMesh).meshSizeZ(), "/MatCodes", &(ptrMesh->M_VoxelCode[0][0][0]),1);
		//writeScalarDatastructure(file_id, (*ptrMesh).numGlobalVertices(), "/PointsX", ptrMesh->M_PointX);
		//writeScalarDatastructure(file_id, (*ptrMesh).numGlobalVertices(), "/PointsY", ptrMesh->M_PointY);
		//writeScalarDatastructure(file_id, (*ptrMesh).numGlobalVertices(), "/PointsZ", ptrMesh->M_PointZ);
		//writeScalarDatastructure(file_id, (*ptrMesh).numGlobalElements(), "/BaricX", ptrMesh->M_BarX);
		//writeScalarDatastructure(file_id, (*ptrMesh).numGlobalElements(), "/BaricY", ptrMesh->M_BarY);
		//writeScalarDatastructure(file_id, (*ptrMesh).numGlobalElements(), "/BaricZ", ptrMesh->M_BarZ);
		status = H5Fclose(file_id);

	}
	void writeH5::writeSolution()
	{
		// Write the data file.
		hid_t     dataset_id, dataspace_id, group_id;
		hsize_t   dims[2];
		herr_t    status;

		writeStructuredVectorDatastructure(file_id, (*ptrMesh).meshSizeX(), (*ptrMesh).meshSizeY(),(*ptrMesh).meshSizeZ(), "/MM", &(ptrSlz->M_M[0][0][0][0]));
		writeStructuredVectorDatastructure(file_id, (*ptrMesh).meshSizeX(), (*ptrMesh).meshSizeY(), (*ptrMesh).meshSizeZ(), "/B_Field", &(ptrSlz->M_Hext[0][0][0][0]));
		//writeVectorDatastructure(file_id, (*ptrMesh).numGlobalElements(), "/Emme", &(ptrSlz->M_M[0][0]));
		//writeVectorDatastructure(file_id, (*ptrMesh).numGlobalElements(), "/Heff", &(ptrSlz->M_Heff[0][0]));
		//writeVectorDatastructure(file_id, (*ptrMesh).numGlobalElements(), "/Hms", &(ptrSlz->M_Hms[0][0]));
		//writeVectorDatastructure(file_id, (*ptrMesh).numGlobalElements(), "/Hexch", &(ptrSlz->M_Hexch[0][0]));
		//writeHext();
		status = H5Fclose(file_id);

	}

	void writeH5::writeHext()
	{
		// Write the data file.
		hid_t     dataset_id, dataspace_id, group_id;
		hsize_t   dims[2];
		herr_t    status;
		real_t* M_Buf_Hext;
		real_t** M_Hext;

		//M_Buf_Hext = new real_t[ptrSlz->M_numGlobalElements * 3];
		//M_Hext = new real_t*[ptrSlz->M_numGlobalElements];
		//for (int i = 0; i < ptrSlz->M_numGlobalElements; i++) M_Hext[i] = M_Buf_Hext + i * 3;
		//for (int i = 0; i < ptrSlz->M_numGlobalElements; i++){
		//	M_Hext[i][0] = ptrSlz->M_Hext[0][0];
		//	M_Hext[i][1] = ptrSlz->M_Hext[0][1];
		//	M_Hext[i][2] = ptrSlz->M_Hext[0][2];
		//}
		//writeVectorDatastructure(file_id, (*ptrMesh).numGlobalElements(), "/Hext", &(M_Hext[0][0]));

	}
	void writeH5::writeMean()
	{
		// Write the data file.
		herr_t    status;
		writeVariableDatastructure(file_id,1, 3, "/MeanM", &(ptrSlz->meanEmme[0]));
		status = H5Fclose(file_id);

	}
	void writeH5::writeTimeStamp(const void* ptrData, int numTime)
	{
		// Write the data file.
		herr_t    status;
		writeScalarDatastructure(file_id, numTime,  "/Timestamps", ptrData, 0);
		status = H5Fclose(file_id);

	}

	void writeH5::writeMatProp(MatData* ptrMat)
	{
		// Write the data file.
		herr_t    status;
		writeScalarDatastructure(file_id, (*ptrMesh).numGlobalElements(), "/Saturation", ptrMat->M_Elem_SatMag, 0);
		writeScalarDatastructure(file_id, (*ptrMesh).numGlobalElements(), "/ExchLength", ptrMat->M_Elem_AlEx, 0);
		writeScalarDatastructure(file_id, (*ptrMesh).numGlobalElements(), "/Alfa", ptrMat->M_Elem_Alfa, 0);
		writeScalarDatastructure(file_id, (*ptrMesh).numGlobalElements(), "/Material", ptrMat->M_Elem_MatType,1);
		writeScalarDatastructure(file_id, (*ptrMesh).numGlobalElements(), "/Anisotropy", ptrMat->M_Elem_AniType,1);
		status = H5Fclose(file_id);

	}
	void writeH5::writeScalarDatastructure(hid_t file_id, int num, std::string datasetName, const void* ptrData, int var_type)
	{
		// Write the data file.
		hid_t     dataset_id, dataspace_id, group_id;
		hsize_t   dims[2];
		herr_t    status;
		const char *dataName = datasetName.c_str();
		dims[0] = num;
		dims[1] = 1;
		dataspace_id = H5Screate_simple(2, dims, NULL);

		if (var_type == 0){
			dataset_id = H5Dcreate(file_id, dataName,  H5_real_type, dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

			status = H5Dwrite(dataset_id,  H5_real_type, H5S_ALL, H5S_ALL, H5P_DEFAULT, ptrData);
		}
		else{
			dataset_id = H5Dcreate(file_id, dataName, H5T_NATIVE_INT, dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

			status = H5Dwrite(dataset_id, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, ptrData);
		}

		

		status = H5Dclose(dataset_id);

		status = H5Sclose(dataspace_id);


	}
	void writeH5::writeVectorDatastructure(hid_t file_id,int num, std::string Name, const void* ptrData)
	{
		// Write the data file.
		hid_t     dataset_id, dataspace_id, group_id;
		hsize_t   dims[2];
		herr_t    status;
		Name  += std::to_string(ptrSlz->M_numTimeInstant);
		const char *groupName = Name.c_str();
		group_id = H5Gcreate(file_id, groupName, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
		dims[0] = num;
		dims[1] = 3;
		dataspace_id = H5Screate_simple(2, dims, NULL);
		Name += "/Val";
		const char *datasetName = Name.c_str();
		dataset_id = H5Dcreate(file_id, datasetName,  H5_real_type, dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

		status = H5Dwrite(dataset_id,  H5_real_type, H5S_ALL, H5S_ALL, H5P_DEFAULT, ptrData);

		status = H5Dclose(dataset_id);

		status = H5Sclose(dataspace_id);


	}

	void writeH5::writeStructuredVectorDatastructure(hid_t file_id, int nx, int ny, int nz, std::string Name, const void* ptrData)
	{
		// Write the data file.
		hid_t     dataset_id, dataspace_id, group_id;
		hsize_t   dims[4];
		herr_t    status;
		Name += std::to_string(ptrSlz->M_numTimeInstant);
		const char* groupName = Name.c_str();
		group_id = H5Gcreate(file_id, groupName, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
		dims[0] = nx;
		dims[1] = ny;
		dims[2] = nz;
		dims[3] = 3;
		dataspace_id = H5Screate_simple(4, dims, NULL);
		Name += "/Val";
		const char* datasetName = Name.c_str();
		dataset_id = H5Dcreate(file_id, datasetName, H5_real_type, dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

		status = H5Dwrite(dataset_id, H5_real_type, H5S_ALL, H5S_ALL, H5P_DEFAULT, ptrData);

		status = H5Dclose(dataset_id);

		status = H5Sclose(dataspace_id);


	}

	void writeH5::writeStructuredScalarDatastructure(hid_t file_id, int nx, int ny, int nz, std::string datasetName, const void* ptrData, int var_type)
	{
		// Write the data file.
		hid_t     dataset_id, dataspace_id, group_id;
		hsize_t   dims[4];
		herr_t    status;
		const char* dataName = datasetName.c_str();
		dims[0] = nx;
		dims[1] = ny;
		dims[2] = nz;
		dims[3] = 1;
		dataspace_id = H5Screate_simple(4, dims, NULL);

		if (var_type == 0) {
			dataset_id = H5Dcreate(file_id, dataName, H5_real_type, dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

			status = H5Dwrite(dataset_id, H5_real_type, H5S_ALL, H5S_ALL, H5P_DEFAULT, ptrData);
		}
		else {
			dataset_id = H5Dcreate(file_id, dataName, H5T_NATIVE_INT, dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

			status = H5Dwrite(dataset_id, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, ptrData);
		}



		status = H5Dclose(dataset_id);

		status = H5Sclose(dataspace_id);


	}

	void writeH5::writeVariableDatastructure(hid_t file_id, int num1, int num2, std::string Name, const void* ptrData)
	{
		// Write the data file.
		hid_t     dataset_id, dataspace_id, group_id;
		hsize_t   dims[2];
		herr_t    status;
		Name += std::to_string(ptrSlz->M_numTimeInstant);
		const char *groupName = Name.c_str();
		group_id = H5Gcreate(file_id, groupName, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
		dims[0] = num1;
		dims[1] = num2;
		dataspace_id = H5Screate_simple(2, dims, NULL);
		Name += "/Val";
		const char *datasetName = Name.c_str();
		dataset_id = H5Dcreate(file_id, datasetName, H5_real_type, dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

		status = H5Dwrite(dataset_id,  H5_real_type, H5S_ALL, H5S_ALL, H5P_DEFAULT, ptrData);

		status = H5Dclose(dataset_id);

		status = H5Sclose(dataspace_id);


	}

}
#endif
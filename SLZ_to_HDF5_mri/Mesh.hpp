#ifndef MeshData_HEADER
#define MeshData_HEADER


#include <iostream>
#include <fstream>
#include <math.h>
#include <string>
#include "Mesh.hpp"
#include "precision.hpp"

//this is a class to write the xdmf file complementing the hdf5 file obtained from our simulation
namespace slzToHdf5
{
	class MeshData
	{
	public:

		int  numGlobalElements(){ return M_numGlobalElements; }
		int  numGlobalVertices(){ return M_numGlobalVertices; }
		int	 numLocalVertices(){ return M_numLocalVertices; }
		int  meshTotalSize() { return M_intNx * M_intNy * M_intNz; }
		int  meshSizeX() { return M_intNx; }
		int  meshSizeY() { return M_intNy; }
		int  meshSizeZ() { return M_intNz; }
		real_t  meshDX() { return M_dblDx; }
		real_t  meshDY() { return M_dblDy; }
		real_t  meshDZ() { return M_dblDz; }
		int  meshType() { return M_meshType; }
		MeshData(std::string& caso, int& Dim, int& EdgeNumber);
		void readMeshStructure1();
		void readMeshStructure2();
		void readMeshStructure3();
		void readMeshData2();
		void readMeshPoint();
		void readMeshEdge();
		void readMeshVoxelCode();
		void readMeshBar();
		void readMaterial();
		void readVolElem();
		void initializeMeshEdge();
		void initializeMeshPoint();
		void initializeMeshBar();
		void initializeVolElem();
		void initializeMaterial();
		void initializeVoxelCode();


		real_t* M_PointX;
		real_t* M_PointY;
		real_t* M_PointZ;

		real_t* M_BarX;
		real_t* M_BarY;
		real_t* M_BarZ;

		real_t M_volElem;

		int* M_Buf_Edge;
		int** M_Edge;

		int* M_BufBuf_VoxelCode;
		int** M_Buf_VoxelCode;
		int*** M_VoxelCode;
	
		int* M_Material;
		
		
	protected:

		int  M_numGlobalElements;
		int  M_numElements;
		int  M_numLocalVertices;
		int  M_numGlobalVertices;
		std::string       M_CaseName;

		real_t M_dblXCellaMin;
		real_t M_dblXCellaMax;
			   
		real_t M_dblYCellaMin;
		real_t M_dblYCellaMax;

		real_t M_dblZCellaMin;
		real_t M_dblZCellaMax;
			   
		real_t M_dblXDominioMin;
		real_t M_dblXDominioMax;
			   
		real_t M_dblYDominioMin;
		real_t M_dblYDominioMax;

		real_t M_dblZDominioMin;
		real_t M_dblZDominioMax;
			   
		real_t M_dblSpessore;
			   
		real_t M_dblToll;

		int M_intNumCelle;

		int M_intNx;
		int M_intNy;
		int M_intNz;

		real_t M_dblDx;
		real_t M_dblDy;
		real_t M_dblDz;


	
	private:
		
		std::ifstream   M_Geo;
		std::string     M_meshDir;
		std::string     M_meshFile;
		std::string     M_matFile;
		int             M_meshType;
		int             M_Dim;
		int             M_EdgeNumber;
		std::string     M_order;
		std::streampos  M_EdgeStreamPos;
		std::streampos  M_NodeStreamPos;
		std::streampos  M_VoxelCodeStreamPos;

	};

	//Constructors
	MeshData::MeshData(std::string& caso, int& Dim, int& EdgeNumber):
		M_meshFile(caso+".GEO"),
		//M_matFile(caso + ".MAT"),
		M_Dim(Dim),
		M_EdgeNumber(EdgeNumber),
		M_meshType(Dim*10+EdgeNumber)
		
	{	
		M_Geo.open(M_meshFile, std::ios::in | std::ios::binary | std::ios::ate);
	}

	//Methods

	void MeshData::readMeshStructure1(){

		int posizione;
		posizione = 0;

		if (!M_Geo.is_open()){ M_Geo.open(M_meshFile, std::ios::in | std::ios::binary | std::ios::ate); }

		if (M_Geo.is_open()){

			M_Geo.seekg(0, M_Geo.end);
			
			M_Geo.seekg(4, M_Geo.beg);

			std::cout << M_Geo.is_open() << " M_GEO aperto \n";
			M_Geo.read((char *)& M_dblXCellaMin, sizeof(real_t));
			M_Geo.read((char *)& M_dblXCellaMax, sizeof(real_t));

			std::cout << posizione << "\n";
			M_Geo.seekg(8, M_Geo.cur);
			M_Geo.read((char *)& M_dblYCellaMin, sizeof(real_t));
			M_Geo.read((char *)& M_dblYCellaMax, sizeof(real_t));
			M_Geo.seekg(8, M_Geo.cur);
			M_Geo.read((char *)& M_dblXDominioMin, sizeof(real_t));
			M_Geo.read((char *)& M_dblXDominioMax, sizeof(real_t));
			M_Geo.seekg(8, M_Geo.cur);
			M_Geo.read((char *)& M_dblYDominioMin, sizeof(real_t));
			M_Geo.read((char *)& M_dblYDominioMax, sizeof(real_t));
			M_Geo.seekg(8, M_Geo.cur);
			M_Geo.read((char *)& M_dblSpessore, sizeof(real_t));
			M_Geo.seekg(8, M_Geo.cur);
			M_Geo.read((char *)& M_dblToll, sizeof(real_t));
			M_Geo.seekg(8, M_Geo.cur);
			M_Geo.read((char *)& M_intNumCelle, sizeof(int));
			M_Geo.seekg(8, M_Geo.cur);
			M_Geo.read((char *)& M_numGlobalVertices, sizeof(int));
			M_Geo.seekg(8, M_Geo.cur);
			M_EdgeStreamPos = M_Geo.tellg();
		}
		
	}

	void MeshData::readMeshStructure2(){

		int posizione;
		posizione = 0;

		if (!M_Geo.is_open()){ M_Geo.open(M_meshFile, std::ios::in | std::ios::binary | std::ios::ate); }

		if (M_Geo.is_open()){

			M_Geo.seekg(0, M_Geo.end);

			M_Geo.seekg(0, M_Geo.beg);

			std::cout << M_Geo.is_open() << " M_GEO aperto \n";
			M_Geo.read((char *)& M_intNx, sizeof(int));
			M_Geo.read((char *)& M_intNy, sizeof(int));
			M_Geo.read((char *)& M_intNz, sizeof(int));
			
			M_Geo.read((char *)& M_dblDx, sizeof(real_t));
			M_Geo.read((char *)& M_dblDy, sizeof(real_t));
			M_Geo.read((char *)& M_dblDz, sizeof(real_t));
	
			M_Geo.read((char *)& M_intNumCelle, sizeof(int));
			M_Geo.read((char *)& M_numGlobalVertices, sizeof(int));
			M_Geo.read((char *)& M_numElements, sizeof(int));
			M_Geo.read((char *)& M_numGlobalElements, sizeof(int));
			M_Geo.read((char *)& M_volElem, sizeof(real_t));

			M_Geo.read((char *)& M_numLocalVertices, sizeof(int));

			M_EdgeStreamPos = M_Geo.tellg();
		}

	}
	void MeshData::readMeshStructure3() {

		int posizione;
		posizione = 0;

		if (!M_Geo.is_open()) { M_Geo.open(M_meshFile, std::ios::in | std::ios::binary | std::ios::ate); }

		if (M_Geo.is_open()) {

			M_Geo.seekg(0, M_Geo.end);

			M_Geo.seekg(0, M_Geo.beg);

			std::cout << M_Geo.is_open() << " M_GEO aperto \n";
			M_Geo.read((char*)&M_intNx, sizeof(int));
			M_Geo.read((char*)&M_intNy, sizeof(int));
			M_Geo.read((char*)&M_intNz, sizeof(int));

			M_Geo.read((char*)&M_dblDx, sizeof(real_t));
			M_Geo.read((char*)&M_dblDy, sizeof(real_t));
			M_Geo.read((char*)&M_dblDz, sizeof(real_t));

			M_volElem = M_dblDx * M_dblDy * M_dblDz;
			//M_Geo.read((char*)&M_intNumCelle, sizeof(int));
			//M_Geo.read((char*)&M_numGlobalVertices, sizeof(int));
			//M_Geo.read((char*)&M_numElements, sizeof(int));
			//M_Geo.read((char*)&M_numGlobalElements, sizeof(int));f
			//M_Geo.read((char*)&M_volElem, sizeof(real_t));
			//
			//M_Geo.read((char*)&M_numLocalVertices, sizeof(int));
			//
			M_VoxelCodeStreamPos = M_Geo.tellg();
		}
		else { std::cout << "Could not open: " << M_meshFile << "\n"; }

	}
	void MeshData::initializeMeshPoint(){
		//M_PointX = (real_t *)malloc(numGlobalVertices() * sizeof(real_t));
		//M_PointY = (real_t *)malloc(numGlobalVertices() * sizeof(real_t));
		//M_PointZ = (real_t *)malloc(numGlobalVertices() * sizeof(real_t));
		M_PointX= new real_t[M_numGlobalVertices];
		M_PointY= new real_t[M_numGlobalVertices];
		M_PointZ= new real_t[M_numGlobalVertices];

	}
	void MeshData::initializeMaterial(){
		M_Material = new int[M_intNx * M_intNy * M_intNz];

	}
	void MeshData::initializeVoxelCode() {
		M_BufBuf_VoxelCode = new int[M_intNx * M_intNy * M_intNz];
		M_Buf_VoxelCode = new int* [M_intNx * M_intNy];
		M_VoxelCode = new int** [M_intNx];
		for (int k = 0; k < M_intNx * M_intNy; k++) M_Buf_VoxelCode[k] = M_BufBuf_VoxelCode + k * M_intNz;
		for (int j = 0; j < M_intNx; j++) M_VoxelCode[j] = M_Buf_VoxelCode + j* M_intNy;


	}
	

	/*void MeshData::initializeVolElem(){
		real_t *M_VolElem = (real_t *)malloc(M_numGlobalElements * sizeof(real_t));

	}*/

	void MeshData::initializeMeshEdge(){
		M_Buf_Edge = new int[M_numGlobalElements*M_numLocalVertices];
		M_Edge = new int*[M_numGlobalElements];
		for (int i = 0; i < M_numGlobalElements; i++) M_Edge[i] = M_Buf_Edge + i * M_numLocalVertices;

	}
	void MeshData::initializeMeshBar(){
		M_BarX = new real_t[M_numGlobalElements];
		M_BarY = new real_t[M_numGlobalElements];
		M_BarZ = new real_t[M_numGlobalElements];
	}
	void MeshData::readMeshPoint(){
		int Vertices;
		initializeMeshPoint();
		M_Geo.seekg(M_NodeStreamPos, M_Geo.beg);
		M_Geo.read((char *)& Vertices, sizeof(int));
		if (Vertices == M_numGlobalVertices){
			for (int i = 0; i < M_numGlobalVertices; i++) {
				M_Geo.read((char *)& M_PointX[i], sizeof(real_t));
				M_Geo.read((char *)& M_PointY[i], sizeof(real_t));
				M_Geo.read((char *)& M_PointZ[i], sizeof(real_t));
				//std::cout << M_PointX[i] << " " << M_PointY[i] << " " << M_PointZ[i] << " " << "\n";
			}
		}
		else{ std::cout << "wrong number of vertex" << "\n"; }
	}

	/*void MeshData::readMeshBar(){
		M_Geo.read((char *)& M_BarX[i], sizeof(real_t));
		M_Geo.read((char *)& M_BarY[i], sizeof(real_t));
		M_Geo.read((char *)& M_BarZ[i], sizeof(real_t));
	}*/

	void MeshData::readMeshEdge(){

		initializeMaterial();
		initializeMeshEdge();
		initializeMeshBar();
		M_Geo.seekg(M_EdgeStreamPos, M_Geo.beg);
		

		for (int ii = 0; ii < M_numGlobalElements; ii++) {

			for (int jj = 0; jj < M_numLocalVertices; jj++) {
				M_Geo.read((char *)& M_Edge[ii][jj], sizeof(int));
				M_NodeStreamPos = M_Geo.tellg();
			}

			M_Geo.read((char *)& M_Material[ii], sizeof(int));

			M_Geo.read((char *)& M_BarX[ii], sizeof(real_t));
			M_Geo.read((char *)& M_BarY[ii], sizeof(real_t));
			M_Geo.read((char *)& M_BarZ[ii], sizeof(real_t));

			M_NodeStreamPos = M_Geo.tellg();
		}
	}
	void MeshData::readMeshVoxelCode() {

		initializeVoxelCode();

		M_Geo.seekg(M_VoxelCodeStreamPos, M_Geo.beg);


		for (int kk = 0; kk < M_intNz; kk++) {
			//std::cout << "ii: " << ii << "\n";
			for (int jj = 0; jj < M_intNy; jj++) {
				//std::cout << "jj: " << jj << "\n";
				for (int ii = 0; ii < M_intNx; ii++) {
//					M_BufBuf_VoxelCode[ii * M_intNy * M_intNz + jj * M_intNz + kk]= ii * M_intNy * M_intNz + jj * M_intNz + kk;
					M_Geo.read((char*)&M_VoxelCode[ii][jj][kk], sizeof(int));
					//std::cout << "Voxel value: M_VoxelCode" << M_VoxelCode[ii][jj][kk] << ii << jj << kk <<"\n";
//					std::cout << "Voxel value: M_Buf_VoxelCode" << M_Buf_VoxelCode[ii* M_intNy+jj][kk] << "\n";
					M_VoxelCodeStreamPos = M_Geo.tellg();
				}
			}

			M_VoxelCodeStreamPos = M_Geo.tellg();
		}
	}
	/*void MeshData::readMeshMaterial(){
			M_Geo.read((char *)& M_Material[i], sizeof(int));
	}*/
	/*void MeshData::readMeshData2(){
		for (int i = 0; i < M_numGlobalElements; i++) {
			MeshData::readMeshEdge(i);
			MeshData::readMeshMat(i);
			fileGeo.seekg(8, fileGeo.cur);
			MeshData::readMeshBaric(i);
			fileGeo.seekg(8, fileGeo.cur);
		}
		for (int i = 0; i < M_numGlobalVertices; i++) {
			MeshData::readMeshPoint(i);
		}


	}*/
	class MatData{
	public:
		
		MatData(std::string& caso, MeshData* Mesh);
		void readMatData();
		void readMaterial();
		void initializeMatProp();
		void generateElemProp();

		int* M_Buf_Edge;
		int** M_Edge;

		int* M_Material;
		double M_GammaG;
		int M_matMax;
		double* M_matVol;
		double* M_SatMag;
		double* M_Alfa;
		double* M_AlEx;
		int* M_AniType;
		double* M_Anisotropy;
		int* M_MatType;
		real_t * M_Elem_SatMag;
		real_t * M_Elem_Alfa;
		real_t * M_Elem_AlEx;
		int* M_Elem_AniType;
		int* M_Elem_MatType;


	protected:
		MeshData* ptrMesh;
		std::string       M_CaseName;

	private:

		std::ifstream   M_Mat;
		std::string     M_meshDir;
		std::string     M_meshFile;
		std::string     M_matFile;
		std::streampos  M_EdgeStreamPos;
		std::streampos  M_NodeStreamPos;

	};
	//Constructors
	MatData::MatData(std::string& caso, MeshData* Mesh) :
		M_matFile(caso + ".MAT"),
		ptrMesh(Mesh)

	{
		M_Mat.open(M_matFile, std::ios::in);
	}
	//Methods
	void MatData::readMatData(){

		int posizione;
		std::string line;
		posizione = 0;

		if (!M_Mat.is_open()){ M_Mat.open(M_matFile, std::ios::in); }

		if (M_Mat.is_open()){
			getline(M_Mat, line);
			M_GammaG = std::stold(line);
			getline(M_Mat, line);
			M_matMax = std::stoi(line);
			//std::string orbits("365.24 29.53");
			std::string::size_type sz;     // alias of size_t
			//
			//real_t earth = std::stod(orbits, &sz);
			//real_t moon = std::stod(orbits.substr(sz));
			M_matMax = M_matMax + 1;
			M_matVol= new double[M_matMax];
			M_SatMag = new double[M_matMax];
			M_Alfa = new double[M_matMax];
			M_AniType = new int[M_matMax];
			M_Anisotropy = new double[M_matMax];
			M_AlEx = new double[M_matMax];
			M_MatType= new int[M_matMax];
			int i = 0;
			while (i<M_matMax)
			{
				getline(M_Mat, line);
				M_MatType[i] = std::stoi(line, &sz);
				M_matVol[i] = std::stod(line.substr(sz));
				getline(M_Mat, line);
				M_SatMag[i] = std::stod(line, &sz);
				std::cout << M_SatMag[i];
				M_Alfa[i] = std::stod(line.substr(sz));
				getline(M_Mat, line);
				M_AlEx[i] = std::stod(line);
				getline(M_Mat, line);
				M_AniType[i] = std::stoi(line);
				if (M_AniType[i]==0){
				}
				else if (M_AniType[i] == 1){
					getline(M_Mat, line);
				}
				else if (M_AniType[i] == 3){
					getline(M_Mat, line);
					getline(M_Mat, line);
					getline(M_Mat, line);
				}
				i++;
			}
			M_Mat.close();
		}

	}
	void  MatData::generateElemProp(){
		M_Elem_SatMag = new real_t[(*ptrMesh).numGlobalElements()];
		M_Elem_Alfa = new real_t[(*ptrMesh).numGlobalElements()];
		M_Elem_AlEx = new real_t[(*ptrMesh).numGlobalElements()];
		M_Elem_AniType = new int[(*ptrMesh).numGlobalElements()];
		M_Elem_MatType = new int[(*ptrMesh).numGlobalElements()];
		for (int i = 0; i < (*ptrMesh).numGlobalElements(); i++){
			M_Elem_SatMag[i] = M_SatMag[ptrMesh->M_Material[i]];
			M_Elem_Alfa[i] = M_Alfa[ptrMesh->M_Material[i]];
			M_Elem_AlEx[i] = M_AlEx[ptrMesh->M_Material[i]];
			M_Elem_AniType[i] = M_AniType[ptrMesh->M_Material[i]];
			M_Elem_MatType[i] = ptrMesh->M_Material[i];
		}

	}
};

#endif
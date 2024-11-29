#ifndef SlzData_HEADER
#define SlzData_HEADER

#include <iostream>
#include <fstream>
#include <math.h>
#include <string>
#include "Mesh.hpp"
#include "precision.hpp"
//this is a class to write the xdmf file complementing the hdf5 file obtained from our simulation
namespace slzToHdf5
{
	class SlzData{
	public:
		void postProcess(real_t& time);
		int  M_numGlobalElements;
		int  M_numGlobalVertices;
		int	 M_numLocalVertices;
		int  M_numTimeInstant;
		int  M_numSavedInstant;
		int  M_intNx;
		int  M_intNy;
		int  M_intNz;
		int numSavedInstant(){ return M_numSavedInstant; }
		real_t M_Time;
		int M_ChangeField;

		SlzData(std::string& caso, MeshData* Mesh);
		real_t* M_BufBufBuf_M;
		real_t** M_BufBuf_M;
		real_t*** M_Buf_M;
		real_t**** M_M;
		real_t* M_BufBufBuf_Hext;
		real_t** M_BufBuf_Hext;
		real_t*** M_Buf_Hext;
		real_t**** M_Hext;
		real_t* M_Buf_Hms;
		real_t** M_Hms;
		real_t* M_Buf_Hexch;
		real_t** M_Hexch;
		real_t* M_Buf_Heff;
		real_t** M_Heff;
		void initializeSlzVectorAttribute(real_t* dblPtr, real_t** dblPtrPtr, real_t*** dblPtrPtrPtr, real_t**** dblPtrPtrPtrPtr);
		void initializeSlzScalarAttribute(real_t* dblPtr, real_t** dblPtrPtr, real_t*** dblPtrPtrPtr);
		void readSlzVectorAttribute(real_t**** dblPtrPtrPtrPtr, int i, int j, int k);
		void readSlzData();
		void readSlz();
		void InitializeData();
		void writeSlzData();
		void getTimeStep();
		real_t meanEmme[3];
		real_t meanHeff[3];
		real_t meanHexch[3];
		real_t meanHms[3];
		
		void computeMean(real_t** Qty, real_t* Mean);
		void computeEnergy(real_t** Qty, real_t* Mean);

		std::streampos  M_StreamPos;

	protected:
		MeshData* M_mesh;
		
		std::string       M_CaseName;

	private:

		std::ifstream   M_slz;
		std::string     M_slzDir;
		std::string     M_slzFile;


	};

	//Constructors
	SlzData::SlzData(std::string& caso, MeshData* Mesh) :
		M_slzFile(caso + ".SLZ"),
		M_mesh(Mesh)
	{	
		const char *dataName = M_slzFile.c_str();
		M_slz.open(dataName, std::ios::in | std::ios::binary );
		if (!M_slz.is_open()){
			std::cout << "could not open file" << std::endl;
			//size_t errmsglen = strerrorlen_s(errno) + 1;
			char errmsg[200];
			strerror_s(errmsg, 200, errno);
			std::cerr << "Error: " << errmsg;
		}
		M_slz.seekg(0, M_slz.end);
		M_StreamPos = M_slz.tellg();
		std::cout << M_StreamPos << "\n";
		M_slz.seekg(0, M_slz.beg);
		M_StreamPos = M_slz.tellg();
		std::cout << M_StreamPos << "\n";
		M_numGlobalElements=(*M_mesh).numGlobalElements();
		M_numGlobalVertices = (*M_mesh).numGlobalVertices();
		M_intNx = M_mesh->meshSizeX();
		M_intNy = M_mesh->meshSizeY();
		M_intNz = M_mesh->meshSizeZ();

		M_numTimeInstant = 0;

	}

	//Methods
	void SlzData::initializeSlzVectorAttribute(real_t* dblPtr, real_t** dblPtrPtr, real_t*** dblPtrPtrPtr, real_t**** dblPtrPtrPtrPtr){

		for (int i = 0; i < M_intNx*M_intNy*M_intNz; i++) dblPtrPtr[i] = dblPtr + i * (int)3;
		for (int i = 0; i < M_intNx * M_intNy; i++) dblPtrPtrPtr[i] = dblPtrPtr + i * M_intNz;
		for (int i = 0; i < M_intNx; i++) dblPtrPtrPtrPtr[i] = dblPtrPtrPtr + i * M_intNy;

	}
	void SlzData::initializeSlzScalarAttribute(real_t* dblPtr, real_t** dblPtrPtr, real_t*** dblPtrPtrPtr){
		dblPtr = new real_t  [M_intNx * M_intNy * M_intNz];
		dblPtrPtr = new real_t * [M_intNx * M_intNy];
		dblPtrPtrPtr = new real_t **[M_intNx];
		for (int i = 0; i < M_intNx * M_intNy; i++) dblPtrPtr[i] = dblPtr + i * M_intNz;
		for (int i = 0; i < M_intNx; i++) dblPtrPtrPtr[i] = dblPtrPtr + i * M_intNy;

	}
// ATTENTION Coordinates are inverted beacuse paraview exchanges X and Z in the mesh s
	void SlzData::readSlzVectorAttribute(real_t**** dblPtrPtrPtrPtr, int i, int j, int k){
		M_slz.seekg(M_StreamPos, M_slz.beg);
		M_slz.read((char *)& dblPtrPtrPtrPtr[i][j][k][2], sizeof(real_t));
		M_slz.read((char *)& dblPtrPtrPtrPtr[i][j][k][1], sizeof(real_t));
		M_slz.read((char *)& dblPtrPtrPtrPtr[i][j][k][0], sizeof(real_t));
		//std::cout << dblPtrPtr[i][0] << " " << dblPtrPtr[i][1] << " " << dblPtrPtr[i][2] << std::endl;

		
		M_StreamPos = M_slz.tellg();
	}
	void SlzData::readSlzData(){
		if (!M_slz.is_open()){
			M_slz.open(M_slzFile, std::iostream::in | std::iostream::binary | std::iostream::ate);
		if (!M_slz.is_open()){
			std::cout << "could not open file" << std::endl;
			//size_t errmsglen = strerrorlen_s(errno) + 1;
			char errmsg[200];
			strerror_s(errmsg, 200, errno);
			std::cerr << "Error: " << errmsg;
		}
		
		}

		if (M_slz.is_open()){
			M_slz.seekg(0, M_slz.cur);
			for (int k = 0; k < M_intNz; k++) {
				for (int j = 0; j < M_intNy; j++) {
					for (int i = 0; i < M_intNx; i++) {
						M_StreamPos = M_slz.tellg();
						readSlzVectorAttribute(M_M, i, j, k);
						readSlzVectorAttribute(M_Hext, i, j, k);
						//std::cout << M_M[k][j][i][0] << std::endl;
						//std::cout << M_M[k][j][i][1] << std::endl;
						//std::cout << M_M[k][j][i][2] << std::endl;

						
					}
				}				
				//readSlzVectorAttribute(M_Buf_Emme, M_Emme, i);
				//readSlzVectorAttribute(M_Buf_Heff, M_Heff, i);
				//readSlzVectorAttribute(M_Buf_Hms, M_Hms, i);
				//readSlzVectorAttribute(M_Buf_Hexch, M_Hexch, i);
				//int j;
				//M_slz.read((char *)& j, sizeof(int));
				//std::cout << j << std::endl;
				//M_slz.seekg(8, M_slz.cur);
				
			}
		}
	}
	void SlzData::readSlz(){
		if (!M_slz.is_open()){ M_slz.open(M_slzFile, std::iostream::in | std::iostream::binary | std::iostream::ate); }
		if (M_slz.is_open()){
			M_slz.read((char *)& M_numTimeInstant, sizeof(int));
			if (M_numTimeInstant != -1){
				//M_slz.read((char *)& M_Time, sizeof(real_t));
				//M_slz.read((char *)& M_ChangeField, sizeof(int));
				//std::cout << M_Time << std::endl;
				//M_slz.seekg(8, M_slz.cur);
				M_slz.read((char *)& M_Hext[0][0][0][0], sizeof(real_t));
				M_slz.read((char *)& M_Hext[0][0][0][1], sizeof(real_t));
				M_slz.read((char *)& M_Hext[0][0][0][2], sizeof(real_t));
				M_StreamPos = M_slz.tellg();
				//M_slz.seekg(8, M_slz.cur);
				//std::cout << M_Hext[0][0][0][0] << std::endl;
				//std::cout << M_Hext[0][0][0][1] << std::endl;
				//std::cout << M_Hext[0][0][0][2] << std::endl;
				readSlzData();
			}
			
		}

	}
	void SlzData::InitializeData(){
		M_BufBufBuf_M = new real_t[M_intNx * M_intNy * M_intNz * 3];
		M_BufBuf_M = new real_t * [M_intNx * M_intNy * M_intNz];
		M_Buf_M = new real_t * *[M_intNx * M_intNy];
		M_M = new real_t * **[M_intNx];
		initializeSlzVectorAttribute(M_BufBufBuf_M, M_BufBuf_M, M_Buf_M, M_M);
		M_BufBufBuf_Hext = new real_t[M_intNx * M_intNy * M_intNz * 3];
		M_BufBuf_Hext = new real_t * [M_intNx * M_intNy * M_intNz];
		M_Buf_Hext = new real_t * *[M_intNx * M_intNy];
		M_Hext = new real_t * **[M_intNx];
		initializeSlzVectorAttribute(M_BufBufBuf_Hext, M_BufBuf_Hext, M_Buf_Hext, M_Hext);
		//M_Buf_Heff = new real_t[M_numGlobalElements * 3];
		//M_Heff = new real_t*[M_numGlobalElements];
		//initializeSlzVectorAttribute(M_Buf_Heff, M_Heff);
		//M_Buf_Hms = new real_t[M_numGlobalElements * 3];
		//M_Hms = new real_t*[M_numGlobalElements];
		//initializeSlzVectorAttribute(M_Buf_Hms, M_Hms);
		//M_Buf_Hexch = new real_t[M_numGlobalElements * 3];
		//M_Hexch = new real_t*[M_numGlobalElements];
		//initializeSlzVectorAttribute(M_Buf_Hexch, M_Hexch);

	}
	void SlzData::getTimeStep(){
		int lastInstpos;
		int numTime;
		std::streampos size = 0;
		//lastInstpos = M_numGlobalElements*(16 + 4 * 3 * sizeof(real_t)) + 40;
		lastInstpos = M_numGlobalElements*(4*sizeof(int) + 4 * 3 * sizeof(real_t)) + 2*sizeof(int)+4*sizeof(real_t);
		M_slz.seekg(0, M_slz.end);
		size = M_slz.tellg();
		M_numSavedInstant =  size/ lastInstpos;
		
		M_slz.seekg(M_StreamPos, M_slz.beg);
	}
	void SlzData::computeMean(real_t** Qty, real_t* Mean){
		for (int i = 0; i < M_numGlobalElements; i++){
			Mean[0] = Mean[0] + Qty[i][0];
			Mean[1] = Mean[1] + Qty[i][1];
			Mean[2] = Mean[2] + Qty[i][2];
			
		}
		Mean[0] = Mean[0]/M_numGlobalElements;
		Mean[1] = Mean[1]/M_numGlobalElements;
		Mean[2] = Mean[2]/M_numGlobalElements;



	}
	void SlzData::computeEnergy(real_t** Qty, real_t* Mean) {
		for (int i = 0; i < M_numGlobalElements; i++) {
			Mean[0] = Mean[0] + Qty[i][0];
			Mean[1] = Mean[1] + Qty[i][1];
			Mean[2] = Mean[2] + Qty[i][2];

		}
		Mean[0] = Mean[0] / M_numGlobalElements;
		Mean[1] = Mean[1] / M_numGlobalElements;
		Mean[2] = Mean[2] / M_numGlobalElements;



	}
};
#endif
// SLZ_to_HDF5_mri.cpp : Defines the entry point for the application.
//

#include "SLZ_to_HDF5_mri.h"
#include <iostream>
#include <fstream>
#include <math.h>
#include <hdf5.h>
#include <string>
#include <vector>
#include <iterator>
#include <list>
#include "Mesh.hpp"
#include "slz.hpp"
#include "xdmf_writer.hpp"
#include "h5_writer.hpp"
#include "precision.hpp"


static void show_usage(std::string name)
{
	std::cerr << "Usage: " << name[0] << " <option(s)> SOURCES"
		<< "Options:\n"
		<< "\t-h,--help\t\tShow this help message\n"
		<< "\t-f,--filename input filename"
		<< "\t-t,--timesteps  number of time steps"
		<< "\t-mean, --mean requires argument(1 for computing the mean, 0 otherwise)"
		<< "\t-mat, --material requires argument(1 to save material properties, 0 otherwise)"
		<< "\t -equi, --equilibrium requires argument (1 for saving only equilibrium points, 0 otherwise)"
		<< std::endl;
}

int main(int argc, char* argv[])
{
	std::string fileName;
	int numTimestep;
	int computeMean;
	int saveMat;
	int saveEqui;
	real_t* dblTimeStamp;
	std::list <real_t> dblListTime;
	computeMean = 0;
	saveMat = 0;
	saveEqui = 1;
	numTimestep = 0;
	if (argc < 2) {
		show_usage(argv[0]);
		std::cout << "An input file is required \n Insert file name (without extension) \n";
		std::cin >> fileName;
	}

	for (int iarg = 1; iarg < argc; ++iarg) {
		std::string arg = argv[iarg];
		if ((arg == "-h") || (arg == "--help")) {
			show_usage(argv[0]);
			return 1;
		}
		else if ((arg == "-f") || (arg == "--filename")) {
			if (iarg + 1 < argc) { // Make sure we aren't at the end of argv!
				fileName = argv[++iarg]; // Increment 'i' so we don't get the argument as the next argv[i].
			}
			else { // Uh-oh, there was no argument to the destination option.
				std::cerr << "--filename option requires one argument." << std::endl;
				return 1;
			}
		}
		else if ((arg == "-t") || (arg == "--timesteps")) {
			if (iarg + 1 < argc) { // Make sure we aren't at the end of argv!
				numTimestep = atoi(argv[++iarg]); // Increment 'i' so we don't get the argument as the next argv[i].
			}
			else { // Uh-oh, there was no argument to the destination option.
				std::cerr << "--timesteps requires argument" << std::endl;
				return 1;
			}
		}
		else if ((arg == "-mean") || (arg == "--mean")) {
			if (iarg + 1 < argc) { // Make sure we aren't at the end of argv!
				computeMean = atoi(argv[++iarg]); // Increment 'i' so we don't get the argument as the next argv[i].
			}
			else { // Uh-oh, there was no argument to the destination option.
				std::cerr << "--mean requires argument (1 for computing the mean, 0 otherwise)" << std::endl;
				return 1;
			}
		}
		else if ((arg == "-mat") || (arg == "--material")) {
			if (iarg + 1 < argc) { // Make sure we aren't at the end of argv!
				saveMat = atoi(argv[++iarg]); // Increment 'i' so we don't get the argument as the next argv[i].
			}
			else { // Uh-oh, there was no argument to the destination option.
				std::cerr << "--material requires argument (1 for saving materials, 0 otherwise)" << std::endl;
				return 1;
			}
		}
		else if ((arg == "-equi") || (arg == "--equilibrium")) {
			if (iarg + 1 < argc) { // Make sure we aren't at the end of argv!
				saveEqui = atoi(argv[++iarg]); // Increment 'i' so we don't get the argument as the next argv[i].
			}
			else { // Uh-oh, there was no argument to the destination option.
				std::cerr << "--equi requires argument (1 for saving only equilibrium points, 0 otherwise)" << std::endl;
				return 1;
			}
		}

	}

	int i;
	i = 3;
	int j;
	j = 9;
	int MeshType = i * 10 + j;

	slzToHdf5::MeshData* Mesh = new slzToHdf5::MeshData(fileName, i, j);
	//p->DoSomething();
	//slzToHdf5::MeshData Mesh(pluto, i,j);
	Mesh->readMeshStructure3();
	//Mesh->readMeshEdge();
	Mesh->readMeshVoxelCode();
	slzToHdf5::MatData* Mat = new slzToHdf5::MatData(fileName, Mesh);
	//Mat->readMatData();
	slzToHdf5::SlzData* Slz = new slzToHdf5::SlzData(fileName, Mesh);
	int plutarco = Slz->M_numGlobalElements;
	if (numTimestep == 0) {
		//Slz->getTimeStep();
		//numTimestep = Slz->M_numSavedInstant;
		numTimestep = 0;
		std::cout << "Please insert the number of timesteps to convert to .H5 and press enter\n";
		std::cin >> numTimestep;
	}
	std::cout << Slz->M_numSavedInstant << " " << numTimestep << " Number of saved instant \n";
	std::cout << Mesh->meshDX() << " " << Mesh->meshDY() << " " << Mesh->meshDZ() << " " << " DX DY DZ \n";
	std::cout << Mesh->meshSizeX() << " " << Mesh->meshSizeY() << " " << Mesh->meshSizeZ() << " " << " NX NY NZ \n";
	slzToHdf5::writeXdmf Xmf(fileName, Mesh, Slz, computeMean, saveMat);
	real_t time = 0;
	Xmf.postProcess(time);

	slzToHdf5::writeH5* H5 = new slzToHdf5::writeH5(fileName, Mesh, Slz);
	H5->writeInitH5();
	H5->writeTopology(H5->file_id);
	H5->writeGeometry();
	Slz->InitializeData();
	if (saveMat == 1) {
		Mat->generateElemProp();
		H5->writeToH5();
		H5->writeMatProp(Mat);
	}
	int kk = 0;
	int kkk = 0;
	while (kk < numTimestep) {
		//std::cout << kk << std::endl;
		Slz->readSlz();
		//if (computeMean == 1) {
		//	Slz->computeMean(Slz->M_Emme, Slz->meanEmme);
		//	H5->writeToH5();
		//	H5->writeMean();
		//}
		//if (saveEqui == 1) {
		//	if (Slz->M_ChangeField == 1) {
		//		H5->writeToH5();
		//		H5->writeSolution();
		//		Xmf.writeToXdmf(kkk, MeshType);
		//		dblListTime.push_back(Slz->M_Time);
		//		kkk++;
		//	}
		//}
		//else {
			dblListTime.push_back(Slz->M_Time);
			H5->writeToH5();
			H5->writeSolution();
			Xmf.writeToXdmf(Slz->M_numTimeInstant, MeshType);
		//}
		kk++;

	}
	dblTimeStamp = new real_t[dblListTime.size()];
	int k = 0;
	for (real_t const& i : dblListTime) {
		dblTimeStamp[k++] = i;
	}
	H5->writeToH5();
	H5->writeTimeStamp(dblTimeStamp, dblListTime.size());
	std::cout << dblTimeStamp[k - 1] << "magnetization evolution  time \n";


	return 0;

}
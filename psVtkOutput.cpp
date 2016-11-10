#include <iostream>
#include <string>
#include <sstream>
#include "vtkVersion.h"
#include "vtkDoubleArray.h"
#include "vtkFloatArray.h"
#include "vtkPointData.h"
#include "vtkPoints.h"
#include "vtkXMLStructuredGridWriter.h"
#include "vtkXMLPStructuredGridWriter.h"
#include "vtkStructuredGrid.h"
#include "vtkSmartPointer.h"
#include "psVtkOutput.h"

void cpphello_()
{
	std::cout<<"hello world from cpp \n";
}

void dwritevtsc_(int* nPartition, int* nSize, int* nFlds, double* dPnts, double* dRFlds,double* dCFlds, char* chFileName, int nLenFN)
{
	int nTotal=(nSize[0]*nSize[1]*nSize[2]);
	//Assign filename to string, trim white spaces and then add .vts
	chFileName[nLenFN--]='\0'; //null terminate string
        std::string stFileName,stPFileName;
	std::stringstream ss;
	for(int i=0;i<nLenFN;i++)
	{
		if(chFileName[i]!=' ')
			stFileName+=chFileName[i];
		else
			break;
	}
	stPFileName=stFileName;
	ss<<*nPartition;
	stFileName+="_";
	stFileName+=ss.str();
	stFileName+=".vts";
	//std::cout<<"nTotal: "<<nTotal<<", "<<*nFlds<<", "<<*nPartition<<" "+stFileName<<"\n";
	//std::cout<<"nSize: "<<nSize[0]<<", "<<nSize[1]<<", "<<nSize[2]<<" "+stFileName<<"\n";

	//Step 1: setup grid, points and variable objects
	vtkSmartPointer<vtkStructuredGrid> sGrid =
		vtkSmartPointer<vtkStructuredGrid>::New();
	vtkSmartPointer<vtkPoints> points =
		vtkSmartPointer<vtkPoints>::New();
	vtkSmartPointer<vtkDoubleArray> R_velocity =
		vtkSmartPointer<vtkDoubleArray>::New();
	vtkSmartPointer<vtkDoubleArray> R_pressure =
		vtkSmartPointer<vtkDoubleArray>::New();
	vtkSmartPointer<vtkDoubleArray> R_temperature =
		vtkSmartPointer<vtkDoubleArray>::New();
	vtkSmartPointer<vtkDoubleArray> C_velocity =
		vtkSmartPointer<vtkDoubleArray>::New();
	vtkSmartPointer<vtkDoubleArray> C_pressure =
		vtkSmartPointer<vtkDoubleArray>::New();
	vtkSmartPointer<vtkDoubleArray> C_temperature =
		vtkSmartPointer<vtkDoubleArray>::New();

	//Step 2: Set Grid dimensions and variable tuples
	sGrid->SetDimensions(nSize[0],nSize[1],nSize[2]);
	points->Allocate(nTotal);
	R_velocity->SetNumberOfComponents(3);
	R_pressure->SetNumberOfComponents(1);
	R_temperature->SetNumberOfComponents(1);
	R_velocity->SetNumberOfTuples(nTotal);
	R_pressure->SetNumberOfTuples(nTotal);
	R_temperature->SetNumberOfTuples(nTotal);
	R_velocity->SetName("R_velocity");
	R_pressure->SetName("R_pressure");
	R_temperature->SetName("R_temperature");
	C_velocity->SetNumberOfComponents(3);
	C_pressure->SetNumberOfComponents(1);
	C_temperature->SetNumberOfComponents(1);
	C_velocity->SetNumberOfTuples(nTotal);
	C_pressure->SetNumberOfTuples(nTotal);
	C_temperature->SetNumberOfTuples(nTotal);
	C_velocity->SetName("C_velocity");
	C_pressure->SetName("C_pressure");
	C_temperature->SetName("C_temperature");
	
	//Step 3: Copy data into vtk Objects
	int j=0;
	for( int i=0;i<nTotal*3;i+=3)
	{
	//	std::cout<<"i="<<i<<"\n";
		points->InsertPoint(j,&dPnts[i]);
		j++;
	}

	//std::cout<<"Points allocated succesfully\n";
	j=0;
	for(int i=0;i<nTotal*(*nFlds);i+=(*nFlds))
	{
	//	std::cout<<"i="<<i<<"\n";
		R_velocity->InsertTuple(j,&dRFlds[i]);
		R_pressure->InsertTuple(j,&dRFlds[i+3]);
		R_temperature->InsertTuple(j,&dRFlds[i+4]);
		C_velocity->InsertTuple(j,&dCFlds[i]);
		C_pressure->InsertTuple(j,&dCFlds[i+3]);
		C_temperature->InsertTuple(j,&dCFlds[i+4]);
		j++;
	}
	//std::cout<<"Fields allocated succesfully\n";

	//Step 4: Assign values to grid
	sGrid->SetPoints(points);
	sGrid->GetPointData()->AddArray(C_velocity);
	sGrid->GetPointData()->AddArray(C_pressure);
	sGrid->GetPointData()->AddArray(C_temperature);
	sGrid->GetPointData()->SetVectors(R_velocity);
	sGrid->GetPointData()->AddArray(R_pressure);
	sGrid->GetPointData()->SetScalars(R_temperature);
	//std::cout<<"Grid allocated succesfully\n";

	//Step 5: Setup writer object and write file
		vtkSmartPointer<vtkXMLStructuredGridWriter> writer =
			vtkSmartPointer<vtkXMLStructuredGridWriter>::New();
		writer->SetFileName(&stFileName[0]);
		writer->SetInputData(sGrid);
		writer->Write();
	//Step 6: If partition zero write pvts vile
	//if(*nPartition==0)
	/*{
		stPFileName+=".pvts";
		vtkSmartPointer<vtkXMLPStructuredGridWriter> pwriter =
			vtkSmartPointer<vtkXMLPStructuredGridWriter>::New();
		pwriter->SetFileName(&stPFileName[0]);
		pwriter->SetInputData(sGrid);
		pwriter->SetNumberOfPieces(8);
		//pwriter->SetUpdateExtent(ext);
		pwriter->Write();
	}*/
}
void dwritevts_(int* nPartition, int* nSize, int* nFlds, double* dPnts, double* dFlds, char* chFileName, int nLenFN)
{
	int nTotal=(nSize[0]*nSize[1]*nSize[2]);
	//Assign filename to string, trim white spaces and then add .vts
	chFileName[nLenFN--]='\0'; //null terminate string
        std::string stFileName,stPFileName;
	std::stringstream ss;
	for(int i=0;i<nLenFN;i++)
	{
		if(chFileName[i]!=' ')
			stFileName+=chFileName[i];
		else
			break;
	}
	stPFileName=stFileName;
	ss<<*nPartition;
	stFileName+="_";
	stFileName+=ss.str();
	stFileName+=".vts";
	//std::cout<<"nTotal: "<<nTotal<<", "<<*nFlds<<", "<<*nPartition<<" "+stFileName<<"\n";
	//std::cout<<"nSize: "<<nSize[0]<<", "<<nSize[1]<<", "<<nSize[2]<<" "+stFileName<<"\n";

	//Step 1: setup grid, points and variable objects
	vtkSmartPointer<vtkStructuredGrid> sGrid =
		vtkSmartPointer<vtkStructuredGrid>::New();
	vtkSmartPointer<vtkPoints> points =
		vtkSmartPointer<vtkPoints>::New();
	vtkSmartPointer<vtkDoubleArray> velocity =
		vtkSmartPointer<vtkDoubleArray>::New();
	vtkSmartPointer<vtkDoubleArray> pressure =
		vtkSmartPointer<vtkDoubleArray>::New();
	vtkSmartPointer<vtkDoubleArray> temperature =
		vtkSmartPointer<vtkDoubleArray>::New();

	//Step 2: Set Grid dimensions and variable tuples
	sGrid->SetDimensions(nSize[0],nSize[1],nSize[2]);
	points->Allocate(nTotal);
	velocity->SetNumberOfComponents(3);
	pressure->SetNumberOfComponents(1);
	temperature->SetNumberOfComponents(1);
	velocity->SetNumberOfTuples(nTotal);
	pressure->SetNumberOfTuples(nTotal);
	temperature->SetNumberOfTuples(nTotal);
	velocity->SetName("velocity");
	pressure->SetName("pressure");
	temperature->SetName("temperature");
	
	//Step 3: Copy data into vtk Objects
	int j=0;
	for( int i=0;i<nTotal*3;i+=3)
	{
	//	std::cout<<"i="<<i<<"\n";
		points->InsertPoint(j,&dPnts[i]);
		j++;
	}

	//std::cout<<"Points allocated succesfully\n";
	j=0;
	for(int i=0;i<nTotal*(*nFlds);i+=(*nFlds))
	{
	//	std::cout<<"i="<<i<<"\n";
		velocity->InsertTuple(j,&dFlds[i]);
		pressure->InsertTuple(j,&dFlds[i+3]);
		temperature->InsertTuple(j,&dFlds[i+4]);
		j++;
	}
	//std::cout<<"Fields allocated succesfully\n";

	//Step 4: Assign values to grid
	sGrid->SetPoints(points);
	sGrid->GetPointData()->SetVectors(velocity);
	sGrid->GetPointData()->AddArray(pressure);
	sGrid->GetPointData()->SetScalars(temperature);

	//Step 5: Setup writer object and write file
		vtkSmartPointer<vtkXMLStructuredGridWriter> writer =
			vtkSmartPointer<vtkXMLStructuredGridWriter>::New();
		writer->SetFileName(&stFileName[0]);
		writer->SetInputData(sGrid);
		writer->Write();
	//Step 6: If partition zero write pvts vile
	//if(*nPartition==0)
	/*{
		stPFileName+=".pvts";
		vtkSmartPointer<vtkXMLPStructuredGridWriter> pwriter =
			vtkSmartPointer<vtkXMLPStructuredGridWriter>::New();
		pwriter->SetFileName(&stPFileName[0]);
		pwriter->SetInputData(sGrid);
		pwriter->SetNumberOfPieces(8);
		//pwriter->SetUpdateExtent(ext);
		pwriter->Write();
	}*/
}

#pragma once
#ifndef MODEL_H
#define MODEL_H
#include <iostream>
#include <vector>
#include <string>
#include <fstream>
#include <sstream>
#include <memory>
#include "vtkPoints.h"
#include "vtkLine.h"
#include "vtkVertex.h"
#include "vtkTriangle.h"

#include "vtkUnstructuredGrid.h"
#include "vtkAppendFilter.h"

#include "vtkSmartPointer.h"
#include "vtkSphereSource.h"

#include "vtkCellArray.h"
#include "vtkDataSetMapper.h"
#include "vtkActor.h"
#include "vtkProperty.h"

using namespace std;

class ContactTypes
{
private:
	unsigned stepNum;

	typedef vector<pair<unsigned, unsigned>> pairCollect;
	pairCollect node_triangle;
	pairCollect triangle_node;
	pairCollect edge_edge;
	pairCollect node_edge;
	pairCollect edge_node;
	pairCollect node_node;

	vector<vtkSmartPointer<vtkUnstructuredGrid>> unstructuredGrids;

};

class Displacement
{
private:
	vector<vtkSmartPointer<vtkPoints> > pvtkPnts;


};

class Model // base class
{
protected:
	unsigned stepNum, nodeNum, dofs;
	vector<vector<double> > dispvecCollection;
	vector<double> stepCollection;
	vector<double> node0;

	vtkPoints *pvtkPosition;
	vector<vtkSmartPointer<vtkPoints> > pvtkPnts;
public:
	virtual ~Model() {};
	Model() {};
	static shared_ptr<Model> New()
	{
		shared_ptr<Model> nw = make_shared<Model>();
		return nw;
	}

	virtual void readDispfile(const vector<string> & filename);
	virtual void initialize();
	virtual void setVtkpnt0(vtkPoints *input) { pvtkPosition = input; }
	unsigned & getNodenum() { return nodeNum; }

	unsigned & getStepNum() { return stepNum; }
	vector<double> & getStep() { return stepCollection; }
	double getStep(unsigned i) { return stepCollection[i]; }
	
	vtkPoints *getvtkPnts(unsigned i) { return pvtkPnts[i]; }

};

void Model::initialize() 
{
	node0.resize(3* nodeNum);
	unsigned index = 0;
	for (unsigned n = 0; n < nodeNum; ++n) {
		double *xyz = pvtkPosition->GetPoint(n);
		node0[index++] = *xyz;
		node0[index++] = *(xyz + 1);
		node0[index++] = *(xyz + 2);
	}

	pvtkPnts.resize(dispvecCollection.size());
	double position[3];
	for (unsigned s = 0; s < stepNum; ++s){
		pvtkPnts[s] = vtkSmartPointer<vtkPoints>::New();
		pvtkPnts[s]->SetNumberOfPoints(nodeNum);
		index = 0;
		for (unsigned n = 0; n < nodeNum; ++n) {
			position[0] = node0[index] + dispvecCollection[s][index];
			position[1] = node0[index+1] + dispvecCollection[s][index+1];
			position[2] = node0[index+2] + dispvecCollection[s][index+2];
			index += 3;
			pvtkPnts[s]->SetPoint(n, position );
		}
	}

}

void Model::readDispfile(const vector<string> & filename)
{
	ifstream infile;
	infile.open(filename[0], ios::in | ios::binary);
	if (infile.is_open()) {
		nodeNum = pvtkPosition->GetNumberOfPoints();
		unsigned nodedeg = 3;
		if (filename.size() == 2) {
			nodedeg = atoi(filename[1].c_str());
		}
		dofs = nodeNum * nodedeg;
		vector<double> datavec(dofs);
		double steptime;
		while (!infile.eof()) {
			infile.read((char*)&steptime, sizeof(double));
			infile.read((char*)datavec.data(), sizeof(double) * (dofs));
			dispvecCollection.push_back(datavec);
			stepCollection.push_back(steptime);
		}
		stepNum = (unsigned)dispvecCollection.size();

		if (filename.size() == 2) {

			unsigned index;
			vector<double> dataPosition;
			for (unsigned s = 0; s < stepNum; ++s) {
				auto &data = dispvecCollection[s];
				index = 0;
				dataPosition.clear();
				for (unsigned i = 0; i < nodeNum; ++i) {
					dataPosition.push_back(data[index]);
					dataPosition.push_back(data[index+1]);
					dataPosition.push_back(data[index+2]);
					index += nodedeg;
				}
				dispvecCollection[s].swap(dataPosition);
			}

		}
		infile.close();
	}
	else {
		cout << "Can not open dispfiles." << endl;
		exit(0);
	}



}


#endif //MODEL_H
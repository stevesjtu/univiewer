#pragma once
#ifndef MODEL_H
#define MODEL_H
#include <iostream>
#include <vector>
#include <string>
#include <fstream>
#include <memory>
#include "Eigen/eigen"
#include <vtkPoints.h>

using namespace Eigen;
using namespace std;


class Model // base class
{
protected:
	unsigned stepNum, nodeNum, dofs;
	vector<VectorXd> dispvecCollection;
	vector<double> stepCollection;
	VectorXd node0;

	vtkPoints * pvtkPosition;
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
	VectorXd & getNode0() { return node0; }
	VectorXd getNode0(unsigned id) { return node0.segment<3>(id * 3); }
	unsigned & getStepNum() { return stepNum; }
	vector<double> & getStep() { return stepCollection; }
	double getStep(unsigned i) { return stepCollection[i]; }

	vector<VectorXd> & getDispvecCollection() { return dispvecCollection; }
	VectorXd & getDispvecCollection(unsigned i) { return dispvecCollection[i]; }
    
	
	vtkPoints *getvtkPnts(unsigned i) { return pvtkPnts[i]; }

};

void Model::initialize() 
{
	node0.resize(dofs);

	for (unsigned n = 0; n < nodeNum; ++n) {
		double *xyz = pvtkPosition->GetPoint(n);
		node0.segment<3>(3 * n) = Vector3d(*xyz, *(xyz + 1), *(xyz + 2));
	}

	pvtkPnts.resize(dispvecCollection.size());
	Vector3d position;
	for (unsigned s = 0; s < stepNum; ++s){
		pvtkPnts[s] = vtkSmartPointer<vtkPoints>::New();
		pvtkPnts[s]->SetNumberOfPoints(nodeNum);
		for (unsigned n = 0; n < nodeNum; ++n) {
			position = node0.segment<3>(3 * n) + dispvecCollection[s].segment<3>(3* n);
			pvtkPnts[s]->SetPoint(n, position.data() );
		}
	}

}

void Model::readDispfile(const vector<string> & filename)
{
	ifstream infile;
	infile.open(filename[0], ios::in | ios::binary);
	nodeNum = pvtkPosition->GetNumberOfPoints();
	unsigned nodedeg = 3;
	if (filename.size() == 2) {
		nodedeg = atoi(filename[1].c_str());
	}
	dofs = nodeNum * nodedeg;
	VectorXd datavec(dofs);
	double steptime;
	while (!infile.eof()) {
		infile.read((char*)&steptime, sizeof(double));
		infile.read((char*)datavec.data(), sizeof(double) * (dofs));
		dispvecCollection.push_back(datavec);
		stepCollection.push_back(steptime);
	}
	stepNum = (unsigned)dispvecCollection.size();

	if (filename.size() == 2) {
		MatrixXd datamat;
		MatrixXd datahalf;
		for (unsigned s = 0; s < stepNum; ++s) {
			datamat = dispvecCollection[s];
			datamat.resize(nodedeg, dispvecCollection[s].size() / nodedeg);
			datahalf = datamat.topRows(3);
			datahalf.resize(nodeNum * 3, 1);
			dispvecCollection[s] = datahalf;
		}

	}

	infile.close();

}


#endif //MODEL_H
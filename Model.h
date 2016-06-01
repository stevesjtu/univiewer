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
	vector<vtkSmartPointer<vtkPoints>> pvtkPnts;
public:
	virtual ~Model() {};
	Model() {};
	static shared_ptr<Model> New()
	{
		shared_ptr<Model> nw = make_shared<Model>();
		return nw;
	}

	virtual void readDispfile(const std::string & filename);
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
	
	virtual void update(unsigned num);
};

void Model::initialize() {

	
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

void Model::readDispfile(const std::string & filename)
{
	ifstream infile;
	infile.open(filename, ios::in | ios::binary);
	nodeNum = pvtkPosition->GetNumberOfPoints();
	dofs = nodeNum * 3;
	VectorXd datavec(dofs);
	double steptime;
	while (!infile.eof()) {
		infile.read((char*)&steptime, sizeof(double));
		infile.read((char*)datavec.data(), sizeof(double) * (dofs));
		dispvecCollection.push_back(datavec);
		stepCollection.push_back(steptime);
	}

	stepNum = (unsigned)dispvecCollection.size();
	infile.close();

}

void Model::update(unsigned num)
{
	//node = node0 + dispvecCollection[num];

	//for (int n = 0; n < (int)nodeNum; ++n)
	//	pvtkPnts->SetPoint(n, node.data() + (n* 3));

}

#endif //MODEL_H
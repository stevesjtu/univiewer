#pragma once
#ifndef MODEL_H
#define MODEL_H
#include <iostream>
#include <vector>
#include <array>
#include <string>
#include <fstream>
#include <sstream>
#include <memory>
#include "vtkPoints.h"
#include "vtkLine.h"
#include "vtkVertex.h"
#include "vtkTriangle.h"

#include "vtkXMLUnstructuredGridReader.h"
#include "vtkUnstructuredGrid.h"
#include "vtkAppendFilter.h"

#include "vtkSmartPointer.h"
#include "vtkSphereSource.h"

#include "vtkCellArray.h"
#include "vtkDataSetMapper.h"
#include "vtkActor.h"
#include "vtkActor2D.h"
#include "vtkProperty.h"

#include "vtkLabeledDataMapper.h"

using namespace std;

class Model;

////////////////////////////////////////////////////////////////////////////
// labelnode
class Labelnode
{
	friend class Model;
private:
	vtkSmartPointer<vtkActor2D> labelActor;
	vtkSmartPointer<vtkLabeledDataMapper> labelMapper;
public:
	Labelnode() {};
	virtual ~Labelnode() {};

	static shared_ptr<Labelnode> New()
	{
		return make_shared<Labelnode>();
	}

	vtkSmartPointer<vtkActor2D> &getlabelActor() { return labelActor; }
	void setLabelActor(vtkSmartPointer<vtkUnstructuredGrid>&ugrid) {
		labelMapper = vtkSmartPointer<vtkLabeledDataMapper>::New();
		labelMapper->SetInputData(ugrid);
		labelActor = vtkSmartPointer<vtkActor2D>::New();
		labelActor->SetMapper(labelMapper);
	}
};

////////////////////////////////////////////////////////////////////////////
// ConstactInfo
class ContactInfo
{
	friend class Model;
private:
	typedef vector<pair<unsigned, unsigned>> pairCollect;
	pairCollect node_triangle;
	pairCollect triangle_node;
	pairCollect edge_edge;
	pairCollect node_edge;
	pairCollect edge_node;
	pairCollect node_node;

	vector<vtkSmartPointer<vtkUnstructuredGrid>> ntUGrid;
	vector<vtkSmartPointer<vtkUnstructuredGrid>> tnUGrid;
	vector<vtkSmartPointer<vtkUnstructuredGrid>> eeUGrid;
	vector<vtkSmartPointer<vtkUnstructuredGrid>> neUGrid;
	vector<vtkSmartPointer<vtkUnstructuredGrid>> enUGrid;
	vector<vtkSmartPointer<vtkUnstructuredGrid>> nnUGrid;
public:
	ContactInfo() {};
	virtual ~ContactInfo() {};
	static shared_ptr<ContactInfo> New(){ return make_shared<ContactInfo>();}
	void readContfile(const string&);
};

void ContactInfo::readContfile(const string& contfile)
{

}

////////////////////////////////////////////////////////////////////////////
// NodePositions
class NodePositions
{
	friend class Model;
private:
	vector<vtkSmartPointer<vtkPoints> > pvtkPnts;
public:
	NodePositions() {};
	virtual ~NodePositions() {};
	static shared_ptr<NodePositions> New() { return make_shared<NodePositions>(); }
};

////////////////////////////////////////////////////////////////////////////
// Displacements
class Displacements
{
	friend class Model;
private:
	vector<vector<double>> dispvecCollection;
	vector<double> stepCollection;
public:
	Displacements() {};
	virtual ~Displacements(){}
	static shared_ptr<Displacements> New() { return make_shared<Displacements>(); }
	void readDispfile(const vector<string> & filename);
};

////////////////////////////////////////////////////////////////////////////
// Model
class Model // base class
{
protected:
	unsigned ID, nodeNum, elemNum, offset;
	vtkSmartPointer<vtkUnstructuredGrid> ugrid;
	vtkSmartPointer<vtkActor> actor;
	vtkSmartPointer<vtkDataSetMapper> mapper;

	shared_ptr<ContactInfo> conactinfo;
	shared_ptr<NodePositions> nodepositions;
	shared_ptr<Labelnode> labelnode;
public:
	static unsigned count, nodeNums, elemNums, stepNum, step;
	static shared_ptr<Displacements> displacements;
	virtual ~Model() {};
	Model()
	{
		ID = count;
		++count;
	}

	static shared_ptr<Model> New()
	{
		shared_ptr<Model> nw = make_shared<Model>();
		return nw;
	}

	inline unsigned getNodenum() { return nodeNum; }
	static double getStepVal(unsigned i) { return displacements->stepCollection[i]; }
	static void setDisplacement(shared_ptr<Displacements> &disp) { displacements = disp; }

	inline vtkSmartPointer<vtkActor2D> &getLabelactor() { return labelnode->getlabelActor(); }
	inline vtkSmartPointer<vtkActor> &getActor() { return actor; }
	
	void setOffset(unsigned index) { offset = index; }
	void readModel(const string&);
	void initializeDisp();
	void setLabelnode();
	void updateDisp(unsigned);
};
unsigned Model::count = 0;
unsigned Model::nodeNums = 0;
unsigned Model::elemNums = 0;
unsigned Model::stepNum = 0;
unsigned Model::step = 0;
shared_ptr<Displacements> Model::displacements = NULL;

void Model::readModel(const string&file)
{
	// read from xml
	vtkSmartPointer<vtkXMLUnstructuredGridReader> ugridReader = vtkSmartPointer<vtkXMLUnstructuredGridReader>::New();
	ugridReader->SetFileName(file.c_str());
	ugridReader->Update();
	// 
	ugrid = vtkSmartPointer<vtkUnstructuredGrid>::New();
	ugrid = ugridReader->GetOutput();

	mapper = vtkSmartPointer<vtkDataSetMapper>::New();
	mapper->SetInputData(ugrid);
	
	actor = vtkSmartPointer<vtkActor>::New();
	actor->SetMapper(mapper);

	// modify the static variables
	nodeNums += nodeNum = ugrid->GetNumberOfPoints();
	elemNums += elemNum = ugrid->GetNumberOfCells();
}

void Model::initializeDisp()
{
	vector<double> node0;
	node0.resize(nodeNum * 3);
	unsigned index = 0;
	for (unsigned n = 0; n < nodeNum; ++n) {
		double *xyz = ugrid->GetPoint(n);
		node0[index++] = *xyz;
		node0[index++] = *(xyz + 1);
		node0[index++] = *(xyz + 2);
	}
	nodepositions = NodePositions::New();
	nodepositions->pvtkPnts.resize(displacements->dispvecCollection.size());
	double position[3];
	for (unsigned s = 0; s < stepNum; ++s) {
		nodepositions->pvtkPnts[s] = vtkSmartPointer<vtkPoints>::New();
		nodepositions->pvtkPnts[s]->SetNumberOfPoints(nodeNum);
		index = offset;
		for (unsigned n = 0; n < nodeNum; ++n) {
			position[0] = node0[n* 3] + displacements->dispvecCollection[s][index];
			position[1] = node0[n * 3 + 1] + displacements->dispvecCollection[s][index + 1];
			position[2] = node0[n * 3 + 2] + displacements->dispvecCollection[s][index + 2];
			index += 3;
			nodepositions->pvtkPnts[s]->SetPoint(n, position);
		}
	}
}

void Model::updateDisp(unsigned s)
{
	ugrid->SetPoints(nodepositions->pvtkPnts[s]);
}

void Model::setLabelnode()
{
	labelnode = Labelnode::New();
	labelnode->setLabelActor(ugrid);
}

void Displacements::readDispfile(const vector<string> & filename)
{
	ifstream infile;
	infile.open(filename[0], ios::in | ios::binary);
	if (infile.is_open()) {

		unsigned nodedeg = 3;
		if (filename.size() == 2) {
			nodedeg = atoi(filename[1].c_str());
		}
		unsigned dofs = Model::nodeNums * nodedeg;

		vector<double> datavec(dofs);
		double steptime;
		while (!infile.eof()) {
			infile.read((char*)&steptime, sizeof(double));
			infile.read((char*)datavec.data(), sizeof(double) * (dofs));
			dispvecCollection.push_back(datavec);
			stepCollection.push_back(steptime);
		}
		Model::stepNum = (unsigned)dispvecCollection.size();

		if (filename.size() == 2) {
			unsigned index;
			vector<double> dataPosition;
			for (unsigned s = 0; s < Model::stepNum; ++s) {
				auto &data = dispvecCollection[s];
				index = 0;
				dataPosition.clear();
				for (unsigned i = 0; i < Model::nodeNums; ++i) {
					dataPosition.push_back(data[index]);
					dataPosition.push_back(data[index + 1]);
					dataPosition.push_back(data[index + 2]);
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
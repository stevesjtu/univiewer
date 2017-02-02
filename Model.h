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

#define STEP_SEPARATOR 0xffffffff // the largest unsigned int in 32bit machine
#define CONTACT_SEPARATOR 0xfffffffe

#define NODE_TRIANGLE 0x0001
#define TRIANGLE_NODE 0x0002
#define EDGE_EDGE	  0x0004
#define NODE_EDGE	  0x0008
#define EDGE_NODE	  0x0010
#define NODE_NODE	  0x0020

using namespace std;

typedef vector<pair<unsigned, unsigned>> pairCollect;
class Model;

////////////////////////////////////////////////////////////////////////////
// labelnode
class Labelnode
{
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
// feMesh
class FEMesh
{
private:
	bool isCalcEdge;
	vtkSmartPointer<vtkUnstructuredGrid> ugrid;
	vector<vtkSmartPointer<vtkPoints> > pvtkPnts;
	vector<array<unsigned, 2>> edge_bynodes;
public:
	FEMesh() {};
	virtual ~FEMesh() {};
	static shared_ptr<FEMesh> New() { return make_shared<FEMesh>(); }
	inline vtkSmartPointer<vtkUnstructuredGrid> &getUGrid() { return ugrid; }
	vector<vtkSmartPointer<vtkPoints> > &getpvtkPnts() { return pvtkPnts; }
	vtkSmartPointer<vtkPoints> &getpvtkPnts(unsigned i) { return pvtkPnts[i]; }
};


////////////////////////////////////////////////////////////////////////////
// Model
class Model // base class
{
protected:
	unsigned ID, nodeNum, elemNum, offset;
	
	vtkSmartPointer<vtkActor> actor;
	vtkSmartPointer<vtkDataSetMapper> mapper;
	
	shared_ptr<FEMesh> feMesh;
	shared_ptr<Labelnode> labelnode;
public:
	static unsigned count, nodeNums, elemNums, stepNum, step;
	static vector<double> stepCollection;

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
	inline unsigned getOffset() { return offset; }
	inline vtkSmartPointer<vtkActor2D> &getLabelactor() { return labelnode->getlabelActor(); }
	inline vtkSmartPointer<vtkActor> &getActor() { return actor; }
	inline shared_ptr<FEMesh> &getFEMesh() { return feMesh; }
	
	void setOffset(unsigned index) { offset = index; }
	void readModel(const string&);
	void setLabelnode();
	void updateDisp(unsigned);
};

unsigned Model::count = 0;
unsigned Model::nodeNums = 0;
unsigned Model::elemNums = 0;
unsigned Model::stepNum = 0;
unsigned Model::step = 0;
vector<double> Model::stepCollection = vector<double>(0);

void Model::readModel(const string&file)
{

	feMesh = FEMesh::New();

	// read from xml
	vtkSmartPointer<vtkXMLUnstructuredGridReader> ugridReader = vtkSmartPointer<vtkXMLUnstructuredGridReader>::New();
	ugridReader->SetFileName(file.c_str());
	ugridReader->Update();
	// 
	feMesh->getUGrid() = vtkSmartPointer<vtkUnstructuredGrid>::New();
	feMesh->getUGrid() = ugridReader->GetOutput();
	
	mapper = vtkSmartPointer<vtkDataSetMapper>::New();
	mapper->SetInputData(feMesh->getUGrid());
	
	actor = vtkSmartPointer<vtkActor>::New();
	actor->SetMapper(mapper);

	// modify the static variables
	nodeNums += nodeNum = feMesh->getUGrid()->GetNumberOfPoints();
	elemNums += elemNum = feMesh->getUGrid()->GetNumberOfCells();
}

void Model::updateDisp(unsigned s)
{
	feMesh->getUGrid()->SetPoints(feMesh->getpvtkPnts(s));
}

void Model::setLabelnode()
{
	labelnode = Labelnode::New();
	labelnode->setLabelActor(feMesh->getUGrid());
}



////////////////////////////////////////////////////////////////////////////
// ConstactInfo
class ContactData
{
private:
	shared_ptr<Model> pModel, nModel;

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
	ContactData() {};
	virtual ~ContactData() {};
	static shared_ptr<ContactData> New() { return make_shared<ContactData>(); }

	void setModels(shared_ptr<Model>& p, shared_ptr<Model> &n) { pModel = p; nModel = n; }
	pairCollect &getPrimitivesPairs(int type)
	{
		switch (type)
		{
		case NODE_TRIANGLE:
			return node_triangle;
		case TRIANGLE_NODE:
			return triangle_node;
		case EDGE_EDGE:
			return edge_edge;
		case NODE_EDGE:
			return node_edge;
		case EDGE_NODE:
			return edge_node;
		case NODE_NODE:
			return node_node;
		default:
			break;
		}
		return node_triangle;
	}

	void InitializeUGrid();
};

void ContactData::InitializeUGrid()
{

}

//////////////////////////////////////////////////////////////////////////////////////////////
//
// read disp files, and read contact files
//
//
//////////////////////////////////////////////////////////////////////////////////////////////
void readDispfile(const vector<string> & filename, vector<shared_ptr<Model>> &pModels)
{
	ifstream infile;
	infile.open(filename[0], ios::in | ios::binary);
	
	if (infile.is_open()) {
		vector<vector<double>> dispvecCollection;
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
			Model::stepCollection.push_back(steptime);
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

		// initialize the nodes position of Model feMesh
		vector<double> node0;
		for (auto& pmodel : pModels) {
			auto &feMesh = pmodel->getFEMesh();

			node0.resize(pmodel->getNodenum() * 3);
			unsigned index = 0;
			for (unsigned n = 0; n < pmodel->getNodenum(); ++n) {
				double *xyz = feMesh->getUGrid()->GetPoint(n);
				node0[index++] = *xyz;
				node0[index++] = *(xyz + 1);
				node0[index++] = *(xyz + 2);
			}

			feMesh->getpvtkPnts().resize(Model::stepNum);

			double position[3];
			for (unsigned s = 0; s < Model::stepNum; ++s) {
				feMesh->getpvtkPnts(s) = vtkSmartPointer<vtkPoints>::New();
				feMesh->getpvtkPnts(s)->SetNumberOfPoints(pmodel->getNodenum());
				index = pmodel->getOffset();
				for (unsigned n = 0; n < pmodel->getNodenum(); ++n) {
					position[0] = node0[n * 3]     + dispvecCollection[s][index];
					position[1] = node0[n * 3 + 1] + dispvecCollection[s][index + 1];
					position[2] = node0[n * 3 + 2] + dispvecCollection[s][index + 2];
					index += 3;
					feMesh->getpvtkPnts(s)->SetPoint(n, position);
				}
			}
		}
	} // if file opened
	else {
		cout << "Can not open dispfiles." << endl;
		exit(0);
	}

}


inline void readpairs(ifstream &infile, pairCollect &pc)
{
	unsigned num, pairs[2];
	infile.read((char*)&num, sizeof(unsigned int));
	for (unsigned i = 0; i < num; ++i) {
		infile.read((char*)pairs, sizeof(unsigned int) * 2);
		pc.push_back(make_pair(pairs[0], pairs[1]));
	}
}

void readContfile(const string& contfile, vector<vector<shared_ptr<ContactData>>> &pContactss, vector<shared_ptr<Model>> &pModels)
{
	ifstream infile;
	infile.open(contfile, ios::in | ios::binary);
	if (infile.is_open()) {
		pContactss.resize(Model::stepNum);
		unsigned bodyid1, bodyid2, mark;
		shared_ptr<ContactData> contactTemp = ContactData::New();

		for (unsigned s = 0; s < Model::stepNum; ++s) {
			infile.read((char*)&mark, sizeof(unsigned int));
			while (mark == CONTACT_SEPARATOR) {
				infile.read((char*)&bodyid1, sizeof(unsigned int));
				infile.read((char*)&bodyid2, sizeof(unsigned int));
				contactTemp->setModels(pModels[bodyid1], pModels[bodyid2]);
				// primitives
				readpairs(infile, contactTemp->getPrimitivesPairs(NODE_TRIANGLE));
				readpairs(infile, contactTemp->getPrimitivesPairs(TRIANGLE_NODE));
				readpairs(infile, contactTemp->getPrimitivesPairs(EDGE_EDGE));
				readpairs(infile, contactTemp->getPrimitivesPairs(NODE_EDGE));
				readpairs(infile, contactTemp->getPrimitivesPairs(EDGE_NODE));
				readpairs(infile, contactTemp->getPrimitivesPairs(NODE_NODE));
				// add this ContactData
				pContactss[s].push_back(contactTemp);

				infile.read((char*)&mark, sizeof(unsigned int));
				if (mark == STEP_SEPARATOR) break;
			} // loop for ContactData
		} // loop for steps
	}
	else {
		cout << "Can not open contfile." << endl;
		exit(0);
	}

}






#endif //MODEL_H
#pragma once
#ifndef MODEL_H
#define MODEL_H

#include "ParamDefine.h"
#include <vector>
#include <set>
#include <array>
#include <string>
#include <fstream>
#include <sstream>
#include <iomanip>
// primitives
#include "vtkPoints.h"
#include "vtkLine.h"
#include "vtkVertex.h"
#include "vtkTriangle.h"
#include "vtkCellArray.h"
#include "vtkCellData.h"
#include "vtkUnsignedCharArray.h"
// data structure
#include "vtkXMLUnstructuredGridReader.h"
#include "vtkUnstructuredGrid.h"
#include "vtkDataSetMapper.h"
#include "vtkLabeledDataMapper.h"
// actor hold by model
#include "vtkActor.h"
#include "vtkActor2D.h"

using namespace std;

typedef vector<pair<unsigned, unsigned>> pairCollect;
class Model;

////////////////////////////////////////////////////////////////////////////
// labelnode
class Labelnode
{
private:
	vtkSmartPointer<vtkActor2D> labelActor;
public:
	Labelnode() {};
	virtual ~Labelnode() {};

	static shared_ptr<Labelnode> New()
	{
		return make_shared<Labelnode>();
	}

	vtkSmartPointer<vtkActor2D> &getlabelActor() { return labelActor; }
	void setLabelActor(vtkSmartPointer<vtkUnstructuredGrid>&ugrid) {
		vtkSmartPointer<vtkLabeledDataMapper> labelMapper = vtkSmartPointer<vtkLabeledDataMapper>::New();
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
	vector<array<unsigned, 2>> edges_bynodes;
	vector<array<unsigned, 3>> elems_bynodes;
public:
	FEMesh():isCalcEdge(false) {};
	virtual ~FEMesh() {};
	static shared_ptr<FEMesh> New() { return make_shared<FEMesh>(); }
	inline vtkSmartPointer<vtkUnstructuredGrid> &getUGrid() { return ugrid; }
	vector<vtkSmartPointer<vtkPoints> > &getpvtkPnts() { return pvtkPnts; }
	vtkSmartPointer<vtkPoints> &getpvtkPnts(unsigned i) { return pvtkPnts[i]; }
	vector<array<unsigned, 2>> &getEdges() { return edges_bynodes; }
	vector<array<unsigned, 3>> &getElems() { return elems_bynodes; }
	array<unsigned, 2> &getEdges(unsigned i) { return edges_bynodes[i]; }
	array<unsigned, 3> &getElems(unsigned i) { return elems_bynodes[i]; }

	bool getIsCalcEdge() { return isCalcEdge; }
	void makeEdges();
};

void FEMesh::makeEdges()
{
	// edgeSet should be a set of edges without duplicated edges.
	unsigned elem_num = ugrid->GetNumberOfCells();
	elems_bynodes.resize(elem_num);

	vtkSmartPointer<vtkCellArray> cellarray = ugrid->GetCells();
	vtkSmartPointer<vtkIdList> idlist = vtkSmartPointer<vtkIdList>::New();
	for (unsigned e = 0; e < elem_num; ++e) {
		cellarray->GetNextCell(idlist);
		elems_bynodes[e] = { (unsigned)idlist->GetId(0), (unsigned)idlist->GetId(1), (unsigned)idlist->GetId(2) };
	}
	
	typedef set<unsigned> edge;
	edge oneEdge;
	set<edge> edgeSet;
	unsigned ind1;
	for (unsigned e = 0; e < elem_num; ++e) {
		for (unsigned ind = 0; ind < 3; ++ind) {
			ind == 3 - 1 ? ind1 = 0 : ind1 = ind + 1;
			oneEdge.insert(elems_bynodes[e][ind]);
			oneEdge.insert(elems_bynodes[e][ind1]);
			edgeSet.insert(oneEdge);
			oneEdge.clear();
		}
	}

	//// fill the edges_bynodes
	unsigned edge_num = static_cast<unsigned>(edgeSet.size());
	edges_bynodes.resize(edge_num);
	unsigned index = 0;
	for (set<edge>::iterator it = edgeSet.begin(); it != edgeSet.end(); it++) {
		edge::iterator nodeit = it->begin();
		edges_bynodes[index][0] = *nodeit++;
		edges_bynodes[index++][1] = *nodeit;
	}

	isCalcEdge = true;
}


////////////////////////////////////////////////////////////////////////////
// Model
class Model // base class
{
protected:
	unsigned ID, nodeNum, elemNum, offset;
	
	vtkSmartPointer<vtkActor> actor;
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
	
	vtkSmartPointer<vtkDataSetMapper> mapper = vtkSmartPointer<vtkDataSetMapper>::New();
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

	vector<pairCollect> node_triangles;
	vector<pairCollect> triangle_nodes;
	vector<pairCollect> edge_edges;
	vector<pairCollect> node_edges;
	vector<pairCollect> edge_nodes;
	vector<pairCollect> node_nodes;

	vector<vtkSmartPointer<vtkUnstructuredGrid>> pUGrids;
	vector<vtkSmartPointer<vtkUnstructuredGrid>> nUGrids;

	//vtkSmartPointer<vtkUnstructuredGrid > pActiveUGrid;
	//vtkSmartPointer<vtkUnstructuredGrid > nActiveUGrid;

	vector<vtkSmartPointer<vtkDataSetMapper>> pMappers;
	vector<vtkSmartPointer<vtkDataSetMapper>> nMappers;
	
	vtkSmartPointer<vtkActor> pActor;
	vtkSmartPointer<vtkActor> nActor;

	inline void insertNode(	
		unsigned n, 
		vtkSmartPointer<vtkCellArray> &Cells,
		vtkSmartPointer<vtkUnsignedCharArray> &Colors,
		vector<int> &Types);
	inline void insertEdge(
		array<unsigned, 2> edge,
		vtkSmartPointer<vtkCellArray> &Cells,
		vtkSmartPointer<vtkUnsignedCharArray> &Colors,
		vector<int> &Types);
	inline void insertTriangle(
		array<unsigned, 3> elem,
		vtkSmartPointer<vtkCellArray> &Cells,
		vtkSmartPointer<vtkUnsignedCharArray> &Colors,
		vector<int> &Types);

public:
	ContactData() {};
	virtual ~ContactData() {};
	static shared_ptr<ContactData> New() { return make_shared<ContactData>(); }

	void setModels(shared_ptr<Model>& p, shared_ptr<Model> &n) { pModel = p; nModel = n; }
	pairCollect &getPrimitivesPairs(int type, unsigned step)
	{
		switch (type)
		{
		case NODE_TRIANGLE:
			return node_triangles[step];
		case TRIANGLE_NODE:
			return triangle_nodes[step];
		case EDGE_EDGE:
			return edge_edges[step];
		case NODE_EDGE:
			return node_edges[step];
		case EDGE_NODE:
			return edge_nodes[step];
		case NODE_NODE:
			return node_nodes[step];
		default:
			break;
		}
		return node_triangles[step];
	}

	vector<pairCollect> &getPrimitivesPairs(int type)
	{
		switch (type)
		{
		case NODE_TRIANGLE:
			return node_triangles;
		case TRIANGLE_NODE:
			return triangle_nodes;
		case EDGE_EDGE:
			return edge_edges;
		case NODE_EDGE:
			return node_edges;
		case EDGE_NODE:
			return edge_nodes;
		case NODE_NODE:
			return node_nodes;
		default:
			break;
		}
		return node_triangles;
	}

	void InitializeUGrid();
	void UpdateUGrid(unsigned s);
	vtkSmartPointer<vtkActor> &getPActor() { return pActor; }
	vtkSmartPointer<vtkActor> &getNActor() { return nActor; }
};

void ContactData::insertNode(
	unsigned n,
	vtkSmartPointer<vtkCellArray> &Cells,
	vtkSmartPointer<vtkUnsignedCharArray> &Colors,
	vector<int> &Types)
{
	vtkSmartPointer<vtkVertex> node = vtkSmartPointer<vtkVertex>::New();
	node->GetPointIds()->SetId(0, n);
	// cell color type
	Cells->InsertNextCell(node);
	Colors->InsertNextTupleValue(red);
	Types.push_back(VTK_VERTEX);
}

void ContactData::insertEdge(
	array<unsigned, 2> edge,
	vtkSmartPointer<vtkCellArray> &Cells,
	vtkSmartPointer<vtkUnsignedCharArray> &Colors,
	vector<int> &Types)
{
	vtkSmartPointer<vtkLine> line = vtkSmartPointer<vtkLine>::New();
	line->GetPointIds()->SetId(0, edge[0]);
	line->GetPointIds()->SetId(1, edge[1]);
	// cell color type
	Cells->InsertNextCell(line);
	Colors->InsertNextTupleValue(blue);
	Types.push_back(VTK_LINE);
}

void ContactData::insertTriangle(
	array<unsigned, 3> elem,
	vtkSmartPointer<vtkCellArray> &Cells,
	vtkSmartPointer<vtkUnsignedCharArray> &Colors,
	vector<int> &Types)
{
	vtkSmartPointer<vtkTriangle> triangle = vtkSmartPointer<vtkTriangle>::New();
	triangle->GetPointIds()->SetId(0, elem[0]);
	triangle->GetPointIds()->SetId(1, elem[1]);
	triangle->GetPointIds()->SetId(2, elem[2]);
	//nCells->InsertNextCell(nUGrid0->GetCell(t1));
	// cell color type
	Cells->InsertNextCell(triangle);
	Colors->InsertNextTupleValue(green);
	Types.push_back(VTK_TRIANGLE);
	
}

void ContactData::InitializeUGrid()
{
	if (!pModel->getFEMesh()->getIsCalcEdge()) pModel->getFEMesh()->makeEdges();
	if (!nModel->getFEMesh()->getIsCalcEdge()) nModel->getFEMesh()->makeEdges();
		
	auto &pfeMesh = pModel->getFEMesh();
	auto &nfeMesh = nModel->getFEMesh();

	pUGrids.resize(Model::stepNum);
	nUGrids.resize(Model::stepNum);

	pMappers.resize(Model::stepNum);
	nMappers.resize(Model::stepNum);

	vtkSmartPointer<vtkCellArray> pCells;
	vtkSmartPointer<vtkCellArray> nCells;
	vtkSmartPointer<vtkUnsignedCharArray> pColors;
	vtkSmartPointer<vtkUnsignedCharArray> nColors;
	vector<int> pTypes, nTypes;

	for (unsigned s = 0; s < Model::stepNum; ++s) {
		pUGrids[s] = vtkSmartPointer<vtkUnstructuredGrid>::New();
		nUGrids[s] = vtkSmartPointer<vtkUnstructuredGrid>::New();

		pUGrids[s]->SetPoints(pfeMesh->getpvtkPnts(s));
		nUGrids[s]->SetPoints(nfeMesh->getpvtkPnts(s));
		
		// some media variables/ cell, color, type
		pCells = vtkSmartPointer<vtkCellArray>::New();
		nCells = vtkSmartPointer<vtkCellArray>::New();
		pColors = vtkSmartPointer<vtkUnsignedCharArray>::New();
		nColors = vtkSmartPointer<vtkUnsignedCharArray>::New();
		pColors->SetNumberOfComponents(3);
		nColors->SetNumberOfComponents(3);
		
		pTypes.clear();
		nTypes.clear();
		// add the primitives to UnstructuredGurd
		unsigned n0, n1, e0, e1, t0, t1;
		for (unsigned i = 0; i < node_triangles[s].size(); ++i) {
			n0 = node_triangles[s][i].first;
			t1 = node_triangles[s][i].second;
			insertNode(n0, pCells, pColors, pTypes);
			insertTriangle(nfeMesh->getElems(t1), nCells, nColors, nTypes);
		}

		for (unsigned i = 0; i < triangle_nodes[s].size(); ++i) {
			t0 = triangle_nodes[s][i].first;
			n1 = triangle_nodes[s][i].second;
			insertTriangle(pfeMesh->getElems(t0), pCells, pColors, pTypes);
			insertNode(n1, nCells, nColors, nTypes);
		}

		for (unsigned i = 0; i < edge_edges[s].size(); ++i) {
			e0 = edge_edges[s][i].first;
			e1 = edge_edges[s][i].second;
			insertEdge(pfeMesh->getEdges(e0), pCells, pColors, pTypes);
			insertEdge(nfeMesh->getEdges(e1), nCells, nColors, nTypes);
		}

		for (unsigned i = 0; i < node_edges[s].size(); ++i) {
			n0 = node_edges[s][i].first;
			e1 = node_edges[s][i].second;
			insertNode(n0, pCells, pColors, pTypes);
			insertEdge(nfeMesh->getEdges(e1), nCells, nColors, nTypes);
		}

		for (unsigned i = 0; i < edge_nodes[s].size(); ++i) {
			e0 = edge_nodes[s][i].first;
			n1 = edge_nodes[s][i].second;
			insertEdge(pfeMesh->getEdges(e0), pCells, pColors, pTypes);
			insertNode(n1, nCells, nColors, nTypes);
		}

		for (unsigned i = 0; i < node_nodes[s].size(); ++i) {
			n0 = node_nodes[s][i].first;
			n1 = node_nodes[s][i].second;
			insertNode(n0, pCells, pColors, pTypes);
			insertNode(n1, nCells, nColors, nTypes);
		}

		// set the UnstructuredGrid
		pUGrids[s]->SetCells(pTypes.data(), pCells);
		pUGrids[s]->GetCellData()->SetScalars(pColors);

		nUGrids[s]->SetCells(nTypes.data(), nCells);
		nUGrids[s]->GetCellData()->SetScalars(nColors);

		pMappers[s] = vtkSmartPointer<vtkDataSetMapper>::New();
		pMappers[s]->SetInputData(pUGrids[s]);

		nMappers[s] = vtkSmartPointer<vtkDataSetMapper>::New();
		nMappers[s]->SetInputData(nUGrids[s]);
	}

	// add them to pipeline
	//vtkSmartPointer<vtkDataSetMapper> pMapper = vtkSmartPointer<vtkDataSetMapper>::New();
	//pActiveUGrid = vtkSmartPointer<vtkUnstructuredGrid>::New();
 //   pActiveUGrid->DeepCopy(pUGrids[0]);
	//pMapper->SetInputData(pActiveUGrid);
	////pMapper->SetInputData(pUGrids[0]);
	pActor = vtkSmartPointer<vtkActor>::New();
	pActor->SetMapper(pMappers[0]);
	
	//vtkSmartPointer<vtkDataSetMapper> nMapper = vtkSmartPointer<vtkDataSetMapper>::New();
	//nActiveUGrid = vtkSmartPointer<vtkUnstructuredGrid>::New();
	//nActiveUGrid->DeepCopy(nUGrids[0]);
	//nMapper->SetInputData(nActiveUGrid);
	//nMapper->SetInputData(nUGrids[0]);
	nActor = vtkSmartPointer<vtkActor>::New();
	nActor->SetMapper(nMappers[0]);

}

void ContactData::UpdateUGrid(unsigned s)
{

	//pActiveUGrid->Reset();
	//nActiveUGrid->Reset();

	//pActiveUGrid->DeepCopy(pUGrids[s]);
	//nActiveUGrid->DeepCopy(nUGrids[s]);

	pActor->SetMapper(pMappers[s]);
	nActor->SetMapper(nMappers[s]);
	// //////////////////////////////////////
	// debug pUGrids nUGrids point coordinates
	/////////////////////////////////////////
	//////////////////////////////////////////////////////////////////////////
	//ofstream outfile;
	//outfile.open("pnodes.txt", ios::app);
	//outfile << "step: " << Model::step << endl;
	//for (unsigned i = 0; i < pActiveUGrid->GetNumberOfPoints(); ++i) {
	//	outfile << setprecision(15) << pActiveUGrid->GetPoint(i)[0] << '\t' << pActiveUGrid->GetPoint(i)[1] << '\t' << pActiveUGrid->GetPoint(i)[2] << endl;
	//}
	//outfile << endl;
	//outfile.close();

	//outfile.open("nnodes.txt", ios::app);
	//outfile << "step: " << Model::step << endl;
	//for (unsigned i = 0; i < nActiveUGrid->GetNumberOfPoints(); ++i) {
	//	outfile << setprecision(15) << nActiveUGrid->GetPoint(i)[0] << '\t' << nActiveUGrid->GetPoint(i)[1] << '\t' << nActiveUGrid->GetPoint(i)[2] << endl;
	//}
	//outfile << endl;
	//outfile.close();
	//////////////////////////////////////////////////////////////////////////

	//for (unsigned i = 0; i < pActiveUGrid->GetNumberOfPoints(); ++i) {
	//	printf("%f, %f, %f\n", pActiveUGrid->GetPoint(i)[0], pActiveUGrid->GetPoint(i)[1], pActiveUGrid->GetPoint(i)[2]);
	//}
	//printf("\n");

	//for (unsigned i = 0; i < nActiveUGrid->GetNumberOfPoints(); ++i) {
	//	printf("%f, %f, %f\n", nActiveUGrid->GetPoint(i)[0], nActiveUGrid->GetPoint(i)[1], nActiveUGrid->GetPoint(i)[2]);
	//}
	//printf("\n");
	
	//ofstream outfile;
	//outfile.open("ptype.txt", ios::app);
	//outfile << "step: " << Model::step << endl;
	//for (unsigned i = 0; i < pActiveUGrid->GetNumberOfCells(); ++i) {
	//	outfile << pActiveUGrid->GetCellType(i)<< '\t';
	//}
	//outfile << endl;
	//outfile.close();

	//outfile.open("ntype.txt", ios::app);
	//outfile << "step: " << Model::step << endl;
	//for (unsigned i = 0; i < nActiveUGrid->GetNumberOfCells(); ++i) {
	//	outfile << nActiveUGrid->GetCellType(i) << '\t';
	//}
	//outfile << endl;
	//outfile.close();

	//for (unsigned i = 0; i < pActiveUGrid->GetNumberOfCells(); ++i) {
	//	printf("%u\t", pActiveUGrid->GetCellType(i));
	//}
	//printf("\n");

	//for (unsigned i = 0; i < nActiveUGrid->GetNumberOfCells(); ++i) {
	//	printf("%u\t", nActiveUGrid->GetCellType(i));
	//}
	//printf("\n");


	//pActiveUGrid->SetPoints(pUGrids[s]->GetPoints());
	//nActiveUGrid->SetPoints(nUGrids[s]->GetPoints());

	//pActiveUGrid->SetCells(pTypes[s].data(), pUGrids[s]->GetCells());
	//pActiveUGrid->GetCellData()->SetScalars(pUGrids[s]->GetCellData()->GetScalars());

	//nActiveUGrid->SetCells(nTypes[s].data(), nUGrids[s]->GetCells());
	//nActiveUGrid->GetCellData()->SetScalars(nUGrids[s]->GetCellData()->GetScalars());

	//pMapper->SetInputData(pUGrids[s]);
	//nMapper->SetInputData(nUGrids[s]);


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
		dispvecCollection.pop_back();
		Model::stepCollection.pop_back();

		//////////////////////////////////////////////////////
		//dispvecCollection.erase(dispvecCollection.begin(), dispvecCollection.begin()+45);
		//Model::stepCollection.erase(Model::stepCollection.begin(), Model::stepCollection.begin() + 45);

		//dispvecCollection.erase(dispvecCollection.begin()+3, dispvecCollection.end());
		//Model::stepCollection.erase(Model::stepCollection.begin()+3, Model::stepCollection.end());
		//////////////////////////////////////////////////////

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

void readContfile(const string& contfile, vector<shared_ptr<ContactData>> &pContacts, vector<shared_ptr<Model>> &pModels)
{
	ifstream infile;
	infile.open(contfile, ios::in | ios::binary);
	if (infile.is_open()) {
		pContacts.resize(Model::stepNum);
		unsigned bodyid1, bodyid2, mark;

		infile.read((char*)&mark, sizeof(unsigned int));
		pContacts.resize(mark);
		for (unsigned c = 0; c < mark; ++c) {
			pContacts[c] = ContactData::New();
			infile.read((char*)&bodyid1, sizeof(unsigned int));
			infile.read((char*)&bodyid2, sizeof(unsigned int));
			pContacts[c]->setModels(pModels[bodyid1], pModels[bodyid2]);
			pContacts[c]->getPrimitivesPairs(NODE_TRIANGLE).resize(Model::stepNum);
			pContacts[c]->getPrimitivesPairs(TRIANGLE_NODE).resize(Model::stepNum);
			pContacts[c]->getPrimitivesPairs(EDGE_EDGE).resize(Model::stepNum);
			pContacts[c]->getPrimitivesPairs(NODE_EDGE).resize(Model::stepNum);
			pContacts[c]->getPrimitivesPairs(EDGE_NODE).resize(Model::stepNum);
			pContacts[c]->getPrimitivesPairs(NODE_NODE).resize(Model::stepNum);
		}
		
		//////////////////////////////////////////////////////////////
		//pairCollect skipBytes;
		//for (unsigned s = 0; s< 45; ++s){
		//	for (unsigned c = 0; c < pContacts.size(); ++c) {
		//		// primitives
		//		readpairs(infile, skipBytes);
		//		readpairs(infile, skipBytes);
		//		readpairs(infile, skipBytes);
		//		readpairs(infile, skipBytes);
		//		readpairs(infile, skipBytes);
		//		readpairs(infile, skipBytes);
		//	} // loop for ContactData
		//}
		//////////////////////////////////////////////////////////////

		for (unsigned s = 0; s < Model::stepNum; ++s) {
			for (unsigned c = 0; c < pContacts.size(); ++c) {
				// primitives
				readpairs(infile, pContacts[c]->getPrimitivesPairs(NODE_TRIANGLE, s));
				readpairs(infile, pContacts[c]->getPrimitivesPairs(TRIANGLE_NODE, s));
				readpairs(infile, pContacts[c]->getPrimitivesPairs(EDGE_EDGE, s));
				readpairs(infile, pContacts[c]->getPrimitivesPairs(NODE_EDGE, s));
				readpairs(infile, pContacts[c]->getPrimitivesPairs(EDGE_NODE, s));
				readpairs(infile, pContacts[c]->getPrimitivesPairs(NODE_NODE, s));
			} // loop for ContactData
		} // loop for steps
	}
	else {
		cout << "Can not open contfile." << endl;
		exit(0);
	}

}






#endif //MODEL_H
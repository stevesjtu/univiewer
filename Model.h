#pragma once
#ifndef MODEL_H
#define MODEL_H

#include "paramdefine.h"
#include <vector>
#include <set>
#include <array>
#include <string>
#include <fstream>
#include <sstream>
// primitives
#include "vtkPoints.h"
#include "vtkLine.h"
#include "vtkVertex.h"
#include "vtkTriangle.h"
#include "vtkHexahedron.h"
#include "vtkQuadraticHexahedron.h"
#include "vtkCellArray.h"
#include "vtkCellData.h"
#include "vtkPointData.h"
#include "vtkUnsignedCharArray.h"
#include "vtkDoubleArray.h"
// data structure
#include "vtkXMLUnstructuredGridReader.h"
#include "vtkUnstructuredGrid.h"
#include "vtkDataSetMapper.h"
#include "vtkLabeledDataMapper.h"
// actor hold by model
#include "vtkActor.h"
#include "vtkActor2D.h"

using namespace std;

template<typename T>
std::ostream & operator<<(std::ostream & s, const std::vector<T> &vec) {
  if (vec.empty()) return s;
  s << "(";
  for(unsigned i=0; i < vec.size()-1; ++i){
    s << vec[i] << ", ";
  }
  s << vec[vec.size()-1] << ")";
  return s;
}

typedef vector<pair<unsigned, unsigned>> pairCollect;
class Model;

// Mdfile means markdown file, 
// which use # ## ### ... symbol to represent titles
class Mdfile {
public:
  Mdfile() {};
  Mdfile(const std::string &filename) {
    this->Open(filename);
  }
  virtual~Mdfile() {
    this->Close();
  };

  void Open(const std::string &filename ) {
    infile.open(filename, std::ios::in);
    if (!infile.is_open()) std::cout << "Error in open file: " << filename << std::endl;
  }
  void Close() {
    infile.close();
  }

  virtual std::vector<double> GetDoubleArrayFrom(const std::string &title) {
    JumpTo(title);
    std::vector<double> double_list;
    std::string temp_str;
    while (infile >> temp_str) {
      if (temp_str[0] == '#') break;
      double_list.push_back(std::stod(temp_str));
    }
    return double_list;
  }

  virtual std::vector<int> GetIntArrayFrom(const std::string &title) {
    JumpTo(title);
    std::vector<int> int_list;
    std::string temp_str;
    while (infile >> temp_str) {
      if (temp_str[0] == '#') break;
      int_list.push_back(std::stoi(temp_str));
    }
    return int_list;
  }

  virtual std::vector<std::string> GetStringFrom(const std::string &title) {
    JumpTo(title);
    std::vector<std::string> str_list;
    std::string temp_str;
    while (infile >> temp_str) {
      if (temp_str[0] == '#') break;
      str_list.push_back(temp_str);
    }
    return str_list;
  }

  void tryit() {
    infile.seekg(0, std::ios::beg);
    std::string str("");
    infile >> str;

    std::cout<< infile.tellg() <<std::endl;
    std::cout<< str << std::endl;
  }

protected:
  std::ifstream infile;
  std::streampos pos1;

  void JumpTo(const std::string &title) {
    infile.seekg(0, std::ios::beg);
    std::string str("");
    while(std::getline(infile, str)) { 
      if (str[0] == '#') {
        str = str.substr(1);
        trim(str);
        if (str.compare(title) == 0) {
          //pos1 = infile.tellg();
          break;
        }
      }
    } 
  }

  void trim(std::string &s) {
    if (s.empty()) return;
    s.erase(0,s.find_first_not_of(" "));
    s.erase(s.find_last_not_of(" ") + 1);
  }

};


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
	vtkSmartPointer<vtkDataSetMapper> mapper;

	shared_ptr<FEMesh> feMesh;
	shared_ptr<Labelnode> labelnode;

	vector<vtkSmartPointer<vtkDoubleArray>> nodalscalars;
	
public:
	static unsigned count, nodeNums, elemNums, stepNum, step;
	static vector<double> stepCollection;
	static vector<pair<double, double>> ranges;

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
	inline vtkSmartPointer<vtkActor2D> getLabelactor() { return labelnode->getlabelActor(); }
	inline vtkSmartPointer<vtkActor> getActor() { return actor; }
	inline vtkSmartPointer<vtkDataSetMapper> getMapper() { return mapper; }

	inline shared_ptr<FEMesh> getFEMesh() { return feMesh; }
	inline vector<vtkSmartPointer<vtkDoubleArray>> &getNodelScalars() { return nodalscalars; }
	inline vtkSmartPointer<vtkDoubleArray> &getNodelScalars(int s) { return nodalscalars[s]; }

	void setOffset(unsigned index) { offset = index; }
	void ReadXmlModel(const string&);
  void ReadTxtModel(const string&);

  void CreateModel(const vector<int> &modelinfo,
                  const vector<int> &elemlist, const vector<double> &nodelist);

	void setLabelnode();
	void updateDisp(unsigned);
};

unsigned Model::count = 0;
unsigned Model::nodeNums = 0;
unsigned Model::elemNums = 0;
unsigned Model::stepNum = 0;
unsigned Model::step = 0;
vector<double> Model::stepCollection = vector<double>(0);
vector<pair<double, double>> Model::ranges = vector<pair<double, double>>(0);


void Model::ReadXmlModel(const string&file)
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

void Model::ReadTxtModel(const string& file) {
  feMesh = FEMesh::New();

  std::vector<int> modelinfo;
  std::vector<int> elemlist;
  std::vector<double> nodelist;

  if (file.substr(file.size() - 3).compare(".md") == 0) {
	  Mdfile txtfile(file);
	  modelinfo = txtfile.GetIntArrayFrom("model info");
	  elemlist = txtfile.GetIntArrayFrom("element list");
	  nodelist = txtfile.GetDoubleArrayFrom("node list");
	  txtfile.Close();
  } else if (file.substr(file.size() - 3).compare("txt") == 0 ){
    ifstream txtfile(file);
    if (!txtfile.is_open()) std::cout << "Error in open file: " << file << std::endl;
    unsigned num_elem, num_node, nofe;
    txtfile >> num_elem >> num_node >> nofe;
    modelinfo.push_back(num_elem);
    modelinfo.push_back(num_node);
    modelinfo.push_back(nofe);

    elemlist.resize(nofe* num_elem);
    nodelist.resize(3 * num_node);

    for (unsigned eln = 0; eln < num_elem* nofe; ++eln) {
      txtfile >> elemlist[eln];
      elemlist[eln]--;
    }

    for (unsigned ndn = 0; ndn < 3 * num_node; ++ndn) {
      txtfile >> nodelist[ndn];
    }
    txtfile.close();
  }

  this->CreateModel(modelinfo, elemlist, nodelist);

}


void Model::setLabelnode()
{
	labelnode = Labelnode::New();
	labelnode->setLabelActor(feMesh->getUGrid());
}

void Model::updateDisp(unsigned s)
{
	feMesh->getUGrid()->SetPoints(feMesh->getpvtkPnts(s));

	if (!nodalscalars.empty()) {
		feMesh->getUGrid()->GetPointData()->SetScalars(nodalscalars[s]);
		mapper->SetScalarRange(ranges[s].first, ranges[s].second);
	}

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
	Colors->InsertNextTypedTuple(red);
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
	Colors->InsertNextTypedTuple(blue);
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
	Colors->InsertNextTypedTuple(green);
	Types.push_back(VTK_TRIANGLE);
	
}

void ContactData::InitializeUGrid()
{
	if (!pModel->getFEMesh()->getIsCalcEdge()) pModel->getFEMesh()->makeEdges();
	if (!nModel->getFEMesh()->getIsCalcEdge()) nModel->getFEMesh()->makeEdges();
		
	const auto &pfeMesh = pModel->getFEMesh();
	const auto &nfeMesh = nModel->getFEMesh();

	pMappers.resize(Model::stepNum);
	nMappers.resize(Model::stepNum);

	vtkSmartPointer<vtkCellArray> pCells;
	vtkSmartPointer<vtkCellArray> nCells;
	vtkSmartPointer<vtkUnsignedCharArray> pColors;
	vtkSmartPointer<vtkUnsignedCharArray> nColors;
	vector<int> pTypes, nTypes;

	for (unsigned s = 0; s < Model::stepNum; ++s) {
		vtkSmartPointer<vtkUnstructuredGrid> pUGrid = vtkSmartPointer<vtkUnstructuredGrid>::New();
		vtkSmartPointer<vtkUnstructuredGrid> nUGrid = vtkSmartPointer<vtkUnstructuredGrid>::New();
		
		pUGrid->SetPoints(pfeMesh->getpvtkPnts(s));
		nUGrid->SetPoints(nfeMesh->getpvtkPnts(s));
		
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
		pUGrid->SetCells(pTypes.data(), pCells);
		pUGrid->GetCellData()->SetScalars(pColors);

		nUGrid->SetCells(nTypes.data(), nCells);
		nUGrid->GetCellData()->SetScalars(nColors);
		
		pMappers[s] = vtkSmartPointer<vtkDataSetMapper>::New();
		pMappers[s]->SetInputData(pUGrid);

		nMappers[s] = vtkSmartPointer<vtkDataSetMapper>::New();
		nMappers[s]->SetInputData(nUGrid);
	}

	// add them to pipeline
	pActor = vtkSmartPointer<vtkActor>::New();
	pActor->SetMapper(pMappers[0]);
	
	nActor = vtkSmartPointer<vtkActor>::New();
	nActor->SetMapper(nMappers[0]);

}

void ContactData::UpdateUGrid(unsigned s)
{
	pActor->SetMapper(pMappers[s]);
	nActor->SetMapper(nMappers[s]);
}
//////////////////////////////////////////////////////////////////////////////////////////////
//
// read contact files
//
//
//////////////////////////////////////////////////////////////////////////////////////////////


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
		
		for (unsigned s = 0; s < Model::stepNum; ++s) {
			//cout << "step: " << s << "\t eof?: " << infile.eof() << endl;
			if (!infile.eof()) {
				for (unsigned c = 0; c < pContacts.size(); ++c) {
					// primitives
					readpairs(infile, pContacts[c]->getPrimitivesPairs(NODE_TRIANGLE, s));
					readpairs(infile, pContacts[c]->getPrimitivesPairs(TRIANGLE_NODE, s));
					readpairs(infile, pContacts[c]->getPrimitivesPairs(EDGE_EDGE, s));
					readpairs(infile, pContacts[c]->getPrimitivesPairs(NODE_EDGE, s));
					readpairs(infile, pContacts[c]->getPrimitivesPairs(EDGE_NODE, s));
					readpairs(infile, pContacts[c]->getPrimitivesPairs(NODE_NODE, s));
				} // loop for ContactData
			}
		} // loop for steps

		infile.close();
	}
	else {
		cout << "Can not open contfile." << endl;
		exit(0);
	}

}

void readNodeDatafile(const string& nodedatafile, vector<shared_ptr<Model>> &pModels)
{
	fstream infile;
	infile.open(nodedatafile, ios::in | ios::binary);
	
	if (infile.is_open()) {
		double time, min, max;
		vector<vector<double>> dataCollect;
		vector<double> databuffer(Model::nodeNums);

		while (1) {
			infile.read((char*)&time, sizeof(double));
			infile.read((char*)databuffer.data(), sizeof(double)* Model::nodeNums);

			if (infile.fail()) break;

			dataCollect.push_back(databuffer);
			auto rg = std::minmax_element(databuffer.begin(), databuffer.end());
			min = *rg.first;
			max = *rg.second;
			Model::ranges.push_back(make_pair(min, max));

		}
		infile.close();
		// make sure the time step
		if (dataCollect.size() == Model::stepNum) {

			for (auto &pModel : pModels) {
				pModel->getNodelScalars().resize(Model::stepNum);
			}

			for (unsigned s = 0; s < Model::stepNum; ++s) {
				unsigned nodeoffset = 0;

				for (auto &pModel : pModels) {
					pModel->getNodelScalars(s) = vtkSmartPointer<vtkDoubleArray>::New();
					pModel->getNodelScalars(s)->SetNumberOfValues(pModel->getNodenum());
	
					for (unsigned n = 0; n < pModel->getNodenum(); ++n) {
						pModel->getNodelScalars(s)->SetValue(n, dataCollect[s][nodeoffset + n]);
					}
					
					nodeoffset += pModel->getNodenum();
				}
			}

			for (auto &pModel : pModels) {
				pModel->getFEMesh()->getUGrid()->GetPointData()->SetScalars(pModel->getNodelScalars(0));
			}
			
		}
		
	}
	else {
		cout << "Can not open NodeDatafile." << endl;
		exit(0);
	}


}



#endif //MODEL_H
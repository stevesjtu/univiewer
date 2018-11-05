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

  virtual std::vector<unsigned int> GetUintArrayFrom(const std::string &title) {
    JumpTo(title);
    std::vector<unsigned int> int_list;
    std::string temp_str;
    while (infile >> temp_str) {
      if (temp_str[0] == '#') break;
      int_list.push_back((unsigned int)std::stoi(temp_str));
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

  void CreateModel(const vector<unsigned int> &modelinfo,
                  const vector<unsigned int> &elemlist, const vector<double> &nodelist);

	void setLabelnode();
	void updateDisp(unsigned);
};







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

//////////////////////////////////////////////////////////////////////////////////////////////
//
// read contact files
//
//
//////////////////////////////////////////////////////////////////////////////////////////////
void readpairs(ifstream &infile, pairCollect &pc);

void readContfile(const string& contfile, vector<shared_ptr<ContactData>> &pContacts, vector<shared_ptr<Model>> &pModels);

void readNodeDatafile(const string& nodedatafile, vector<shared_ptr<Model>> &pModels);

#endif //MODEL_H
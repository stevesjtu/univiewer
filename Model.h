#pragma once
#ifndef MODEL_H
#define MODEL_H

#include "paramdefine.h"
#include "auxfunc.h"
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

// label
#include "vtkLabeledDataMapper.h"

// actor hold by model
#include "vtkActor.h"
#include "vtkActor2D.h"


namespace univiewer {

////////////////////////////////////////////////////////////////////////////
// feMesh
class FEMesh
{
private:
	bool isCalcEdge;
	vtkSmartPointer<vtkUnstructuredGrid> ugrid;
	std::vector<vtkSmartPointer<vtkPoints> > pvtkPnts;
	std::vector<std::array<unsigned, 2> > edges_bynodes;
	std::vector<std::array<unsigned, 3> > elems_bynodes;
public:
	FEMesh():isCalcEdge(false) {};
	virtual ~FEMesh() {};

	inline vtkSmartPointer<vtkUnstructuredGrid> &getUGrid() { return ugrid; }
	std::vector<vtkSmartPointer<vtkPoints> > &getpvtkPnts() { return pvtkPnts; }
	vtkSmartPointer<vtkPoints> &getpvtkPnts(unsigned i) { return pvtkPnts[i]; }
	std::vector<std::array<unsigned, 2> > &getEdges() { return edges_bynodes; }
	std::vector<std::array<unsigned, 3> > &getElems() { return elems_bynodes; }
	std::array<unsigned, 2> &getEdges(unsigned i) { return edges_bynodes[i]; }
	std::array<unsigned, 3> &getElems(unsigned i) { return elems_bynodes[i]; }

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
	vtkSmartPointer<vtkActor2D> labelActor;

	vtkSmartPointer<vtkDataSetMapper> mapper;

	sptr<FEMesh> feMesh;
	//shared_ptr<Labelnode> labelnode;

	std::vector<vtkSmartPointer<vtkDoubleArray> > nodalscalars;
	
public:
	static unsigned count, nodeNums, elemNums, stepNum, step;
	static std::vector<double> stepCollection;
	static std::vector<std::pair<double, double> > ranges;

	virtual ~Model() {};
	Model()
	{
		ID = count;
		++count;
	}

	inline unsigned getNodenum() { return nodeNum; }
	inline unsigned getOffset() { return offset; }
	inline vtkSmartPointer<vtkActor2D> getLabelactor() { return labelActor; }
	inline vtkSmartPointer<vtkActor> getActor() { return actor; }
	inline vtkSmartPointer<vtkDataSetMapper> getMapper() { return mapper; }
	
	inline sptr<FEMesh> getFEMesh() { return feMesh; }
	inline std::vector<vtkSmartPointer<vtkDoubleArray> > &getNodelScalars() { return nodalscalars; }
	inline vtkSmartPointer<vtkDoubleArray> &getNodelScalars(int s) { return nodalscalars[s]; }

	void setOffset(unsigned index) { offset = index; }
	void ReadXmlModel(const std::string&);
  void ReadTxtModel(const std::string&);

  void CreateModel(const std::vector<unsigned int> &modelinfo,
                  const std::vector<unsigned int> &elemlist, const std::vector<double> &nodelist);

	void setLabelnode();
	
	void updateDisp(unsigned);
};







////////////////////////////////////////////////////////////////////////////

typedef std::vector<std::pair<unsigned, unsigned> > pairCollect;
// ConstactInfo
class ContactData
{
private:
	sptr<Model> pModel, nModel;

	std::vector<pairCollect> node_triangles;
	std::vector<pairCollect> triangle_nodes;
	std::vector<pairCollect> edge_edges;
	std::vector<pairCollect> node_edges;
	std::vector<pairCollect> edge_nodes;
	std::vector<pairCollect> node_nodes;
	
	std::vector<vtkSmartPointer<vtkDataSetMapper> > pMappers;
	std::vector<vtkSmartPointer<vtkDataSetMapper> > nMappers;
	
	vtkSmartPointer<vtkActor> pActor;
	vtkSmartPointer<vtkActor> nActor;

	inline void insertNode(	
		unsigned n, 
		vtkSmartPointer<vtkCellArray> &Cells,
		vtkSmartPointer<vtkUnsignedCharArray> &Colors,
		std::vector<int> &Types);
	inline void insertEdge(
		std::array<unsigned, 2> edge,
		vtkSmartPointer<vtkCellArray> &Cells,
		vtkSmartPointer<vtkUnsignedCharArray> &Colors,
		std::vector<int> &Types);
	inline void insertTriangle(
		std::array<unsigned, 3> elem,
		vtkSmartPointer<vtkCellArray> &Cells,
		vtkSmartPointer<vtkUnsignedCharArray> &Colors,
		std::vector<int> &Types);

public:
	ContactData() {};
	virtual ~ContactData() {};

	void setModels(sptr<Model>& p, sptr<Model> &n) { pModel = p; nModel = n; }
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

	std::vector<pairCollect> &getPrimitivesPairs(int type)
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
void readpairs(std::ifstream &infile, pairCollect &pc);

void readContfile(const std::string& contfile, std::vector<sptr<ContactData> > &pContacts, std::vector<sptr<Model> > &pModels);

void readNodeDatafile(const std::string& nodedatafile, std::vector<sptr<Model> > &pModels);

}

#endif //MODEL_H
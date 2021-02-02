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

// actor_ hold by model
#include "vtkActor.h"
#include "vtkActor2D.h"


namespace univiewer {

////////////////////////////////////////////////////////////////////////////
// fe_mesh_
class FEMesh
{
private:
	bool is_calc_edge;
	vtkSmartPointer<vtkUnstructuredGrid> ugrid_;
	std::vector<vtkSmartPointer<vtkPoints> > pvtk_pnts_;
	std::vector<std::array<Uint, 2> > edges_bynodes_;
	std::vector<std::array<Uint, 3> > elems_bynodes_;
public:
	FEMesh():is_calc_edge(false) {};
	virtual ~FEMesh() {};

	inline vtkSmartPointer<vtkUnstructuredGrid> &GetUGrid() { return ugrid_; }
	std::vector<vtkSmartPointer<vtkPoints> > &GetPvtkpnts() { return pvtk_pnts_; }
	vtkSmartPointer<vtkPoints> &GetPvtkpnts(Uint i) { return pvtk_pnts_[i]; }
	std::vector<std::array<Uint, 2> > &GetEdges() { return edges_bynodes_; }
	std::vector<std::array<Uint, 3> > &GetElems() { return elems_bynodes_; }
	std::array<Uint, 2> &GetEdges(Uint i) { return edges_bynodes_[i]; }
	std::array<Uint, 3> &GetElems(Uint i) { return elems_bynodes_[i]; }

	bool getIsCalcEdge() { return is_calc_edge; }
	void MakeEdges();
};



////////////////////////////////////////////////////////////////////////////
// Model
class Model // base class
{
protected:
	Uint id_, num_node_, num_elem_, offset_;
	
	vtkSmartPointer<vtkActor> actor_;
	vtkSmartPointer<vtkActor2D> label_actor_;

	vtkSmartPointer<vtkDataSetMapper> mapper_;

	sptr<FEMesh> fe_mesh_;

	std::vector<vtkSmartPointer<vtkDoubleArray> > nodal_scalars_;
	
public:
	static Uint count_, num_nodes_, num_elems_, num_step_, step_;
	static std::vector<double> step_collection_;
	static std::vector<std::pair<double, double> > ranges_;

	virtual ~Model() {};
	Model()
	{
		id_ = count_;
		++count_;
	}

	inline Uint GetNumOfNode() { return num_node_; }
	inline Uint GetOffset() { return offset_; }
	inline vtkSmartPointer<vtkActor2D> GetLabelActor() { return label_actor_; }
	inline vtkSmartPointer<vtkActor> GetActor() { return actor_; }
	inline vtkSmartPointer<vtkDataSetMapper> GetMapper() { return mapper_; }
	
	inline sptr<FEMesh> GetFEMesh() { return fe_mesh_; }
	inline std::vector<vtkSmartPointer<vtkDoubleArray> > &GetNodelScalars() { return nodal_scalars_; }
	inline vtkSmartPointer<vtkDoubleArray> &GetNodelScalars(Uint s) { return nodal_scalars_[s]; }

	void SetOffset(Uint index) { offset_ = index; }
	void ReadXmlModel(const std::string&);
  void ReadTxtModel(const std::string&);

  void CreateModel(const std::vector<Uint> &modelinfo,
                  const std::vector<Uint> &elemlist, const std::vector<double> &nodelist);

	void SetLabelnode();
	
	void UpdateDisp(Uint);
};







////////////////////////////////////////////////////////////////////////////

typedef std::vector<std::pair<Uint, Uint> > PairCollect;
// ConstactInfo
class ContactData
{
private:
	sptr<Model> prev_model_, next_model_;

	std::vector<PairCollect> node_triangles_;
	std::vector<PairCollect> triangle_nodes_;
	std::vector<PairCollect> edge_edges_;
	std::vector<PairCollect> node_edges_;
	std::vector<PairCollect> edge_nodes_;
	std::vector<PairCollect> node_nodes_;
	
	std::vector<vtkSmartPointer<vtkDataSetMapper> > prev_mappers_;
	std::vector<vtkSmartPointer<vtkDataSetMapper> > next_mappers_;
	
	vtkSmartPointer<vtkActor> prev_actor_;
	vtkSmartPointer<vtkActor> next_actor_;

	inline void InsertNode(	
		Uint n, 
		vtkSmartPointer<vtkCellArray> &Cells,
		vtkSmartPointer<vtkUnsignedCharArray> &Colors,
		std::vector<int> &Types);
	inline void InsertEdge(
		std::array<Uint, 2> edge,
		vtkSmartPointer<vtkCellArray> &Cells,
		vtkSmartPointer<vtkUnsignedCharArray> &Colors,
		std::vector<int> &Types);
	inline void InsertTriangle(
		std::array<Uint, 3> elem,
		vtkSmartPointer<vtkCellArray> &Cells,
		vtkSmartPointer<vtkUnsignedCharArray> &Colors,
		std::vector<int> &Types);

public:
	ContactData() {};
	virtual ~ContactData() {};

	void SetModels(sptr<Model>& p, sptr<Model> &n) { prev_model_ = p; next_model_ = n; }
	PairCollect &GetPrimitivesPairs(int type, Uint step_)
	{
		switch (type)
		{
		case NODE_TRIANGLE:
			return node_triangles_[step_];
		case TRIANGLE_NODE:
			return triangle_nodes_[step_];
		case EDGE_EDGE:
			return edge_edges_[step_];
		case NODE_EDGE:
			return node_edges_[step_];
		case EDGE_NODE:
			return edge_nodes_[step_];
		case NODE_NODE:
			return node_nodes_[step_];
		default:
			break;
		}
		return node_triangles_[step_];
	}

	std::vector<PairCollect> &GetPrimitivesPairs(int type)
	{
		switch (type)
		{
		case NODE_TRIANGLE:
			return node_triangles_;
		case TRIANGLE_NODE:
			return triangle_nodes_;
		case EDGE_EDGE:
			return edge_edges_;
		case NODE_EDGE:
			return node_edges_;
		case EDGE_NODE:
			return edge_nodes_;
		case NODE_NODE:
			return node_nodes_;
		default:
			break;
		}
		return node_triangles_;
	}

	void InitializeUGrid();
	void UpdateUGrid(Uint s);
	vtkSmartPointer<vtkActor> &GetPrevActor() { return prev_actor_; }
	vtkSmartPointer<vtkActor> &GetNextActor() { return next_actor_; }
};

//////////////////////////////////////////////////////////////////////////////////////////////
//
// read contact files
//
//
//////////////////////////////////////////////////////////////////////////////////////////////
void ReadPairs(std::ifstream &infile, PairCollect &pc);

void ReadContactFile(const std::string& contfile, std::vector<sptr<ContactData> > &contacts_, std::vector<sptr<Model> > &models_);

void ReadNodalDataFile(const std::string& nodedatafile, std::vector<sptr<Model> > &models_);

}

#endif //MODEL_H
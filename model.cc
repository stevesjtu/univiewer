
#include"model.h"

namespace univiewer {

void FEMesh::MakeEdges()
{
	// edgeSet should be a set of edges without duplicated edges.
	Uint elem_num = ugrid_->GetNumberOfCells();
	elems_bynodes_.resize(elem_num);

	vtkSmartPointer<vtkCellArray> cellarray = ugrid_->GetCells();
	vtkSmartPointer<vtkIdList> idlist = vtkSmartPointer<vtkIdList>::New();
	for (Uint e = 0; e < elem_num; ++e) {
		cellarray->GetNextCell(idlist);
		elems_bynodes_[e] = { (Uint)idlist->GetId(0), (Uint)idlist->GetId(1), (Uint)idlist->GetId(2) };
	}
	
	typedef std::set<Uint> edge;
	edge oneEdge;
	std::set<edge> edgeSet;
	Uint ind1;
	for (Uint e = 0; e < elem_num; ++e) {
		for (Uint ind = 0; ind < 3; ++ind) {
			ind == 3 - 1 ? ind1 = 0 : ind1 = ind + 1;
			oneEdge.insert(elems_bynodes_[e][ind]);
			oneEdge.insert(elems_bynodes_[e][ind1]);
			edgeSet.insert(oneEdge);
			oneEdge.clear();
		}
	}

	//// fill the edges_bynodes_
	Uint edge_num = static_cast<Uint>(edgeSet.size());
	edges_bynodes_.resize(edge_num);
	Uint index = 0;
	for (std::set<edge>::iterator it = edgeSet.begin(); it != edgeSet.end(); it++) {
		edge::iterator nodeit = it->begin();
		edges_bynodes_[index][0] = *nodeit++;
		edges_bynodes_[index++][1] = *nodeit;
	}

	is_calc_edge = true;
}


Uint Model::count_ = 0;
Uint Model::num_nodes_ = 0;
Uint Model::num_elems_ = 0;
Uint Model::num_step_ = 0;
Uint Model::step_ = 0;
std::vector<double> Model::step_collection_ = std::vector<double>(0);
std::vector<std::pair<double, double> > Model::ranges_ = std::vector<std::pair<double, double> >(0);

void Model::CreateModel(const std::vector<Uint> &modelinfo,
  const std::vector<Uint> &elemlist, const std::vector<double> &nodelist) {
  
  fe_mesh_ = CreateOneOf<FEMesh>();
  fe_mesh_->GetUGrid() = vtkSmartPointer<vtkUnstructuredGrid>::New();

  vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();
  for (Uint i = 0; i < modelinfo[1]; ++i) {
    double nodes[3] = { nodelist[i * 3], nodelist[i * 3 + 1], nodelist[i * 3 + 2] };
    points->InsertPoint(i, nodes);
  }
  fe_mesh_->GetUGrid()->SetPoints(points);

  int vtktype;
  switch (modelinfo[2]) {
  case 3:
    vtktype = VTK_TRIANGLE;
    break;
  case 8:
    vtktype = VTK_HEXAHEDRON;
    break;
  case 2:
    vtktype = VTK_LINE;
    break;
  case 1:
    vtktype = VTK_VERTEX;
    break;
  default:
    break;
  }
  
  for (Uint i = 0; i < modelinfo[0]; ++i) {
    std::vector<vtkIdType> element(modelinfo[2]);
    for (Uint j = 0; j < modelinfo[2]; ++j) {
      element[j] = (vtkIdType)(*(elemlist.data() + i * modelinfo[2] + j) );
    }

    fe_mesh_->GetUGrid()->InsertNextCell(vtktype, (vtkIdType)modelinfo[2], element.data());

  }

  mapper_ = vtkSmartPointer<vtkDataSetMapper>::New();
  mapper_->SetInputData(fe_mesh_->GetUGrid());

  actor_ = vtkSmartPointer<vtkActor>::New();
  actor_->SetMapper(mapper_);

  // modify the static variables
  num_nodes_ += num_node_ = fe_mesh_->GetUGrid()->GetNumberOfPoints();
  num_elems_ += num_elem_ = fe_mesh_->GetUGrid()->GetNumberOfCells();
  
}

void Model::ReadXmlModel(const std::string&file)
{

	fe_mesh_ = CreateOneOf<FEMesh>();

	// read from xml
	vtkSmartPointer<vtkXMLUnstructuredGridReader> ugridReader = vtkSmartPointer<vtkXMLUnstructuredGridReader>::New();
	ugridReader->SetFileName(file.c_str());
	ugridReader->Update();
	// 
	fe_mesh_->GetUGrid() = vtkSmartPointer<vtkUnstructuredGrid>::New();
	fe_mesh_->GetUGrid() = ugridReader->GetOutput();

	mapper_ = vtkSmartPointer<vtkDataSetMapper>::New();
	mapper_->SetInputData(fe_mesh_->GetUGrid());
	
	actor_ = vtkSmartPointer<vtkActor>::New();
	actor_->SetMapper(mapper_);

	// modify the static variables
	num_nodes_ += num_node_ = fe_mesh_->GetUGrid()->GetNumberOfPoints();
	num_elems_ += num_elem_ = fe_mesh_->GetUGrid()->GetNumberOfCells();
}

void Model::ReadTxtModel(const std::string& file) {
  fe_mesh_ = CreateOneOf<FEMesh>();

  std::vector<Uint> modelinfo;
  std::vector<Uint> elemlist;
  std::vector<double> nodelist;

  if (file.substr(file.size() - 3).compare(".md") == 0) {
	  Mdfile txtfile(file);
	  modelinfo = txtfile.GetUintArrayFrom("elem_node_nodeOfElem");
	  elemlist = txtfile.GetUintArrayFrom("element_list");
	  nodelist = txtfile.GetDoubleArrayFrom("node_coordinate_list");
		//for(Uint &el : elemlist) el--;
	  txtfile.Close();
  } else if (file.substr(file.size() - 3).compare("txt") == 0 ){
    std::ifstream txtfile(file);
    if (!txtfile.is_open()) std::cout << "Error in open file: " << file << std::endl;
    Uint num_elem, num_node, nofe;
    txtfile >> num_elem >> num_node >> nofe;
    modelinfo.push_back(num_elem);
    modelinfo.push_back(num_node);
    modelinfo.push_back(nofe);

    elemlist.resize(nofe* num_elem);
    nodelist.resize(3 * num_node);

    for (Uint eln = 0; eln < num_elem* nofe; ++eln) {
      txtfile >> elemlist[eln];
      //elemlist[eln]--;
    }

    for (Uint ndn = 0; ndn < 3 * num_node; ++ndn) {
      txtfile >> nodelist[ndn];
    }
    txtfile.close();
  }

  this->CreateModel(modelinfo, elemlist, nodelist);

}


void Model::SetLabelnode()
{
	// labelnode = Labelnode::New();
	// labelnode->setLabelActor(fe_mesh_->GetUGrid());
	
	vtkSmartPointer<vtkLabeledDataMapper> labelMapper = vtkSmartPointer<vtkLabeledDataMapper>::New();
	labelMapper->SetInputData(fe_mesh_->GetUGrid());
	label_actor_ = vtkSmartPointer<vtkActor2D>::New();
	label_actor_->SetMapper(labelMapper);

}

void Model::UpdateDisp(Uint s)
{
	fe_mesh_->GetUGrid()->SetPoints(fe_mesh_->GetPvtkpnts(s));

	if (!nodal_scalars_.empty()) {
		fe_mesh_->GetUGrid()->GetPointData()->SetScalars(nodal_scalars_[s]);
		mapper_->SetScalarRange(ranges_[s].first, ranges_[s].second);
	}

}



void ContactData::InsertNode(
	Uint n,
	vtkSmartPointer<vtkCellArray> &Cells,
	vtkSmartPointer<vtkUnsignedCharArray> &Colors,
	std::vector<int> &Types)
{
	vtkSmartPointer<vtkVertex> node = vtkSmartPointer<vtkVertex>::New();
	node->GetPointIds()->SetId(0, n);
	// cell color type
	Cells->InsertNextCell(node);
	Colors->InsertNextTypedTuple(red);
	Types.push_back(VTK_VERTEX);
}

void ContactData::InsertEdge(
	std::array<Uint, 2> edge,
	vtkSmartPointer<vtkCellArray> &Cells,
	vtkSmartPointer<vtkUnsignedCharArray> &Colors,
	std::vector<int> &Types)
{
	vtkSmartPointer<vtkLine> line = vtkSmartPointer<vtkLine>::New();
	line->GetPointIds()->SetId(0, edge[0]);
	line->GetPointIds()->SetId(1, edge[1]);
	// cell color type
	Cells->InsertNextCell(line);
	Colors->InsertNextTypedTuple(blue);
	Types.push_back(VTK_LINE);
}

void ContactData::InsertTriangle(
	std::array<Uint, 3> elem,
	vtkSmartPointer<vtkCellArray> &Cells,
	vtkSmartPointer<vtkUnsignedCharArray> &Colors,
	std::vector<int> &Types)
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
	if (!prev_model_->GetFEMesh()->getIsCalcEdge()) prev_model_->GetFEMesh()->MakeEdges();
	if (!next_model_->GetFEMesh()->getIsCalcEdge()) next_model_->GetFEMesh()->MakeEdges();
		
	const auto &pfeMesh = prev_model_->GetFEMesh();
	const auto &nfeMesh = next_model_->GetFEMesh();

	prev_mappers_.resize(Model::num_step_);
	next_mappers_.resize(Model::num_step_);

	vtkSmartPointer<vtkCellArray> pCells;
	vtkSmartPointer<vtkCellArray> nCells;
	vtkSmartPointer<vtkUnsignedCharArray> pColors;
	vtkSmartPointer<vtkUnsignedCharArray> nColors;
	std::vector<int> pTypes, nTypes;

	for (Uint s = 0; s < Model::num_step_; ++s) {
		vtkSmartPointer<vtkUnstructuredGrid> pUGrid = vtkSmartPointer<vtkUnstructuredGrid>::New();
		vtkSmartPointer<vtkUnstructuredGrid> nUGrid = vtkSmartPointer<vtkUnstructuredGrid>::New();
		
		pUGrid->SetPoints(pfeMesh->GetPvtkpnts(s));
		nUGrid->SetPoints(nfeMesh->GetPvtkpnts(s));
		
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
		Uint n0, n1, e0, e1, t0, t1;
		for (Uint i = 0; i < node_triangles_[s].size(); ++i) {
			n0 = node_triangles_[s][i].first;
			t1 = node_triangles_[s][i].second;
			InsertNode(n0, pCells, pColors, pTypes);
			InsertTriangle(nfeMesh->GetElems(t1), nCells, nColors, nTypes);
		}

		for (Uint i = 0; i < triangle_nodes_[s].size(); ++i) {
			t0 = triangle_nodes_[s][i].first;
			n1 = triangle_nodes_[s][i].second;
			InsertTriangle(pfeMesh->GetElems(t0), pCells, pColors, pTypes);
			InsertNode(n1, nCells, nColors, nTypes);
		}

		for (Uint i = 0; i < edge_edges_[s].size(); ++i) {
			e0 = edge_edges_[s][i].first;
			e1 = edge_edges_[s][i].second;
			InsertEdge(pfeMesh->GetEdges(e0), pCells, pColors, pTypes);
			InsertEdge(nfeMesh->GetEdges(e1), nCells, nColors, nTypes);
		}

		for (Uint i = 0; i < node_edges_[s].size(); ++i) {
			n0 = node_edges_[s][i].first;
			e1 = node_edges_[s][i].second;
			InsertNode(n0, pCells, pColors, pTypes);
			InsertEdge(nfeMesh->GetEdges(e1), nCells, nColors, nTypes);
		}

		for (Uint i = 0; i < edge_nodes_[s].size(); ++i) {
			e0 = edge_nodes_[s][i].first;
			n1 = edge_nodes_[s][i].second;
			InsertEdge(pfeMesh->GetEdges(e0), pCells, pColors, pTypes);
			InsertNode(n1, nCells, nColors, nTypes);
		}

		for (Uint i = 0; i < node_nodes_[s].size(); ++i) {
			n0 = node_nodes_[s][i].first;
			n1 = node_nodes_[s][i].second;
			InsertNode(n0, pCells, pColors, pTypes);
			InsertNode(n1, nCells, nColors, nTypes);
		}

		// set the UnstructuredGrid
		pUGrid->SetCells(pTypes.data(), pCells);
		pUGrid->GetCellData()->SetScalars(pColors);

		nUGrid->SetCells(nTypes.data(), nCells);
		nUGrid->GetCellData()->SetScalars(nColors);
		
		prev_mappers_[s] = vtkSmartPointer<vtkDataSetMapper>::New();
		prev_mappers_[s]->SetInputData(pUGrid);

		next_mappers_[s] = vtkSmartPointer<vtkDataSetMapper>::New();
		next_mappers_[s]->SetInputData(nUGrid);
	}

	// add them to pipeline
	prev_actor_ = vtkSmartPointer<vtkActor>::New();
	prev_actor_->SetMapper(prev_mappers_[0]);
	
	next_actor_ = vtkSmartPointer<vtkActor>::New();
	next_actor_->SetMapper(next_mappers_[0]);

}

void ContactData::UpdateUGrid(Uint s)
{
	prev_actor_->SetMapper(prev_mappers_[s]);
	next_actor_->SetMapper(next_mappers_[s]);
}



//////////////////////////////////////////////////////////////////////////////////////////////
//
// read contact files
//
//
//////////////////////////////////////////////////////////////////////////////////////////////
void ReadPairs(std::ifstream &infile, PairCollect &pc)
{
	Uint num, pairs[2];
	infile.read((char*)&num, sizeof(Uint));
	for (Uint i = 0; i < num; ++i) {
		infile.read((char*)pairs, sizeof(Uint) * 2);
		pc.push_back(std::make_pair(pairs[0], pairs[1]));
	}
}

void ReadContactFile(const std::string& contfile, std::vector<sptr<ContactData> > &contacts_, std::vector<sptr<Model> > &models_)
{
	std::ifstream infile;
	infile.open(contfile, std::ios::in | std::ios::binary);
	if (infile.is_open()) {
		contacts_.resize(Model::num_step_);
		Uint bodyid1, bodyid2, mark;

		infile.read((char*)&mark, sizeof(Uint));
		contacts_.resize(mark);
		for (Uint c = 0; c < mark; ++c) {
			contacts_[c] = CreateOneOf<ContactData>();
			infile.read((char*)&bodyid1, sizeof(Uint));
			infile.read((char*)&bodyid2, sizeof(Uint));
			contacts_[c]->SetModels(models_[bodyid1], models_[bodyid2]);
			contacts_[c]->GetPrimitivesPairs(NODE_TRIANGLE).resize(Model::num_step_);
			contacts_[c]->GetPrimitivesPairs(TRIANGLE_NODE).resize(Model::num_step_);
			contacts_[c]->GetPrimitivesPairs(EDGE_EDGE).resize(Model::num_step_);
			contacts_[c]->GetPrimitivesPairs(NODE_EDGE).resize(Model::num_step_);
			contacts_[c]->GetPrimitivesPairs(EDGE_NODE).resize(Model::num_step_);
			contacts_[c]->GetPrimitivesPairs(NODE_NODE).resize(Model::num_step_);
		}
		
		for (Uint s = 0; s < Model::num_step_; ++s) {
			//cout << "step_: " << s << "\t eof?: " << infile.eof() << endl;
			if (!infile.eof()) {
				for (Uint c = 0; c < contacts_.size(); ++c) {
					// primitives
					ReadPairs(infile, contacts_[c]->GetPrimitivesPairs(NODE_TRIANGLE, s));
					ReadPairs(infile, contacts_[c]->GetPrimitivesPairs(TRIANGLE_NODE, s));
					ReadPairs(infile, contacts_[c]->GetPrimitivesPairs(EDGE_EDGE, s));
					ReadPairs(infile, contacts_[c]->GetPrimitivesPairs(NODE_EDGE, s));
					ReadPairs(infile, contacts_[c]->GetPrimitivesPairs(EDGE_NODE, s));
					ReadPairs(infile, contacts_[c]->GetPrimitivesPairs(NODE_NODE, s));
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

void ReadNodalDataFile(const std::string& nodedatafile, std::vector<sptr<Model> > &models_)
{
	std::fstream infile;
	infile.open(nodedatafile, ios::in | ios::binary);
	
	if (infile.is_open()) {
		double time, min, max;
		std::vector<std::vector<double> > dataCollect;
		std::vector<double> databuffer(Model::num_nodes_);

		while (1) {
			infile.read((char*)&time, sizeof(double));
			infile.read((char*)databuffer.data(), sizeof(double)* Model::num_nodes_);

			if (infile.fail()) break;

			dataCollect.push_back(databuffer);
			auto rg = std::minmax_element(databuffer.begin(), databuffer.end());
			min = *rg.first;
			max = *rg.second;
			Model::ranges_.push_back(std::make_pair(min, max));

		}
		infile.close();
		// make sure the time step_
		if (dataCollect.size() == Model::num_step_) {

			for (auto &prev_model_ : models_) {
				prev_model_->GetNodelScalars().resize(Model::num_step_);
			}

			for (Uint s = 0; s < Model::num_step_; ++s) {
				Uint nodeoffset = 0;

				for (auto &prev_model_ : models_) {
					prev_model_->GetNodelScalars(s) = vtkSmartPointer<vtkDoubleArray>::New();
					prev_model_->GetNodelScalars(s)->SetNumberOfValues(prev_model_->GetNumOfNode());
	
					for (Uint n = 0; n < prev_model_->GetNumOfNode(); ++n) {
						prev_model_->GetNodelScalars(s)->SetValue(n, dataCollect[s][nodeoffset + n]);
					}
					
					nodeoffset += prev_model_->GetNumOfNode();
				}
			}

			for (auto &prev_model_ : models_) {
				prev_model_->GetFEMesh()->GetUGrid()->GetPointData()->SetScalars(prev_model_->GetNodelScalars(0));
			}
			
		}
		
	}
	else {
		cout << "Can not open NodeDatafile." << endl;
		exit(0);
	}


}

}

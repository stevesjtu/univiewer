
#include"model.h"

namespace univiewer {

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


unsigned Model::count = 0;
unsigned Model::nodeNums = 0;
unsigned Model::elemNums = 0;
unsigned Model::stepNum = 0;
unsigned Model::step = 0;
vector<double> Model::stepCollection = vector<double>(0);
vector<pair<double, double> > Model::ranges = vector<pair<double, double> >(0);

void Model::CreateModel(const vector<unsigned int> &modelinfo,
  const vector<unsigned int> &elemlist, const vector<double> &nodelist) {
  
  feMesh = FEMesh::New();
  feMesh->getUGrid() = vtkSmartPointer<vtkUnstructuredGrid>::New();

  vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();
  for (unsigned i = 0; i < modelinfo[1]; ++i) {
    double nodes[3] = { nodelist[i * 3], nodelist[i * 3 + 1], nodelist[i * 3 + 2] };
    points->InsertPoint(i, nodes);
  }
  feMesh->getUGrid()->SetPoints(points);

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
  
  for (unsigned i = 0; i < modelinfo[0]; ++i) {
    std::vector<vtkIdType> element(modelinfo[2]);
    for (unsigned j = 0; j < modelinfo[2]; ++j) {
      element[j] = (vtkIdType)(*(elemlist.data() + i * modelinfo[2] + j) );
    }

    feMesh->getUGrid()->InsertNextCell(vtktype, (vtkIdType)modelinfo[2], element.data());

  }

  mapper = vtkSmartPointer<vtkDataSetMapper>::New();
  mapper->SetInputData(feMesh->getUGrid());

  actor = vtkSmartPointer<vtkActor>::New();
  actor->SetMapper(mapper);

  // modify the static variables
  nodeNums += nodeNum = feMesh->getUGrid()->GetNumberOfPoints();
  elemNums += elemNum = feMesh->getUGrid()->GetNumberOfCells();
  
}

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

  std::vector<unsigned int> modelinfo;
  std::vector<unsigned int> elemlist;
  std::vector<double> nodelist;

  if (file.substr(file.size() - 3).compare(".md") == 0) {
	  Mdfile txtfile(file);
	  modelinfo = txtfile.GetUintArrayFrom("model info");
	  elemlist = txtfile.GetUintArrayFrom("element list");
	  nodelist = txtfile.GetDoubleArrayFrom("node list");
		for(unsigned &el : elemlist) el--;
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
	// labelnode = Labelnode::New();
	// labelnode->setLabelActor(feMesh->getUGrid());
	
	vtkSmartPointer<vtkLabeledDataMapper> labelMapper = vtkSmartPointer<vtkLabeledDataMapper>::New();
	labelMapper->SetInputData(feMesh->getUGrid());
	labelActor = vtkSmartPointer<vtkActor2D>::New();
	labelActor->SetMapper(labelMapper);

}

void Model::updateDisp(unsigned s)
{
	feMesh->getUGrid()->SetPoints(feMesh->getpvtkPnts(s));

	if (!nodalscalars.empty()) {
		feMesh->getUGrid()->GetPointData()->SetScalars(nodalscalars[s]);
		mapper->SetScalarRange(ranges[s].first, ranges[s].second);
	}

}



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
void readpairs(ifstream &infile, pairCollect &pc)
{
	unsigned num, pairs[2];
	infile.read((char*)&num, sizeof(unsigned int));
	for (unsigned i = 0; i < num; ++i) {
		infile.read((char*)pairs, sizeof(unsigned int) * 2);
		pc.push_back(make_pair(pairs[0], pairs[1]));
	}
}

void readContfile(const string& contfile, vector<shared_ptr<ContactData> > &pContacts, vector<shared_ptr<Model> > &pModels)
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

void readNodeDatafile(const string& nodedatafile, vector<shared_ptr<Model> > &pModels)
{
	fstream infile;
	infile.open(nodedatafile, ios::in | ios::binary);
	
	if (infile.is_open()) {
		double time, min, max;
		vector<vector<double> > dataCollect;
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

}

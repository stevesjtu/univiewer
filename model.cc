
#include"model.h"

void Model::CreateModel(const vector<int> &modelinfo,
  const vector<int> &elemlist, const vector<double> &nodelist) {
  
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
      element[j] = (vtkIdType)(*(elemlist.data() + i * modelinfo[2] + j) - 1);
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



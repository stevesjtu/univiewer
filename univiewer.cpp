#include "ControlView.h"
int main(int argc, char *argv[])
{

	shared_ptr<ControlView> pControlView = ControlView::New();
	pControlView->setRender();

	pControlView->inputModelfiles(argc, argv);

	pControlView->setAnimationMethod(DEFAULT_TIMERCALLBACK);
	pControlView->setKeyboardMethod(DEFAULT_KEYPRESSCALLBACK);
	pControlView->setWindowMethod(DEFAULT_WINDOWCALLBACK);
	pControlView->Display();
	
	return EXIT_SUCCESS;

#ifdef HOW_TO_WRITE_UNSTRUCTUREDGRID
	vtkSmartPointer<vtkUnstructuredGrid> ug = vtkSmartPointer<vtkUnstructuredGrid>::New();

	// Create points
	double origin[3] = { 0.0, 0.0, 0.0 };
	double p0[3] = { 1.0, 0.0, 0.0 };
	double p1[3] = { 0.0, 1.0, 0.0 };

	double pnode[3] = { 1.0, 1.0, 0.0 };

	double ptri0[3] = { 0.5, 0.5, 0.0 };
	double ptri1[3] = { 0.1, 0.3, 0.0 };
	double ptri2[3] = { 0.7, 0.1, 0.0 };

	vtkSmartPointer<vtkPoints> pts =
		vtkSmartPointer<vtkPoints>::New();
	pts->InsertNextPoint(origin);
	pts->InsertNextPoint(p0);
	pts->InsertNextPoint(p1);
	pts->InsertNextPoint(pnode);

	pts->InsertNextPoint(ptri0);
	pts->InsertNextPoint(ptri1);
	pts->InsertNextPoint(ptri2);

	ug->SetPoints(pts);

	//////////////////////////////////////////////////////
	vtkSmartPointer<vtkLine> line0 =
		vtkSmartPointer<vtkLine>::New();
	line0->GetPointIds()->SetId(0, 0);
	line0->GetPointIds()->SetId(1, 1);

	vtkSmartPointer<vtkLine> line1 =
		vtkSmartPointer<vtkLine>::New();
	line1->GetPointIds()->SetId(0, 0);
	line1->GetPointIds()->SetId(1, 2);

	vtkSmartPointer<vtkCellArray> cells =
		vtkSmartPointer<vtkCellArray>::New();
	cells->InsertNextCell(line0);
	cells->InsertNextCell(line1);
	////////////////////////////////////////////////////

	vtkSmartPointer<vtkVertex> vertex =
		vtkSmartPointer<vtkVertex>::New();
	vertex->GetPointIds()->SetId(0, 3);

	cells->InsertNextCell(vertex);
	/////////////////////////////////////////////////////
	vtkSmartPointer<vtkTriangle> tri =
		vtkSmartPointer<vtkTriangle>::New();
	tri->GetPointIds()->SetId(0, 4);
	tri->GetPointIds()->SetId(1, 5);
	tri->GetPointIds()->SetId(2, 6);

	cells->InsertNextCell(tri);

	/////////////////////////////////////////////////////
	unsigned char red[3] = { 255, 0, 0 };
	unsigned char green[3] = { 0, 255, 0 };
	unsigned char blue[3] = { 0, 0, 255 };
	unsigned char cyan[3] = { 0, 255, 255 };
	vtkSmartPointer<vtkUnsignedCharArray> colors =
		vtkSmartPointer<vtkUnsignedCharArray>::New();
	colors->SetNumberOfComponents(3);
	colors->InsertNextTupleValue(red);
	colors->InsertNextTupleValue(green);
	colors->InsertNextTupleValue(blue);
	colors->InsertNextTupleValue(cyan);

	int type[4] = { VTK_LINE, VTK_LINE, VTK_VERTEX, VTK_TRIANGLE };

	ug->SetCells(type, cells);

	ug->GetCellData()->SetScalars(colors);

	// Setup the visualization pipeline
	//vtkSmartPointer<vtkPolyDataMapper> mapper =
	//	vtkSmartPointer<vtkPolyDataMapper>::New();
	vtkSmartPointer<vtkDataSetMapper> mapper = vtkSmartPointer<vtkDataSetMapper>::New();

#if VTK_MAJOR_VERSION <= 5
	mapper->SetInput(linesPolyData);
#else
	mapper->SetInputData(ug);
#endif

	vtkSmartPointer<vtkActor> actor =
		vtkSmartPointer<vtkActor>::New();
	actor->SetMapper(mapper);
	
	actor->GetProperty()->SetLineWidth(2);
	actor->GetProperty()->SetPointSize(10);

	actor->GetProperty()->SetEdgeColor(1.0, 0.0, 0.0);
	actor->GetProperty()->EdgeVisibilityOff();
	
	vtkSmartPointer<vtkRenderer> renderer =
		vtkSmartPointer<vtkRenderer>::New();
	renderer->AddActor(actor);

	vtkSmartPointer<vtkRenderWindow> window =
		vtkSmartPointer<vtkRenderWindow>::New();
	window->AddRenderer(renderer);

	vtkSmartPointer<vtkRenderWindowInteractor> interactor =
		vtkSmartPointer<vtkRenderWindowInteractor>::New();
	interactor->SetRenderWindow(window);
	vtkSmartPointer<vtkInteractorStyleTrackballCamera> style =
		vtkSmartPointer<vtkInteractorStyleTrackballCamera>::New();
	interactor->SetInteractorStyle(style);

	// Visualize
	window->Render();
	interactor->Start();

	return EXIT_SUCCESS;
#endif
}



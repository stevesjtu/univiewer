#pragma once
#ifndef CONTROLVIEW_H
#define CONTROLVIEW_H
#include "vtkXMLUnstructuredGridReader.h"
#include "vtkUnstructuredGrid.h"
#include "vtkAppendFilter.h"
#include "vtkGeometryFilter.h"
#include "vtkUnstructuredGridGeometryFilter.h"
#include "vtkSmartPointer.h"
#include "vtkSphereSource.h"
#include "vtkPolyDataMapper.h"
#include "vtkLabeledDataMapper.h"
#include "vtkPointData.h"
#include "vtkDataSetMapper.h"
#include "vtkActor.h"
#include "vtkProperty.h"
#include "vtkExtractEdges.h"
#include "vtkRenderWindow.h"
//#include "vtkCocoaRenderWindow.h"
#include "vtkRenderer.h"
#include "vtkRendererCollection.h"
#include "vtkActorCollection.h"
#include "vtkAxesActor.h"
#include "vtkOrientationMarkerWidget.h"
#include "vtkRenderWindowInteractor.h"
#include "vtkCamera.h"
#include "vtkProgrammableFilter.h"
#include "vtkCallbackCommand.h"
#include "vtkCommand.h"
#include "vtkInteractorStyleTrackballCamera.h"

#include "vtkTriangle.h"
#include "vtkCellArray.h"
#include "vtkCellData.h"
#include "vtkUnsignedCharArray.h"
#include "vtkPoints.h"
#include "vtkLine.h"
#include "vtkPolyData.h"
#include "vtkPolyDataMapper.h"
#include "vtkImageData.h"
#include "vtkCoordinate.h"
#include "vtkButtonWidget.h"
#include "vtkTexturedButtonRepresentation2D.h"

#include "vtkTextActor.h"
#include "vtkTextProperty.h"
#include "vtkTextRepresentation.h"

#include "vtkProperty2D.h"
#include "vtkWidgetEvent.h"
#include "vtkWidgetEventTranslator.h"
#include "vtkSliderRepresentation2D.h"
#include "vtkSliderWidget.h"

#include <iostream>
#include <vector>
#include <string>
#include <fstream>
#include <memory>
#include "Eigen/eigen"
#include <vtkPoints.h>

using namespace Eigen;
using namespace std;

#define DEFAULT_TIMERCALLBACK TimerCallback
#define DEFAULT_KEYPRESSCALLBACK KeypressCallbackFunction
#define DEFAULT_WINDOWCALLBACK WindowModifiedCallback

// Timer callback proxy
void TimerCallbackFunction(vtkObject* caller, long unsigned int vtkNotUsed(eventId), void* clientData, void* vtkNotUsed(callData));

// real callback function
void TimerCallback(void* arguments);
void WindowModifiedCallback(vtkObject* caller, long unsigned int eventId, void* clientData, void* callData);
void KeypressCallbackFunction(vtkObject* caller, long unsigned int eventId, void* clientData, void* callData);

struct TriangleMesh
{
	TriangleMesh() {};
	TriangleMesh(MatrixXu &elem, MatrixXd &node) {
		pelem = &elem;
		pnode = &node;
	}
	MatrixXu *pelem;
	MatrixXd *pnode;
};

class ControlView
{
protected:
	vtkSmartPointer<vtkActor> actor;
	vtkSmartPointer<vtkActor> axesActor;
	vtkSmartPointer<vtkTextActor> textActor;
    vtkSmartPointer<vtkActor2D> labelActor;

    vector<vtkSmartPointer<vtkXMLUnstructuredGridReader> > ugridReaders;
    vtkSmartPointer<vtkLabeledDataMapper> labelMapper;
	vtkSmartPointer<vtkDataSetMapper> mapper;

	vtkSmartPointer<vtkRenderWindow> renderWindow;
	vtkSmartPointer<vtkRenderer> renderer;
	vtkSmartPointer<vtkRenderWindowInteractor> renderWindowInteractor;
	vtkSmartPointer<vtkInteractorStyleTrackballCamera> style;
	vtkSmartPointer<vtkAxesActor> axes;
	vtkSmartPointer<vtkOrientationMarkerWidget> widget;
	vtkSmartPointer<vtkPolyData> linesPolyData;
	vtkSmartPointer<vtkPoints> pts;
	vtkSmartPointer<vtkLine> line0, line1, line2;
	vtkSmartPointer<vtkCellArray> lines;
	vtkSmartPointer<vtkUnsignedCharArray> colors;
	vtkSmartPointer<vtkPolyDataMapper> axesMapper;

	//vtkSmartPointer<vtkAppendFilter> appendFilter;
	vtkSmartPointer<vtkProgrammableFilter> programmableFilter;
	vtkSmartPointer<vtkCallbackCommand> timerCallback;
	vtkSmartPointer<vtkCallbackCommand> keypressCallback;
	vtkSmartPointer<vtkCallbackCommand> windowCallback;

	//////////////////////////////////
	//////////////////////////////////
	unsigned stepNum;
	unsigned step;
	bool play;
	bool ShowMarker;
	bool ShowMesh;
    bool ShowLabel;

	virtual void inputModel_unit(MatrixXu &elem, MatrixXd &node,
					vtkSmartPointer<vtkPoints> &points, vtkSmartPointer<vtkCellArray> &cellArray)
	{
		unsigned nodeNum = node.cols();
		;
		for (unsigned i = 0; i< nodeNum; ++i) {
			points->InsertNextPoint(node.col(i).x(), node.col(i).y(), node.col(i).z());
		}

		unsigned elemNum = elem.cols();
		vector<vtkSmartPointer<vtkTriangle> > triangle(elemNum);
		

		for (unsigned e = 0; e< elemNum; ++e) {
			triangle[e] = vtkSmartPointer<vtkTriangle>::New();
			triangle[e]->GetPointIds()->SetId(0, elem(0, e));
			triangle[e]->GetPointIds()->SetId(1, elem(1, e));
			triangle[e]->GetPointIds()->SetId(2, elem(2, e));
			cellArray->InsertNextCell(triangle[e]);
		}
	}

public:
	virtual ~ControlView() {}
	ControlView(): stepNum(0), step(0),  play(false), ShowMarker(true), ShowMesh(true), ShowLabel(false) {};
	static shared_ptr<ControlView> New()
	{
		shared_ptr<ControlView> nw = make_shared<ControlView>();
		return nw;
	}

	unsigned &getStepNum() { return stepNum; }
	unsigned & getStep() { return step; }
	bool & IsPlay() { return play; }
	bool & IsShowMarker() { return ShowMarker; }
	bool & IsShowMesh() { return ShowMesh; }
    bool & IsShowLabel() { return ShowLabel; }

	vtkSmartPointer<vtkProgrammableFilter> & getProgrammableFilter() { return programmableFilter; }
	//vtkSmartPointer<vtkAppendFilter> &getAppendFilter() { return appendFilter; }
	vtkSmartPointer<vtkActor> &getActor() { return actor; }
	vtkSmartPointer<vtkTextActor> &getTextActor() {return textActor;}
	vtkSmartPointer<vtkActor> &getAxesActor() { return axesActor; }
    vtkSmartPointer<vtkActor2D> &getLabelActor() { return labelActor; }
    vtkSmartPointer<vtkRenderer> &getRenderer(){return renderer;}
    vtkSmartPointer<vtkRenderWindowInteractor> &getRenderWindowInteractor(){return renderWindowInteractor;}
    vtkSmartPointer<vtkCallbackCommand> &getTimerCallback(){return timerCallback;}

	///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	//input model for elem node structure
	///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    virtual void inputModel(MatrixXu &elem, MatrixXd &node)
    {

		vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();
		vtkSmartPointer<vtkCellArray> cellArray = vtkSmartPointer<vtkCellArray>::New();
		inputModel_unit(elem, node, points, cellArray);

		vtkSmartPointer<vtkUnstructuredGrid> unstructuredGrid = vtkSmartPointer<vtkUnstructuredGrid>::New();
        unstructuredGrid->SetPoints(points);
        unstructuredGrid->SetCells(VTK_TRIANGLE, cellArray);
        
        programmableFilter = vtkSmartPointer<vtkProgrammableFilter>::New();
		programmableFilter->AddInputData(unstructuredGrid);

    }

	///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	//input model for <vector> mesh structure
	///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	virtual void inputModel(vector<TriangleMesh> &mesh)
	{

		vtkSmartPointer<vtkAppendFilter> appendFilter = vtkSmartPointer<vtkAppendFilter>::New();
		vtkSmartPointer<vtkPoints> points;
		vtkSmartPointer<vtkCellArray> cellArray;

		for (unsigned i = 0; i < mesh.size(); ++i) {
			points = vtkSmartPointer<vtkPoints>::New();
			cellArray = vtkSmartPointer<vtkCellArray>::New();
			inputModel_unit(*mesh[i].pelem, *mesh[i].pnode, points, cellArray);
			vtkSmartPointer<vtkUnstructuredGrid> unstructuredGrid = vtkSmartPointer<vtkUnstructuredGrid>::New();
			unstructuredGrid->SetPoints(points);
			unstructuredGrid->SetCells(VTK_TRIANGLE, cellArray);
			appendFilter->AddInputData(unstructuredGrid);
		}

		programmableFilter = vtkSmartPointer<vtkProgrammableFilter>::New();
		programmableFilter->AddInputConnection(appendFilter->GetOutputPort());

	}

	///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	//input model for files
	///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	virtual void inputModelfiles(vector<string> &modelFiles, const int& argc,  char* argv[]) 
	{
		vector<string> *ptr_vector_name = NULL;
		for (int i = 1; i<argc; ++i) {

			if ((*argv[i] == '-') || (*argv[i] == '/')) {
				switch (*(argv[i] + 1)) {
				case 'm':
					ptr_vector_name = &modelFiles;
					break;
				default:
					break;
				}
				continue;
			}
			ptr_vector_name->push_back(argv[i]);
		}

		ugridReaders.resize(modelFiles.size());
		vtkSmartPointer<vtkAppendFilter> appendFilter = vtkSmartPointer<vtkAppendFilter>::New();
		for (unsigned i = 0; i< (unsigned)modelFiles.size(); ++i) {
			
			ugridReaders[i] = vtkSmartPointer<vtkXMLUnstructuredGridReader>::New();
			ugridReaders[i]->SetFileName(modelFiles[i].c_str());
			ugridReaders[i]->Update();
			appendFilter->AddInputConnection(ugridReaders[i]->GetOutputPort());
			//appendFilter->AddInputData(ugridReaders[i]->GetOutput());
		}
		appendFilter->Update();

		programmableFilter = vtkSmartPointer<vtkProgrammableFilter>::New();
		programmableFilter->AddInputConnection(appendFilter->GetOutputPort());
	}

	virtual void setKeyboardMethod(void (*f)(vtkObject* , long unsigned int vtkNotUsed(eventId), void* clientData, void* vtkNotUsed(callData))){
		keypressCallback = vtkSmartPointer<vtkCallbackCommand>::New();
		keypressCallback->SetCallback(f);
		keypressCallback->SetClientData(this);
		renderWindowInteractor->AddObserver(vtkCommand::KeyPressEvent, keypressCallback);
	}

	virtual void setWindowMethod(void (*f)(vtkObject* caller, long unsigned int vtkNotUsed(eventId), void* clientData, void* vtkNotUsed(callData))){
		windowCallback = vtkSmartPointer<vtkCallbackCommand>::New();
		windowCallback->SetCallback(f);
		windowCallback->SetClientData(this);
		renderWindow->AddObserver(vtkCommand::ModifiedEvent, windowCallback);
	}

	virtual void setMainActor() {
		// Create a mapper and actor
		mapper = vtkSmartPointer<vtkDataSetMapper>::New();
		mapper->SetInputConnection(programmableFilter->GetOutputPort());
		actor = vtkSmartPointer<vtkActor>::New();
		actor->SetMapper(mapper);
		actor->GetProperty()->SetColor(0.9, 0.9, 0.9);
		actor->GetProperty()->SetEdgeColor(0.0, 0.0, 0.0);
		actor->GetProperty()->EdgeVisibilityOn();
		//actor->SetScale(0.1);
	}

    virtual void setLabelActor(){
        labelMapper = vtkSmartPointer<vtkLabeledDataMapper>::New();

        for (unsigned i = 0; i< (unsigned)ugridReaders.size(); ++i) {
            labelMapper->AddInputConnection(ugridReaders[i]->GetOutputPort());
        }
        labelActor = vtkSmartPointer<vtkActor2D>::New();
        labelActor->SetMapper(labelMapper);

    }

	virtual void setRender() {
		// Create a renderer, render window, and interactor
		renderer = vtkSmartPointer<vtkRenderer>::New();
		renderer->GetActiveCamera()->SetClippingRange(0.01, 1000);
		renderer->SetAmbient(1.0, 1.0, 1.0);
		renderer->SetLightFollowCamera(1);

		renderWindow = vtkSmartPointer<vtkRenderWindow>::New();
		renderWindow->AddRenderer(renderer);
		renderWindow->SetSize(800, 640);
		renderWindow->SetWindowName("MLV 2.0");
        renderWindow->Render();

		renderWindowInteractor = vtkSmartPointer<vtkRenderWindowInteractor>::New();
		renderWindowInteractor->SetRenderWindow(renderWindow);
		style = vtkSmartPointer<vtkInteractorStyleTrackballCamera>::New();
		renderWindowInteractor->SetInteractorStyle(style);

		// Initialize must be called prior to creating timer events.
		renderWindowInteractor->Initialize();

		axes = vtkSmartPointer<vtkAxesActor>::New();
		widget = vtkSmartPointer<vtkOrientationMarkerWidget>::New();
		widget->SetOutlineColor(0.9300, 0.5700, 0.1300);
		widget->SetOrientationMarker(axes);
		widget->SetInteractor(renderWindowInteractor);
		widget->SetViewport(0.0, 0.0, 0.2, 0.2);
		widget->SetEnabled(1);
		widget->InteractiveOn();
		widget->SetInteractive(0);

		int* windowSize = renderWindow->GetSize();
		textActor->SetPosition(windowSize[0] - 170, windowSize[1] - 25);
	}

	virtual void setAxesActor() {
		linesPolyData = vtkSmartPointer<vtkPolyData>::New();

		// Create three points
		double px0[3] = { -1e2, 0.0, 0.0 };
		double px1[3] = { 1e2, 0.0, 0.0 };

		double py0[3] = { 0.0, -1e2, 0.0 };
		double py1[3] = { 0.0, 1e2, 0.0 };

		double pz0[3] = { 0.0, 0.0, -1e2 };
		double pz1[3] = { 0.0, 0.0, 1e2 };

		// Create a vtkPoints container and store the points in it
		pts = vtkSmartPointer<vtkPoints>::New();
		pts->InsertNextPoint(px0);
		pts->InsertNextPoint(px1);

		pts->InsertNextPoint(py0);
		pts->InsertNextPoint(py1);

		pts->InsertNextPoint(pz0);
		pts->InsertNextPoint(pz1);
		// Add the points to the polydata container
		linesPolyData->SetPoints(pts);

		// Create the first line (between Origin and P0)
		line0 = vtkSmartPointer<vtkLine>::New();
		line0->GetPointIds()->SetId(0, 0); // the second 0 is the index of the Origin in linesPolyData's points
		line0->GetPointIds()->SetId(1, 1); // the second 1 is the index of P0 in linesPolyData's points

										   // Create the second line (between Origin and P1)
		line1 = vtkSmartPointer<vtkLine>::New();
		line1->GetPointIds()->SetId(0, 2); // the second 0 is the index of the Origin in linesPolyData's points
		line1->GetPointIds()->SetId(1, 3); // 2 is the index of P1 in linesPolyData's points

	    line2 = vtkSmartPointer<vtkLine>::New();
		line2->GetPointIds()->SetId(0, 4); // the second 0 is the index of the Origin in linesPolyData's points
		line2->GetPointIds()->SetId(1, 5); // 3 is the index of P1 in linesPolyData's points


										   // Create a vtkCellArray container and store the lines in it
		lines = vtkSmartPointer<vtkCellArray>::New();
		lines->InsertNextCell(line0);
		lines->InsertNextCell(line1);
		lines->InsertNextCell(line2);

		// Add the lines to the polydata container
		linesPolyData->SetLines(lines);

		// Create two colors - one for each line
		unsigned char red[3] = { 255, 0, 0 };
		unsigned char green[3] = { 0, 255, 0 };
		unsigned char blue[3] = { 0, 0, 255 };
		// Create a vtkUnsignedCharArray container and store the colors in it
		colors = vtkSmartPointer<vtkUnsignedCharArray>::New();
		colors->SetNumberOfComponents(3);
		colors->InsertNextTupleValue(red);
		colors->InsertNextTupleValue(green);
		colors->InsertNextTupleValue(blue);
		// Color the lines.
		linesPolyData->GetCellData()->SetScalars(colors);

		// Setup the visualization pipeline
		axesMapper = vtkSmartPointer<vtkPolyDataMapper>::New();
		axesMapper->SetInputData(linesPolyData);
		axesActor = vtkSmartPointer<vtkActor>::New();
		axesActor->SetMapper(axesMapper);

	}

	virtual void setTextActor() {
		textActor = vtkSmartPointer<vtkTextActor>::New();
		//textActor->SetInput("Hello world");
		//textActor->SetPosition(80, 40);
		//textActor->GetTextProperty()->SetFontSize(24);
		//textActor->GetTextProperty()->SetColor(1.0, 0.0, 0.0);

		////////////////////

		textActor->GetTextProperty()->SetFontSize(16);
		textActor->GetTextProperty()->SetColor(1.0, 1.0, 1.0);

	}

	virtual void Display() {
		// Add the actor to the scene
		renderer->AddActor(actor);
		renderer->AddActor(axesActor);
		renderer->AddActor2D(textActor);
        renderer->AddActor2D(labelActor);
        labelActor->VisibilityOff();

		renderer->GradientBackgroundOn();
		renderer->SetBackground2(13.0 / 255.0, 71.0 / 255.0, 161.0 / 255.0);
		renderer->SetBackground(144.0 / 255.0, 202.0 / 255.0, 249.0 / 255.0);
		// Render and interact

		renderWindowInteractor->Start();
	}

};

void KeypressCallbackFunction(vtkObject* caller, long unsigned int vtkNotUsed(eventId), void* clientData, void* vtkNotUsed(callData))
{

	vtkRenderWindowInteractor *iren =
		static_cast<vtkRenderWindowInteractor*>(caller);

	ControlView* pCtr = static_cast<ControlView*>(clientData);

	if (strcmp(iren->GetKeySym(), "i") == 0) {
		//vtkRendererCollection* rendercollction = iren->GetRenderWindow()->GetRenderers();
		if (pCtr->IsShowMarker())
			//rendercollction->GetFirstRenderer()->RemoveActor(pCtr->getAxesActor());
            pCtr->getAxesActor()->VisibilityOff();
		else
//			rendercollction->GetFirstRenderer()->AddActor(pCtr->getAxesActor());
            pCtr->getAxesActor()->VisibilityOn();
		pCtr->IsShowMarker() = !pCtr->IsShowMarker();
	}

	if (strcmp(iren->GetKeySym(), "u") == 0) {
		if (pCtr->IsShowMesh())
			pCtr->getActor()->GetProperty()->EdgeVisibilityOff();
		else
			pCtr->getActor()->GetProperty()->EdgeVisibilityOn();

		pCtr->IsShowMesh() = !pCtr->IsShowMesh();
	}

	if (strcmp(iren->GetKeySym(), "l") == 0) {
        if (pCtr->IsShowLabel())
            pCtr->getLabelActor()->VisibilityOff();
        else
            pCtr->getLabelActor()->VisibilityOn();

        pCtr->IsShowLabel() = !pCtr->IsShowLabel();
	}

	iren->Render();
}
#endif //CONTROLVIEW_H
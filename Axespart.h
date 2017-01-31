#pragma once
#ifndef AXESPART_H
#define AXESPART_H
#include "vtkSmartPointer.h"
#include "vtkAxesActor.h"
#include "vtkOrientationMarkerWidget.h"
#include "vtkRenderWindowInteractor.h"

#include "vtkActor.h"
#include "vtkPolyData.h"
#include "vtkPolyDataMapper.h"
#include "vtkPoints.h"
#include "vtkLine.h"
#include "vtkCellArray.h"
#include "vtkCellData.h"
#include "vtkUnsignedCharArray.h"

#include<memory>
using namespace std;

class Axesline
{
private:

	// actor
	vtkSmartPointer<vtkActor> axesActor;
	vtkSmartPointer<vtkPolyDataMapper> axesMapper;

	vtkSmartPointer<vtkPolyData> linesPolyData;
	vtkSmartPointer<vtkPoints> pts;
	vtkSmartPointer<vtkLine> line0;
	vtkSmartPointer<vtkLine> line1;
	vtkSmartPointer<vtkLine> line2;
	vtkSmartPointer<vtkCellArray> lines;
	vtkSmartPointer<vtkUnsignedCharArray> colors;

public:
	Axesline() {};
	virtual ~Axesline() {};
	static shared_ptr<Axesline> New()
	{
		return make_shared<Axesline>();
	}

	vtkSmartPointer<vtkActor> &getAxesActor() { return axesActor; }

	void setAxesActor()
	{
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

};


class Axesframe
{
private:
	// widget
	vtkSmartPointer<vtkAxesActor> axes;
	vtkSmartPointer<vtkOrientationMarkerWidget> widget;

public:
	Axesframe() {};
	virtual ~Axesframe() {};
	static shared_ptr<Axesframe> New()
	{
		return make_shared<Axesframe>();
	}
	void setAxesWidget(vtkSmartPointer<vtkRenderWindowInteractor> &renderWindowInteractor)
	{

		axes = vtkSmartPointer<vtkAxesActor>::New();
		widget = vtkSmartPointer<vtkOrientationMarkerWidget>::New();
		widget->SetOutlineColor(0.9300, 0.5700, 0.1300);
		widget->SetOrientationMarker(axes);
		widget->SetInteractor(renderWindowInteractor);
		widget->SetViewport(0.0, 0.0, 0.2, 0.2);
		widget->EnabledOn();
		widget->InteractiveOff();
		
	}
};

#endif
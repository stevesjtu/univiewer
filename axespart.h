#pragma once
#ifndef AXESPART_H
#define AXESPART_H

#include "paramdefine.h"

#include "vtkAxesActor.h"
#include "vtkOrientationMarkerWidget.h"
#include "vtkRenderWindowInteractor.h"
#include "vtkPolyData.h"
#include "vtkPolyDataMapper.h"
#include "vtkCellData.h"
#include "vtkUnsignedCharArray.h"


namespace univiewer {

class Axesline
{
private:

	// actor_
	vtkSmartPointer<vtkActor> axes_actor_;

public:
	Axesline() {};
	virtual ~Axesline() {};

	vtkSmartPointer<vtkActor> &GetAxesActor() { return axes_actor_; }

	void SetAxesActor()
	{
		vtkSmartPointer<vtkPolyData> linesPolyData = vtkSmartPointer<vtkPolyData>::New();

		// Create three points
		double px0[3] = { -1e2, 0.0, 0.0 };
		double px1[3] = { 1e2, 0.0, 0.0 };

		double py0[3] = { 0.0, -1e2, 0.0 };
		double py1[3] = { 0.0, 1e2, 0.0 };

		double pz0[3] = { 0.0, 0.0, -1e2 };
		double pz1[3] = { 0.0, 0.0, 1e2 };

		// Create a vtkPoints container and store the points in it
		vtkSmartPointer<vtkPoints> pts = vtkSmartPointer<vtkPoints>::New();
		pts->InsertNextPoint(px0);
		pts->InsertNextPoint(px1);

		pts->InsertNextPoint(py0);
		pts->InsertNextPoint(py1);

		pts->InsertNextPoint(pz0);
		pts->InsertNextPoint(pz1);
		// Add the points to the polydata container
		linesPolyData->SetPoints(pts);

		// Create the first line (between Origin and P0)
		vtkSmartPointer<vtkLine> line0 = vtkSmartPointer<vtkLine>::New();
		line0->GetPointIds()->SetId(0, 0); // the second 0 is the index of the Origin in linesPolyData's points
		line0->GetPointIds()->SetId(1, 1); // the second 1 is the index of P0 in linesPolyData's points

										   // Create the second line (between Origin and P1)
		vtkSmartPointer<vtkLine> line1 = vtkSmartPointer<vtkLine>::New();
		line1->GetPointIds()->SetId(0, 2); // the second 0 is the index of the Origin in linesPolyData's points
		line1->GetPointIds()->SetId(1, 3); // 2 is the index of P1 in linesPolyData's points

		vtkSmartPointer<vtkLine> line2 = vtkSmartPointer<vtkLine>::New();
		line2->GetPointIds()->SetId(0, 4); // the second 0 is the index of the Origin in linesPolyData's points
		line2->GetPointIds()->SetId(1, 5); // 3 is the index of P1 in linesPolyData's points

		// Create a vtkCellArray container and store the lines in it
		vtkSmartPointer<vtkCellArray> lines = vtkSmartPointer<vtkCellArray>::New();
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
		vtkSmartPointer<vtkUnsignedCharArray> colors = vtkSmartPointer<vtkUnsignedCharArray>::New();
		colors->SetNumberOfComponents(3);
		colors->InsertNextTypedTuple(red);
		colors->InsertNextTypedTuple(green);
		colors->InsertNextTypedTuple(blue);
		// Color the lines.
		linesPolyData->GetCellData()->SetScalars(colors);
		
		// Setup the visualization pipeline
		vtkSmartPointer<vtkPolyDataMapper> axesMapper = vtkSmartPointer<vtkPolyDataMapper>::New();
		axesMapper->SetInputData(linesPolyData);
		axes_actor_ = vtkSmartPointer<vtkActor>::New();
		axes_actor_->SetMapper(axesMapper);

	}

};


class Axesframe
{
private:
	// widget_
	vtkSmartPointer<vtkOrientationMarkerWidget> widget_;

public:
	Axesframe() {};
	virtual ~Axesframe() {};
	
	void SetAxesWidget(vtkSmartPointer<vtkRenderWindowInteractor> &render_window_interactor_)
	{

		vtkSmartPointer<vtkAxesActor> axes = vtkSmartPointer<vtkAxesActor>::New();
		widget_ = vtkSmartPointer<vtkOrientationMarkerWidget>::New();
		widget_->SetOutlineColor(0.9300, 0.5700, 0.1300);
		widget_->SetOrientationMarker(axes);
		widget_->SetInteractor(render_window_interactor_);
		widget_->SetViewport(0.0, 0.0, 0.2, 0.2);
		widget_->EnabledOn();
		widget_->InteractiveOff();
	}
};

}

#endif
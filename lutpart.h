#pragma once
#ifndef LUTPART_H
#define LUTPART_H

#include "paramdefine.h"

#include "vtkDoubleArray.h"
#include "vtkLookupTable.h"
#include "vtkPointData.h"
#include "vtkScalarBarActor.h"

#include "model.h"

namespace univiewer {

class LookUpTable
{
private:

	vtkSmartPointer<vtkScalarBarActor> scalarBar;
	vtkSmartPointer<vtkLookupTable> hueLut;
public:
	LookUpTable() {};
	virtual ~LookUpTable() {};

	vtkSmartPointer<vtkScalarBarActor> &getScalarBar() { return scalarBar; }

	void setScalars( std::vector<sptr<Model> > pModels)
	{

		hueLut = vtkSmartPointer<vtkLookupTable>::New();
		hueLut->SetHueRange(1.85 / 3.0, 0);
		hueLut->Build();

		for (auto& pmodel : pModels) {
			pmodel->getMapper()->SetLookupTable(hueLut);
		}
		//mapper->SetScalarRange(0.0, 600.0);
		//mapper->ScalarVisibilityOn();
		//mapper->SetScalarModeToUsePointData();
		//mapper->SetColorModeToMapScalars();
		
		scalarBar = vtkSmartPointer<vtkScalarBarActor>::New();
		scalarBar->SetLookupTable(hueLut);
		scalarBar->SetTitle("Data");
		scalarBar->SetNumberOfLabels(8);

		scalarBar->GetPositionCoordinate()->SetCoordinateSystemToDisplay();
		scalarBar->GetPosition2Coordinate()->SetCoordinateSystemToDisplay();
		
		scalarBar->SetPosition(800 - 70 - 10, 10);
		scalarBar->SetPosition2(70, 300);
	}
};


}



#endif // !LUTPART_H

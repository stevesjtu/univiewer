#pragma once
#ifndef LUTPART_H
#define LUTPART_H

#include "ParamDefine.h"

#include "vtkDoubleArray.h"
#include "vtkLookupTable.h"
#include "vtkPointData.h"
#include "vtkScalarBarActor.h"

#include "Model.h"

class LookUpTable
{
private:

	vtkSmartPointer<vtkScalarBarActor> scalarBar;
	vtkSmartPointer<vtkLookupTable> hueLut;
public:
	LookUpTable() {};
	virtual ~LookUpTable() {};
	static shared_ptr<LookUpTable> New()
	{
		return make_shared<LookUpTable>();
	}

	vtkSmartPointer<vtkScalarBarActor> &getScalarBar() { return scalarBar; }

	void setScalars(vtkSmartPointer<vtkDataSetMapper> &mapper, vector<shared_ptr<Model>> &pModels)
	{

		hueLut = vtkSmartPointer<vtkLookupTable>::New();
		hueLut->SetTableRange(0, 1);
		hueLut->Build();

		mapper->SetScalarRange(0.0, 1.0);
		mapper->ScalarVisibilityOn();
		mapper->SetScalarModeToUsePointData();
		mapper->SetColorModeToMapScalars();
		mapper->SetLookupTable(hueLut);

		//scalarBar = vtkSmartPointer<vtkScalarBarActor>::New();
		//scalarBar->SetLookupTable(mapper->GetLookupTable());
		//scalarBar->SetTitle("Data");
		//scalarBar->SetNumberOfLabels(6);

		//scalarBar->GetPositionCoordinate()->SetCoordinateSystemToDisplay();
		//scalarBar->GetPosition2Coordinate()->SetCoordinateSystemToDisplay();
		//scalarBar->SetPosition(10, 10);
		//scalarBar->SetPosition2(50, 300);
		//
		//scalarBar->SetLookupTable(hueLut);
	}
};






#endif // !LUTPART_H

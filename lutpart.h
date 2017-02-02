#pragma once
#ifndef LUTPART_H
#define LUTPART_H

#include "ParamDefine.h"

#include "vtkDoubleArray.h"
#include "vtkLookupTable.h"
#include "vtkPointData.h"
#include "vtkScalarBarActor.h"

class LookUpTable
{
private:
	vtkSmartPointer<vtkDoubleArray> scalars;
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
	vtkSmartPointer<vtkDoubleArray> &getScalars() { return scalars; }
	void setScalars(vtkSmartPointer<vtkDataSetMapper> &mapper, const unsigned nodenum)
	{
		scalars = vtkSmartPointer<vtkDoubleArray>::New();
		scalars->SetNumberOfValues(nodenum);

		mapper->ScalarVisibilityOn();
		mapper->SetScalarModeToUsePointData();
		mapper->SetColorModeToMapScalars();
		mapper->Update();

		scalarBar = vtkSmartPointer<vtkScalarBarActor>::New();
		scalarBar->SetLookupTable(mapper->GetLookupTable());
		scalarBar->SetTitle("Title");
		scalarBar->SetNumberOfLabels(6);

		hueLut = vtkSmartPointer<vtkLookupTable>::New();
		hueLut->SetTableRange(0, 1);
		hueLut->SetHueRange(0, 1);
		hueLut->SetSaturationRange(1, 1);
		hueLut->SetValueRange(1, 1);
		hueLut->Build();

		mapper->SetLookupTable(hueLut);
		scalarBar->SetLookupTable(hueLut);
	}
};






#endif // !LUTPART_H

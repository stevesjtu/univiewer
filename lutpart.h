#pragma once
#ifndef LUTPART_H
#define LUTPART_H

#include "vtkDoubleArray.h"
#include "vtkLookupTable.h"
#include "vtkPointData.h"
#include "vtkScalarBarActor.h"
#include "vtkSmartPointer.h"
#include "vtkActor.h"
#include "vtkDataSetMapper.h"
#include "vtkUnstructuredGrid.h"
#include "vtkDataSetMapper.h"

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

	void setScalars(vtkSmartPointer<vtkDataSetMapper> &mapper, vtkUnstructuredGrid* usgrid)
	{
		
		scalars = vtkSmartPointer<vtkDoubleArray>::New();
		unsigned numPts = usgrid->GetNumberOfPoints();
		scalars->SetNumberOfValues(numPts);
		for (unsigned i = 0; i < numPts; ++i) {
			scalars->SetValue(i, static_cast<double>(i) / (double)numPts);
		}
		usgrid->GetPointData()->SetScalars(scalars);

		mapper->ScalarVisibilityOn();
		mapper->SetScalarModeToUsePointData();
		mapper->SetColorModeToMapScalars();

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

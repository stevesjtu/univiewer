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

	vtkSmartPointer<vtkScalarBarActor> scalar_bar_;
	vtkSmartPointer<vtkLookupTable> hue_lut_;
public:
	LookUpTable() {};
	virtual ~LookUpTable() {};

	vtkSmartPointer<vtkScalarBarActor> &GetScalarBar() { return scalar_bar_; }

	void SetScalars( std::vector<sptr<Model> > models_)
	{

		hue_lut_ = vtkSmartPointer<vtkLookupTable>::New();
		hue_lut_->SetHueRange(1.85 / 3.0, 0);
		hue_lut_->Build();

		for (auto& pmodel : models_) {
			pmodel->GetMapper()->SetLookupTable(hue_lut_);
		}
		//mapper_->SetScalarRange(0.0, 600.0);
		//mapper_->ScalarVisibilityOn();
		//mapper_->SetScalarModeToUsePointData();
		//mapper_->SetColorModeToMapScalars();
		
		scalar_bar_ = vtkSmartPointer<vtkScalarBarActor>::New();
		scalar_bar_->SetLookupTable(hue_lut_);
		scalar_bar_->SetTitle("Data");
		scalar_bar_->SetNumberOfLabels(8);

		scalar_bar_->GetPositionCoordinate()->SetCoordinateSystemToDisplay();
		scalar_bar_->GetPosition2Coordinate()->SetCoordinateSystemToDisplay();
		
		scalar_bar_->SetPosition(800 - 70 - 10, 10);
		scalar_bar_->SetPosition2(70, 300);
	}
};


}



#endif // !LUTPART_H

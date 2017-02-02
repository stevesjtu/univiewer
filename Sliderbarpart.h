#pragma once
#ifndef SLIDERBARPART_H
#define SLIDERBARPART_H

#include "ParamDefine.h"

#include "vtkSliderRepresentation2D.h"
#include "vtkSliderWidget.h"

#include "vtkImageData.h"
#include "vtkCoordinate.h"
#include "vtkTexturedButtonRepresentation2D.h"
#include "vtkButtonWidget.h"

#include "vtkProperty2D.h"
#include "vtkWidgetEvent.h"
#include "vtkWidgetEventTranslator.h"

#include "vtkCommand.h"

#include "vtkRenderer.h"
#include "vtkRenderWindowInteractor.h"
using namespace std;

void CreateImagePause(vtkSmartPointer<vtkImageData> image)
{
	// Specify the size of the image data
	image->SetDimensions(10, 10, 1);
	image->AllocateScalars(VTK_UNSIGNED_CHAR, 4);
	int* dims = image->GetDimensions();

	// Fill the image with
	for (int y = 0; y < dims[1]; y++)
	{
		for (int x = 0; x < dims[0]; x++)
		{
			unsigned char* pixel =
				static_cast<unsigned char*>(image->GetScalarPointer(x, y, 0));
			if (x < 4)
			{
				pixel[0] = 255;
				pixel[1] = 255;
				pixel[2] = 255;
				pixel[3] = 50;
			}
			else if (x > 5)
			{
				pixel[0] = 255;
				pixel[1] = 255;
				pixel[2] = 255;
				pixel[3] = 50;
			}
			else {
				pixel[3] = 0;
			}
		}
	}
}

void CreateImagePlay(vtkSmartPointer<vtkImageData> image)
{
	// Specify the size of the image data
	image->SetDimensions(200, 200, 1);
	image->AllocateScalars(VTK_UNSIGNED_CHAR, 4);
	int* dims = image->GetDimensions();

	// Fill the image with
	for (int x = 0; x < dims[0]; x++) {
		for (int y = 0; y< x / 2; y++) {
			unsigned char* pixel = static_cast<unsigned char*>(image->GetScalarPointer(x, y, 0));
			pixel[3] = 0;
		}
		for (int y = dims[1] - x / 2; y< dims[1]; y++) {
			unsigned char* pixel = static_cast<unsigned char*>(image->GetScalarPointer(x, y, 0));
			pixel[3] = 0;
		}
		for (int y = x / 2; y < dims[1] - x / 2; y++) {
			unsigned char* pixel = static_cast<unsigned char*>(image->GetScalarPointer(x, y, 0));
			pixel[0] = 255;
			pixel[1] = 255;
			pixel[2] = 255;
			pixel[3] = 50;
		}
	}
}



class vtkSliderCallback : public vtkCommand
{
protected:

public:

	static vtkSliderCallback *New()
	{
		return new vtkSliderCallback;
	}
	virtual void Execute(vtkObject *caller, unsigned long, void*)
	{
		vtkSliderWidget *sliderWidget = reinterpret_cast<vtkSliderWidget*>(caller);
		*pst = (unsigned)static_cast<vtkSliderRepresentation *>(sliderWidget->GetRepresentation())->GetValue();
	}
	unsigned *pst;
};


class Sliderbar
{
private:

	vtkSmartPointer<vtkSliderRepresentation2D> sliderRep;
	vtkSmartPointer<vtkSliderWidget> sliderWidget;

	vtkSmartPointer<vtkSliderCallback> SliderCallback;

	vtkSmartPointer<vtkTexturedButtonRepresentation2D> buttonRepresentation;
	vtkSmartPointer<vtkButtonWidget> buttonWidget;

public:

	Sliderbar() {};
	virtual ~Sliderbar() {};

	static shared_ptr<Sliderbar> New()
	{
		return make_shared<Sliderbar>();
	}

	vtkSmartPointer<vtkSliderRepresentation2D> &getSliderRep() { return sliderRep; }
	vtkSmartPointer<vtkTexturedButtonRepresentation2D> &getButtonRepresentation() { return buttonRepresentation; }
	void setSliderBar(vtkSmartPointer<vtkRenderer> &renderer, 
					  vtkSmartPointer<vtkRenderWindowInteractor> &renderWindowInteractor, 
					  unsigned &step,
					  unsigned &stepNum)
	{

		sliderRep = vtkSmartPointer<vtkSliderRepresentation2D>::New();

		sliderRep->SetMinimumValue(0);
		sliderRep->SetMaximumValue(stepNum - 1);
		sliderRep->SetValue(0);

		//sliderRep->SetTitleText("");
		sliderRep->SetShowSliderLabel(0);

		// Set color properties:
		// Change the color of the knob that slides
		sliderRep->GetSliderProperty()->SetColor(1.0, 1.0, 1.0);//red
		sliderRep->GetSliderProperty()->SetOpacity(0.8);
		// Change the color of the bar
		sliderRep->GetTubeProperty()->SetColor(1.0, 1.0, 1.0);
		sliderRep->GetTubeProperty()->SetOpacity(0.2);
		// Change the color of the text indicating what the slider controls
		//sliderRep->GetTitleProperty()->SetColor(1,0,0);//red
		//        sliderRep->GetTitleProperty()->SetBold(0);
		//        sliderRep->GetTitleProperty()->SetShadow(0);
		// Change the color of the text displaying the value
		//sliderRep->GetLabelProperty()->SetColor(1,0,0);//red
		// Change the color of the knob when the mouse is held on it
		sliderRep->GetSelectedProperty()->SetColor(0.0, 0.0, 0.0);
		sliderRep->GetSelectedProperty()->SetOpacity(0.9);
		// Change the color of the ends of the bar
		//sliderRep->GetCapProperty()->SetColor(1,1,0);//yellow

		sliderRep->GetPoint1Coordinate()->SetCoordinateSystemToDisplay();
		sliderRep->GetPoint2Coordinate()->SetCoordinateSystemToDisplay();
		sliderRep->GetPoint1Coordinate()->SetValue(800 - 400, 640 - 15);
		sliderRep->GetPoint2Coordinate()->SetValue(800 - 200, 640 - 15);

		sliderRep->SetEndCapLength(0);

		sliderWidget = vtkSmartPointer<vtkSliderWidget>::New();
		sliderWidget->SetInteractor(renderWindowInteractor);
		sliderWidget->SetRepresentation(sliderRep);
		sliderWidget->SetAnimationModeToAnimate();
		sliderWidget->EnabledOn();

		SliderCallback = vtkSmartPointer<vtkSliderCallback>::New();
		SliderCallback->pst = &step;
		sliderWidget->AddObserver(vtkCommand::InteractionEvent, SliderCallback);

		// button widget
		vtkSmartPointer<vtkImageData> imagePlay = vtkSmartPointer<vtkImageData>::New();
		vtkSmartPointer<vtkImageData> imagePause = vtkSmartPointer<vtkImageData>::New();
		CreateImagePlay(imagePlay);
		CreateImagePause(imagePause);
		buttonRepresentation = vtkSmartPointer<vtkTexturedButtonRepresentation2D>::New();
		buttonRepresentation->SetNumberOfStates(2);
		buttonRepresentation->SetButtonTexture(0, imagePlay);
		buttonRepresentation->SetButtonTexture(1, imagePause);

		buttonWidget = vtkSmartPointer<vtkButtonWidget>::New();
		buttonWidget->SetInteractor(renderWindowInteractor);
		buttonWidget->SetRepresentation(buttonRepresentation);

		vtkSmartPointer<vtkCoordinate> upperRight = vtkSmartPointer<vtkCoordinate>::New();
		upperRight->SetCoordinateSystemToNormalizedDisplay();
		//upperRight->SetCoordinateSystemToDisplay();
		upperRight->SetValue(1.0, 1.0);

		double bds[6];
		double sz = 20.0;

		bds[0] = upperRight->GetComputedDisplayValue(renderer)[0] - sz - 175;
		bds[1] = bds[0] + sz;
		bds[2] = upperRight->GetComputedDisplayValue(renderer)[1] - sz - 5;
		bds[3] = bds[2] + sz;
		bds[4] = bds[5] = 0.0;

		// Scale to 1, default is .5
		buttonRepresentation->SetPlaceFactor(1.0);
		buttonRepresentation->PlaceWidget(bds);

		buttonWidget->On();
		
	}
};



#endif
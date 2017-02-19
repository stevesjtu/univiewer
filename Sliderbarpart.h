#pragma once
#ifndef SLIDERBARPART_H
#define SLIDERBARPART_H

#include "ParamDefine.h"

#include "vtkSliderRepresentation2D.h"
#include "vtkSliderWidget.h"

#include "vtkImageData.h"
#include "vtkTexturedButtonRepresentation2D.h"
#include "vtkButtonWidget.h"

#include "vtkProperty2D.h"
#include "vtkWidgetEvent.h"
#include "vtkWidgetEventTranslator.h"

#include "vtkCommand.h"
#include "vtkCallbackCommand.h"

#include "vtkRenderWindow.h"
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

void CreateImageNextStep(vtkSmartPointer<vtkImageData> image)
{
	// Specify the size of the image data
	image->SetDimensions(200, 200, 1);
	image->AllocateScalars(VTK_UNSIGNED_CHAR, 4);
	int* dims = image->GetDimensions();

	for (int x = 0; x < 50; x++) {
		for (int y = 0; y < dims[1]; y++) {
			unsigned char* pixel = static_cast<unsigned char*>(image->GetScalarPointer(x, y, 0));
			pixel[0] = 255;
			pixel[1] = 255;
			pixel[2] = 255;
			pixel[3] = 50;
		}
	}
	for (int x = 50; x < 100; x++) {
		for (int y = 0; y < dims[1]; y++) {
			unsigned char* pixel = static_cast<unsigned char*>(image->GetScalarPointer(x, y, 0));
			pixel[3] = 0;
		}
	}
	// Fill the image with
	for (int x = 100; x < dims[0]; x++) {
		for (int y = 0; y < x - 100; y++) {
			unsigned char* pixel = static_cast<unsigned char*>(image->GetScalarPointer(x, y, 0));
			pixel[3] = 0;
		}
		for (int y = dims[1] + 100 - x; y< dims[1]; y++) {
			unsigned char* pixel = static_cast<unsigned char*>(image->GetScalarPointer(x, y, 0));
			pixel[3] = 0;
		}
		for (int y = x - 100; y < dims[1] + 100 - x; y++) {
			unsigned char* pixel = static_cast<unsigned char*>(image->GetScalarPointer(x, y, 0));
			pixel[0] = 255;
			pixel[1] = 255;
			pixel[2] = 255;
			pixel[3] = 50;
		}
	}
}

void CreateImagePrevStep(vtkSmartPointer<vtkImageData> image)
{
	// Specify the size of the image data
	image->SetDimensions(200, 200, 1);
	image->AllocateScalars(VTK_UNSIGNED_CHAR, 4);
	int* dims = image->GetDimensions();

	for (int x = 150; x < dims[0]; x++) {
		for (int y = 0; y < dims[1]; y++) {
			unsigned char* pixel = static_cast<unsigned char*>(image->GetScalarPointer(x, y, 0));
			pixel[0] = 255;
			pixel[1] = 255;
			pixel[2] = 255;
			pixel[3] = 50;
		}
	}

	for (int x = 100; x < 150; x++) {
		for (int y = 0; y < dims[1]; y++) {
			unsigned char* pixel = static_cast<unsigned char*>(image->GetScalarPointer(x, y, 0));
			pixel[3] = 0;
		}
	}
	// Fill the image with
	for (int x = 0; x < dims[0] / 2; x++) {
		for (int y = 0; y < -x + 100; y++) {
			unsigned char* pixel = static_cast<unsigned char*>(image->GetScalarPointer(x, y, 0));
			pixel[3] = 0;
		}
		for (int y = dims[1] - 100 + x; y< dims[1]; y++) {
			unsigned char* pixel = static_cast<unsigned char*>(image->GetScalarPointer(x, y, 0));
			pixel[3] = 0;
		}
		for (int y = -x + 100; y < dims[1] - 100 + x; y++) {
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

	static vtkSliderCallback *New(){return new vtkSliderCallback;}
	virtual void Execute(vtkObject *caller, unsigned long, void*)
	{
		vtkSliderWidget *sliderWidget = reinterpret_cast<vtkSliderWidget*>(caller);
		*step = (unsigned)static_cast<vtkSliderRepresentation *>(sliderWidget->GetRepresentation())->GetValue();
		if (!(*play))
			*stepPlay = true;
	}
	unsigned *step;
	bool *play, *stepPlay;
};

class NextPressedCallback : public vtkCommand
{
public:
	static NextPressedCallback *New(){return new NextPressedCallback;}
	NextPressedCallback() {}
	virtual void Execute(vtkObject *caller, unsigned long, void*)
	{
		*step = (*step == *stepNum - 1) ? 0 : *step + 1;
		*stepPlay = true;
	}
	unsigned *step, *stepNum;
	bool *stepPlay;

};


class PrevPressedCallback : public vtkCommand
{
public:
	static PrevPressedCallback *New() { return new PrevPressedCallback; }
	PrevPressedCallback() {}
	virtual void Execute(vtkObject *caller, unsigned long, void*)
	{
		*step = (*step == 0) ? 0 : *step - 1;
		*stepPlay = true;
	}
	unsigned *step;
	bool *stepPlay;
};

class playCallback : public vtkCommand
{
public:
	static playCallback *New() { return new playCallback; }
	playCallback() {}
	virtual void Execute(vtkObject *caller, unsigned long, void*)
	{
		vtkButtonWidget *thisButton = reinterpret_cast<vtkButtonWidget*>(caller);
		vtkTexturedButtonRepresentation2D* thisButtonRep = 
			reinterpret_cast<vtkTexturedButtonRepresentation2D*>(thisButton->GetRepresentation());
		
		if (thisButtonRep->GetState()==0)
			*play = false;
		if (thisButtonRep->GetState() == 1)
			*play = true;

		*stepPlay = false;
	}
	bool *play, *stepPlay;
};


class CommandButton
{
private:
	vtkSmartPointer<vtkTexturedButtonRepresentation2D> buttonRepresentation;
	vtkSmartPointer<vtkButtonWidget> buttonWidget;

public:
	CommandButton() {};
	virtual ~CommandButton() {};
	static shared_ptr<CommandButton> New() { return make_shared<CommandButton>(); }
	void setOneStatesButtonContent(vtkSmartPointer<vtkImageData> image0,
				vtkSmartPointer<vtkRenderWindowInteractor> renderWindowInteractor)
	{
		buttonRepresentation = vtkSmartPointer<vtkTexturedButtonRepresentation2D>::New();
		buttonRepresentation->SetNumberOfStates(1);
		buttonRepresentation->SetButtonTexture(0, image0);

		buttonWidget = vtkSmartPointer<vtkButtonWidget>::New();
		buttonWidget->SetInteractor(renderWindowInteractor);
		buttonWidget->SetRepresentation(buttonRepresentation);

		buttonWidget->On();
	}
	
	void setTwoStatesButtonContent(vtkSmartPointer<vtkImageData> image0, vtkSmartPointer<vtkImageData> image1,
						  vtkSmartPointer<vtkRenderWindowInteractor> renderWindowInteractor)
	{
		buttonRepresentation = vtkSmartPointer<vtkTexturedButtonRepresentation2D>::New();
		buttonRepresentation->SetNumberOfStates(2);
		buttonRepresentation->SetButtonTexture(0, image0);
		buttonRepresentation->SetButtonTexture(1, image1);

		buttonWidget = vtkSmartPointer<vtkButtonWidget>::New();
		buttonWidget->SetInteractor(renderWindowInteractor);
		buttonWidget->SetRepresentation(buttonRepresentation);

		buttonWidget->On();
	}

	void setButtonSizePosition(int width, int height, int bottomLeftX, int bottomLeftY)
	{
		double bds[6];
		bds[0] = bottomLeftX;
		bds[1] = bottomLeftX + width;
		bds[2] = bottomLeftY;
		bds[3] = bottomLeftY + height;
		bds[4] = bds[5] = 0.0;

		// Scale to 1, default is .5
		buttonRepresentation->SetPlaceFactor(1.0);
		buttonRepresentation->PlaceWidget(bds);
	}
	
	void setButtonPressCallBack(vtkSmartPointer<vtkCommand> callback)
	{
		buttonWidget->AddObserver(vtkCommand::StateChangedEvent, callback);
	}

	//void setButtonReleaseCallBack(vtkSmartPointer<vtkCommand> callback)
	//{
	//	buttonWidget->AddObserver(vtkCommand::EndInteractionEvent, callback);
	//}
	int getState() { return buttonRepresentation->GetState(); }
	void setState(int input) { buttonRepresentation->SetState(input); }

};


class Sliderbar
{
private:

	vtkSmartPointer<vtkSliderRepresentation2D> sliderRep;
	vtkSmartPointer<vtkSliderWidget> sliderWidget;

	vtkSmartPointer<vtkSliderCallback> SliderCallback;
	vtkSmartPointer<vtkRenderWindowInteractor> renderWindowInteractor;

	shared_ptr<CommandButton> playButton;
	shared_ptr<CommandButton> prevButton;
	shared_ptr<CommandButton> nextButton;


public:

	Sliderbar() {};
	virtual ~Sliderbar() {};

	static shared_ptr<Sliderbar> New()
	{
		return make_shared<Sliderbar>();
	}

	vtkSmartPointer<vtkSliderRepresentation2D> &getSliderRep() { return sliderRep; }
	shared_ptr<CommandButton> &getPlayButton() { return playButton; }
	shared_ptr<CommandButton> &getNextButton() { return nextButton; }
	shared_ptr<CommandButton> &getPrevButton() { return prevButton; }
	
	void setSliderBar(vtkSmartPointer<vtkRenderWindowInteractor> renderWindowInteractor, 
					  unsigned &step,
					  unsigned &stepNum,
					  bool &play, bool &stepPlay)
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
		int *winsize = renderWindowInteractor->GetRenderWindow()->GetSize();
		sliderRep->GetPoint1Coordinate()->SetCoordinateSystemToDisplay();
		sliderRep->GetPoint2Coordinate()->SetCoordinateSystemToDisplay();

		//sliderRep->GetPoint1Coordinate()->SetValue(winsize[0] - 260 - winsize[0] / 4.0, winsize[1] - 15);
		//sliderRep->GetPoint2Coordinate()->SetValue(winsize[0] - 260, winsize[1] - 15);

		sliderRep->SetEndCapLength(0);

		sliderWidget = vtkSmartPointer<vtkSliderWidget>::New();
		sliderWidget->SetInteractor(renderWindowInteractor);
		sliderWidget->SetRepresentation(sliderRep);
		sliderWidget->SetAnimationModeToAnimate();
		sliderWidget->EnabledOn();

		SliderCallback = vtkSmartPointer<vtkSliderCallback>::New();
		SliderCallback->step = &step;
		SliderCallback->play = &play;
		SliderCallback->stepPlay = &stepPlay;
		sliderWidget->AddObserver(vtkCommand::InteractionEvent, SliderCallback);

		//////////////////////////////////////////////////////////////////////////////////////////////
		// play button widget ////////////////////////////////////////////////////////////////////////
		vtkSmartPointer<vtkImageData> imagePlay = vtkSmartPointer<vtkImageData>::New();
		vtkSmartPointer<vtkImageData> imagePause = vtkSmartPointer<vtkImageData>::New();
		CreateImagePlay(imagePlay);
		CreateImagePause(imagePause);
		playButton = CommandButton::New();
		playButton->setTwoStatesButtonContent(imagePlay, imagePause, renderWindowInteractor);
		//int width = 20, height = 20;
		//playButton->setButtonSizePosition(width, height, winsize[0] - 235 - width, winsize[1] - 5 - height);
		vtkSmartPointer<playCallback> playButtonCallback = vtkSmartPointer<playCallback>::New();
		playButtonCallback->play = &play;
		playButtonCallback->stepPlay = &stepPlay;
		playButton->setButtonPressCallBack(playButtonCallback);
		//////////////////////////////////////////////////////////////////////////////////////////////
		// Single Step button widget /////////////////////////////////////////////////////////////////
		vtkSmartPointer<vtkImageData> imageNext = vtkSmartPointer<vtkImageData>::New();
		vtkSmartPointer<vtkImageData> imagePrev = vtkSmartPointer<vtkImageData>::New();
		CreateImageNextStep(imageNext);
		CreateImagePrevStep(imagePrev);
		nextButton = CommandButton::New();
		prevButton = CommandButton::New();
		nextButton->setOneStatesButtonContent(imageNext, renderWindowInteractor);
		prevButton->setOneStatesButtonContent(imagePrev, renderWindowInteractor);
	
		//nextButton->setButtonSizePosition(width, height, winsize[0] - 210 - width, winsize[1] - 5 - height);
		//prevButton->setButtonSizePosition(width, height, winsize[0] - 185 - width, winsize[1] - 5 - height);
		// call back
		vtkSmartPointer<NextPressedCallback> nextPressedCallback = vtkSmartPointer<NextPressedCallback>::New();
		nextPressedCallback->step = &step;
		nextPressedCallback->stepNum = &stepNum;
		nextPressedCallback->stepPlay = &stepPlay;

		nextButton->setButtonPressCallBack(nextPressedCallback);

		vtkSmartPointer<PrevPressedCallback> prevPressedCallback = vtkSmartPointer<PrevPressedCallback>::New();
		prevPressedCallback->step = &step;
		prevPressedCallback->stepPlay = &stepPlay;

		prevButton->setButtonPressCallBack(prevPressedCallback);
	}
};



#endif
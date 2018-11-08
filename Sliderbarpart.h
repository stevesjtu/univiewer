#pragma once
#ifndef SLIDERBARPART_H
#define SLIDERBARPART_H

#include "paramdefine.h"

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

namespace univiewer {

void CreateImagePause(vtkSmartPointer<vtkImageData> image);

void CreateImagePlay(vtkSmartPointer<vtkImageData> image);

void CreateImageNextStep(vtkSmartPointer<vtkImageData> image);

void CreateImagePrevStep(vtkSmartPointer<vtkImageData> image);


class vtkSliderCallback : public vtkCommand
{
protected:

public:

	VTKSubClass(vtkSliderCallback)

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
	VTKSubClass(NextPressedCallback)
	NextPressedCallback() {}
	virtual void Execute(vtkObject *caller, unsigned long, void*)
	{
		*step = (*step == *stepNum - 1) ? *stepNum - 1 : *step + 1;
		*stepPlay = true;
	}
	unsigned *step, *stepNum;
	bool *stepPlay;

};


class PrevPressedCallback : public vtkCommand
{
public:
	VTKSubClass(PrevPressedCallback)
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
	VTKSubClass(playCallback)
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

	sptr<CommandButton> playButton;
	sptr<CommandButton> prevButton;
	sptr<CommandButton> nextButton;


public:

	Sliderbar() {};
	virtual ~Sliderbar() {};

	vtkSmartPointer<vtkSliderRepresentation2D> &getSliderRep() { return sliderRep; }
	sptr<CommandButton> &getPlayButton() { return playButton; }
	sptr<CommandButton> &getNextButton() { return nextButton; }
	sptr<CommandButton> &getPrevButton() { return prevButton; }
	
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
		playButton = CreateOneOf<CommandButton>();
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
		nextButton = CreateOneOf<CommandButton>();
		prevButton = CreateOneOf<CommandButton>();
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

}

#endif
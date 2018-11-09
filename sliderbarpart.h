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
		vtkSliderWidget *slider_widget_ = reinterpret_cast<vtkSliderWidget*>(caller);
		*step_ = (unsigned)static_cast<vtkSliderRepresentation *>(slider_widget_->GetRepresentation())->GetValue();
		if (!(*play_))
			*step_play_ = true;
	}
	unsigned *step_;
	bool *play_, *step_play_;
};

class NextPressedCallback : public vtkCommand
{
public:
	VTKSubClass(NextPressedCallback)
	NextPressedCallback() {}
	virtual void Execute(vtkObject *caller, unsigned long, void*)
	{
		*step_ = (*step_ == *num_step_ - 1) ? *num_step_ - 1 : *step_ + 1;
		*step_play_ = true;
	}
	unsigned *step_, *num_step_;
	bool *step_play_;

};


class PrevPressedCallback : public vtkCommand
{
public:
	VTKSubClass(PrevPressedCallback)
	PrevPressedCallback() {}
	virtual void Execute(vtkObject *caller, unsigned long, void*)
	{
		*step_ = (*step_ == 0) ? 0 : *step_ - 1;
		*step_play_ = true;
	}
	unsigned *step_;
	bool *step_play_;
};

class PlayCallback : public vtkCommand
{
public:
	VTKSubClass(PlayCallback)
	PlayCallback() {}
	virtual void Execute(vtkObject *caller, unsigned long, void*)
	{
		vtkButtonWidget *thisButton = reinterpret_cast<vtkButtonWidget*>(caller);
		vtkTexturedButtonRepresentation2D* thisButtonRep = 
			reinterpret_cast<vtkTexturedButtonRepresentation2D*>(thisButton->GetRepresentation());
		
		if (thisButtonRep->GetState()==0)
			*play_ = false;
		if (thisButtonRep->GetState() == 1)
			*play_ = true;

		*step_play_ = false;
	}
	bool *play_, *step_play_;
};


class CommandButton
{
private:
	vtkSmartPointer<vtkTexturedButtonRepresentation2D> button_representation_;
	vtkSmartPointer<vtkButtonWidget> button_widget_;

public:
	CommandButton() {};
	virtual ~CommandButton() {};

	void SetOneStatesButtonContent(vtkSmartPointer<vtkImageData> image0,
				vtkSmartPointer<vtkRenderWindowInteractor> render_window_interactor_)
	{
		button_representation_ = vtkSmartPointer<vtkTexturedButtonRepresentation2D>::New();
		button_representation_->SetNumberOfStates(1);
		button_representation_->SetButtonTexture(0, image0);

		button_widget_ = vtkSmartPointer<vtkButtonWidget>::New();
		button_widget_->SetInteractor(render_window_interactor_);
		button_widget_->SetRepresentation(button_representation_);

		button_widget_->On();
	}
	
	void SetTwoStatesButtonContent(vtkSmartPointer<vtkImageData> image0, vtkSmartPointer<vtkImageData> image1,
						  vtkSmartPointer<vtkRenderWindowInteractor> render_window_interactor_)
	{
		button_representation_ = vtkSmartPointer<vtkTexturedButtonRepresentation2D>::New();
		button_representation_->SetNumberOfStates(2);
		button_representation_->SetButtonTexture(0, image0);
		button_representation_->SetButtonTexture(1, image1);

		button_widget_ = vtkSmartPointer<vtkButtonWidget>::New();
		button_widget_->SetInteractor(render_window_interactor_);
		button_widget_->SetRepresentation(button_representation_);

		button_widget_->On();
	}

	void SetButtonSizePosition(int width, int height, int bottomLeftX, int bottomLeftY)
	{
		double bds[6];
		bds[0] = bottomLeftX;
		bds[1] = bottomLeftX + width;
		bds[2] = bottomLeftY;
		bds[3] = bottomLeftY + height;
		bds[4] = bds[5] = 0.0;

		// Scale to 1, default is .5
		button_representation_->SetPlaceFactor(1.0);
		button_representation_->PlaceWidget(bds);
	}
	
	void SetButtonPressCallBack(vtkSmartPointer<vtkCommand> callback)
	{
		button_widget_->AddObserver(vtkCommand::StateChangedEvent, callback);
	}

	//void setButtonReleaseCallBack(vtkSmartPointer<vtkCommand> callback)
	//{
	//	button_widget_->AddObserver(vtkCommand::EndInteractionEvent, callback);
	//}
	int GetState() { return button_representation_->GetState(); }
	void SetState(int input) { button_representation_->SetState(input); }

};


class SliderBar
{
private:

	vtkSmartPointer<vtkSliderRepresentation2D> slider_rep_;
	vtkSmartPointer<vtkSliderWidget> slider_widget_;

	vtkSmartPointer<vtkSliderCallback> slider_callback_;
	vtkSmartPointer<vtkRenderWindowInteractor> render_window_interactor_;

	sptr<CommandButton> play_button_;
	sptr<CommandButton> prev_button_;
	sptr<CommandButton> next_button_;


public:

	SliderBar() {};
	virtual ~SliderBar() {};

	vtkSmartPointer<vtkSliderRepresentation2D> &GetSliderRep() { return slider_rep_; }
	sptr<CommandButton> &GetPlayButton() { return play_button_; }
	sptr<CommandButton> &GetNextButton() { return next_button_; }
	sptr<CommandButton> &GetPrevButton() { return prev_button_; }
	
	void SetSliderBar(vtkSmartPointer<vtkRenderWindowInteractor> render_window_interactor_, 
					  unsigned &step_,
					  unsigned &num_step_,
					  bool &play_, bool &step_play_)
	{
		
		slider_rep_ = vtkSmartPointer<vtkSliderRepresentation2D>::New();

		slider_rep_->SetMinimumValue(0);
		slider_rep_->SetMaximumValue(num_step_ - 1);
		slider_rep_->SetValue(0);

		//slider_rep_->SetTitleText("");
		slider_rep_->SetShowSliderLabel(0);

		// Set color properties:
		// Change the color of the knob that slides
		slider_rep_->GetSliderProperty()->SetColor(1.0, 1.0, 1.0);//red
		slider_rep_->GetSliderProperty()->SetOpacity(0.8);
		// Change the color of the bar
		slider_rep_->GetTubeProperty()->SetColor(1.0, 1.0, 1.0);
		slider_rep_->GetTubeProperty()->SetOpacity(0.2);
		// Change the color of the text indicating what the slider controls
		//slider_rep_->GetTitleProperty()->SetColor(1,0,0);//red
		//        slider_rep_->GetTitleProperty()->SetBold(0);
		//        slider_rep_->GetTitleProperty()->SetShadow(0);
		// Change the color of the text displaying the value
		//slider_rep_->GetLabelProperty()->SetColor(1,0,0);//red
		// Change the color of the knob when the mouse is held on it
		slider_rep_->GetSelectedProperty()->SetColor(0.0, 0.0, 0.0);
		slider_rep_->GetSelectedProperty()->SetOpacity(0.9);
		// Change the color of the ends of the bar
		//slider_rep_->GetCapProperty()->SetColor(1,1,0);//yellow
		int *winsize = render_window_interactor_->GetRenderWindow()->GetSize();
		slider_rep_->GetPoint1Coordinate()->SetCoordinateSystemToDisplay();
		slider_rep_->GetPoint2Coordinate()->SetCoordinateSystemToDisplay();

		//slider_rep_->GetPoint1Coordinate()->SetValue(winsize[0] - 260 - winsize[0] / 4.0, winsize[1] - 15);
		//slider_rep_->GetPoint2Coordinate()->SetValue(winsize[0] - 260, winsize[1] - 15);

		slider_rep_->SetEndCapLength(0);

		slider_widget_ = vtkSmartPointer<vtkSliderWidget>::New();
		slider_widget_->SetInteractor(render_window_interactor_);
		slider_widget_->SetRepresentation(slider_rep_);
		slider_widget_->SetAnimationModeToAnimate();
		slider_widget_->EnabledOn();

		slider_callback_ = vtkSmartPointer<vtkSliderCallback>::New();
		slider_callback_->step_ = &step_;
		slider_callback_->play_ = &play_;
		slider_callback_->step_play_ = &step_play_;
		slider_widget_->AddObserver(vtkCommand::InteractionEvent, slider_callback_);

		//////////////////////////////////////////////////////////////////////////////////////////////
		// play_ button widget_ ////////////////////////////////////////////////////////////////////////
		vtkSmartPointer<vtkImageData> imagePlay = vtkSmartPointer<vtkImageData>::New();
		vtkSmartPointer<vtkImageData> imagePause = vtkSmartPointer<vtkImageData>::New();
		CreateImagePlay(imagePlay);
		CreateImagePause(imagePause);
		play_button_ = CreateOneOf<CommandButton>();
		play_button_->SetTwoStatesButtonContent(imagePlay, imagePause, render_window_interactor_);
		//int width = 20, height = 20;
		//play_button_->SetButtonSizePosition(width, height, winsize[0] - 235 - width, winsize[1] - 5 - height);
		vtkSmartPointer<PlayCallback> playButtonCallback = vtkSmartPointer<PlayCallback>::New();
		playButtonCallback->play_ = &play_;
		playButtonCallback->step_play_ = &step_play_;
		play_button_->SetButtonPressCallBack(playButtonCallback);
		//////////////////////////////////////////////////////////////////////////////////////////////
		// Single Step button widget_ /////////////////////////////////////////////////////////////////
		vtkSmartPointer<vtkImageData> imageNext = vtkSmartPointer<vtkImageData>::New();
		vtkSmartPointer<vtkImageData> imagePrev = vtkSmartPointer<vtkImageData>::New();
		CreateImageNextStep(imageNext);
		CreateImagePrevStep(imagePrev);
		next_button_ = CreateOneOf<CommandButton>();
		prev_button_ = CreateOneOf<CommandButton>();
		next_button_->SetOneStatesButtonContent(imageNext, render_window_interactor_);
		prev_button_->SetOneStatesButtonContent(imagePrev, render_window_interactor_);
	
		//next_button_->SetButtonSizePosition(width, height, winsize[0] - 210 - width, winsize[1] - 5 - height);
		//prev_button_->SetButtonSizePosition(width, height, winsize[0] - 185 - width, winsize[1] - 5 - height);
		// call back
		vtkSmartPointer<NextPressedCallback> nextPressedCallback = vtkSmartPointer<NextPressedCallback>::New();
		nextPressedCallback->step_ = &step_;
		nextPressedCallback->num_step_ = &num_step_;
		nextPressedCallback->step_play_ = &step_play_;

		next_button_->SetButtonPressCallBack(nextPressedCallback);

		vtkSmartPointer<PrevPressedCallback> prevPressedCallback = vtkSmartPointer<PrevPressedCallback>::New();
		prevPressedCallback->step_ = &step_;
		prevPressedCallback->step_play_ = &step_play_;

		prev_button_->SetButtonPressCallBack(prevPressedCallback);
	}
};

}

#endif
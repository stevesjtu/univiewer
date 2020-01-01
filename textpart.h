#pragma once
#ifndef TEXTPART_H
#define TEXTPART_H

#include "paramdefine.h"
// text 
#include "vtkTextActor.h"
#include "vtkTextProperty.h"
#include "vtkRenderWindow.h"
#include "vtkRenderWindowInteractor.h"

#include "vtkTextWidget.h"
#include "vtkTextRepresentation.h"
#include "vtkCoordinate.h"
#include "vtkCommand.h"

#include "model.h"
#include "auxfunc.h"


namespace univiewer {

// base press
class PressTextCallback : public vtkCommand
{
public:
	VTKSubClass(PressTextCallback)
	PressTextCallback() {}

	virtual void Execute(vtkObject *caller, unsigned long, void*)
	{
		vtkTextWidget *text_widget_ = reinterpret_cast<vtkTextWidget*>(caller);
		text_widget_->GetTextActor()->GetTextProperty()->SetOpacity(0.2);
	}

};
// base release
class ReleaseTextCallback : public vtkCommand
{
protected:
	void releaseCommon(vtkObject *caller)
	{
		vtkTextWidget *text_widget_ = reinterpret_cast<vtkTextWidget*>(caller);
		text_widget_->GetTextActor()->GetTextProperty()->SetOpacity(0.5);
	}
public:
	VTKSubClass(ReleaseTextCallback)
	ReleaseTextCallback() {}

	virtual void Execute(vtkObject *caller, unsigned long, void*)
	{
		releaseCommon(caller);
		// something to do
	}
	std::vector<vtkSmartPointer<vtkTextWidget> > text_widget_leaf_;
};

////////////////////////////////////////////////////////////////////////////////////////
// derived class ///////////////////////////////////////////////////////////////////////
class SystemReleaseTextCallback : public ReleaseTextCallback
{
public:
	VTKSubClass(SystemReleaseTextCallback)
	SystemReleaseTextCallback() {}

	virtual void Execute(vtkObject *caller, unsigned long, void*)
	{
		releaseCommon(caller);
		// something to do
		for (auto &tw : text_widget_leaf_) {
			tw->SetEnabled(!tw->GetEnabled());
		}
		render_window_interactor_->Render();

#ifdef WIN32
    std::string fpathname, fname;
    OpenFileDlg(fpathname, fname);
    std::cout << fpathname << " " << fname << std::endl;
#endif
	}

	vtkSmartPointer<vtkRenderWindowInteractor> render_window_interactor_;
};


class CommandText
{
private:
	vtkSmartPointer<vtkTextWidget> text_widget_;
	vtkSmartPointer<vtkTextActor> text_actor_;
	vtkSmartPointer<vtkTextRepresentation> text_representation_;
	vtkSmartPointer<ReleaseTextCallback> release_callback_;

	vtkSmartPointer<vtkRenderWindowInteractor> render_window_interactor_;

	std::string name_;
public:
	CommandText() {}
	virtual ~CommandText(){}

	vtkSmartPointer<vtkTextWidget> GetTextWidget() { return text_widget_; }
	vtkSmartPointer<vtkTextRepresentation> GetTextRepresentation() { return text_representation_; }
	
	void SetName(const std::string &str){name_ = str;}
	std::string GetName() { return name_; }

	void AddLeafNode(vtkSmartPointer<vtkTextWidget> leaf)
	{
		release_callback_->text_widget_leaf_.push_back(leaf);
		leaf->Off();
	}

	void SetRenderWindowInteractor(vtkSmartPointer<vtkRenderWindowInteractor> iren)
	{
		render_window_interactor_ = iren;
	}

	void SetCommandTextContent(const std::string &str, const double r, const double g, const double b, const double a, const int isbold)
	{
		// Create the widget_
		text_actor_ = vtkSmartPointer<vtkTextActor>::New();
		text_actor_->SetInput(str.c_str());
		text_actor_->GetTextProperty()->SetColor(r, g, b);
		text_actor_->GetTextProperty()->SetOpacity(a);
		text_actor_->GetTextProperty()->SetBold(isbold);
		//text_actor_->GetTextProperty()->SetFontFamilyToCourier();
    //text_actor_->GetTextProperty()->SetFontFamilyToArial();
    text_actor_->GetTextProperty()->SetFontFamilyToTimes();

    //text_actor_->GetTextProperty()->SetFontFamilyAsString("consolas");

    //text_actor_->GetTextProperty()->SetFontFamilyAsString("����");


		text_actor_->GetTextProperty()->SetJustificationToCentered();
		text_actor_->GetTextProperty()->SetVerticalJustificationToCentered();

		text_representation_ = vtkSmartPointer<vtkTextRepresentation>::New();
		text_representation_->SetShowBorderToOff();
		text_representation_->GetPositionCoordinate()->SetCoordinateSystemToDisplay();
		text_representation_->GetPosition2Coordinate()->SetCoordinateSystemToDisplay();	

		text_widget_ = vtkSmartPointer<vtkTextWidget>::New();
		text_widget_->SetRepresentation(text_representation_);
		text_widget_->SetInteractor(render_window_interactor_);
		text_widget_->SetTextActor(text_actor_);
		text_widget_->SelectableOn();
		text_widget_->ResizableOff();
		text_widget_->On();
	}

	void SetTextSizePosition(int leftbuttomX, int leftbuttomY, int height)
	{
		text_representation_->GetPositionCoordinate()->SetValue(leftbuttomX, leftbuttomY); // absolute value
		std::string str = text_representation_->GetText();
		text_representation_->GetPosition2Coordinate()->SetValue(str.length() * height / 2 , height); //relative value
		
	}

	void SetTextCallback(vtkSmartPointer<ReleaseTextCallback> releasetext)
	{
		vtkSmartPointer<PressTextCallback> pressCallback = vtkSmartPointer<PressTextCallback>::New();
		text_widget_->AddObserver(vtkCommand::StartInteractionEvent, pressCallback);
		release_callback_ = releasetext;
		text_widget_->AddObserver(vtkCommand::EndInteractionEvent, release_callback_);
	}

};




class CurrentTimer
{
private:
	vtkSmartPointer<vtkTextActor> text_actor_;

public:
	CurrentTimer() {};
	virtual ~CurrentTimer() {
  };

	vtkSmartPointer<vtkTextActor> &GetTextActor() { return text_actor_; }
	void SetTextActor(vtkSmartPointer<vtkRenderWindow> &render_window_) {
		text_actor_ = vtkSmartPointer<vtkTextActor>::New();
		//text_actor_->SetInput("Hello world");
		//text_actor_->SetPosition(80, 40);
		//text_actor_->GetTextProperty()->SetFontSize(24);
		//text_actor_->GetTextProperty()->SetColor(1.0, 0.0, 0.0);

		////////////////////
		text_actor_->GetTextProperty()->SetFontSize(16);
		text_actor_->GetTextProperty()->SetColor(1.0, 1.0, 1.0);
		text_actor_->SetInput("Current Time = 0");

    text_actor_->GetTextProperty()->SetFontFamilyToTimes();

		int *winsize = render_window_->GetSize();
		text_actor_->SetPosition(winsize[0] - 180, winsize[1] - 25);
		
	}
};

}


#endif
#pragma once
#ifndef TEXTPART_H
#define TEXTPART_H

#include "ParamDefine.h"
// text 
#include "vtkTextActor.h"
#include "vtkTextProperty.h"
#include "vtkRenderWindow.h"
#include "vtkRenderWindowInteractor.h"

#include "vtkTextWidget.h"
#include "vtkTextRepresentation.h"
#include "vtkCoordinate.h"
#include "vtkCommand.h"

using namespace std;

class pressTextCallback : public vtkCommand
{
public:
	static pressTextCallback *New()
	{
		return new pressTextCallback;
	}
	pressTextCallback() {}

	virtual void Execute(vtkObject *caller, unsigned long, void*)
	{
		vtkTextWidget *textWidget = reinterpret_cast<vtkTextWidget*>(caller);
		textWidget->GetTextActor()->GetTextProperty()->SetOpacity(0.2);
	}

};

class releaseTextCallback : public vtkCommand
{
public:
	static releaseTextCallback *New()
	{
		return new releaseTextCallback;
	}
	releaseTextCallback() {}

	virtual void Execute(vtkObject *caller, unsigned long, void*)
	{
		vtkTextWidget *textWidget = reinterpret_cast<vtkTextWidget*>(caller);
		textWidget->GetTextActor()->GetTextProperty()->SetOpacity(0.5);
		// something to do
	}
};

class CommandText
{
private:
	vtkSmartPointer<vtkTextWidget> textWidget;
	vtkSmartPointer<vtkTextActor> textActor;
	vtkSmartPointer<vtkTextRepresentation> textRepresentation;

	vtkSmartPointer<vtkRenderWindowInteractor> renderWindowInteractor;

public:
	CommandText() {}
	virtual ~CommandText(){}
	static shared_ptr<CommandText> New() { return make_shared<CommandText>(); }

	void setRenderWindowInteractor(vtkSmartPointer<vtkRenderWindowInteractor> iren)
	{
		renderWindowInteractor = iren;
	}

	void setCommandTextContent(const string &str)
	{
		// Create the widget
		textActor = vtkSmartPointer<vtkTextActor>::New();
		textActor->SetInput(str.c_str());
		textActor->GetTextProperty()->SetColor(1.0, 1.0, 1.0);
		textActor->GetTextProperty()->SetOpacity(0.5);
		textActor->GetTextProperty()->SetBold(1);

		textActor->GetTextProperty()->SetJustificationToCentered();
		textActor->GetTextProperty()->SetVerticalJustificationToCentered();

		textRepresentation = vtkSmartPointer<vtkTextRepresentation>::New();
		//textRepresentation->SetShowBorderToActive();
		textRepresentation->SetShowBorderToOn();
		
		textWidget = vtkSmartPointer<vtkTextWidget>::New();
		textWidget->SetRepresentation(textRepresentation);

		textWidget->SetInteractor(renderWindowInteractor);
		textWidget->SetTextActor(textActor);
		textWidget->SelectableOn();
		textWidget->ResizableOff();
		
		textWidget->On();
	}

	void setTextSizePosition()
	{
		int *winSize = renderWindowInteractor->GetRenderWindow()->GetSize();
		textRepresentation->GetPositionCoordinate()->SetCoordinateSystemToDisplay();
		textRepresentation->GetPosition2Coordinate()->SetCoordinateSystemToDisplay();
		
		textRepresentation->GetPositionCoordinate()->SetValue(10, winSize[1]-50); // absolute value
		textRepresentation->GetPosition2Coordinate()->SetValue(210, 40); //relative value
	}

	void setTextCallback()
	{
		vtkSmartPointer<pressTextCallback> pressCallback = vtkSmartPointer<pressTextCallback>::New();
		textWidget->AddObserver(vtkCommand::StartInteractionEvent, pressCallback);

		vtkSmartPointer<releaseTextCallback> releaseCallback = vtkSmartPointer<releaseTextCallback>::New();
		textWidget->AddObserver(vtkCommand::EndInteractionEvent, releaseCallback);
	}

};




class CurrentTimer
{
private:
	vtkSmartPointer<vtkTextActor> textActor;

public:
	CurrentTimer() {};
	virtual ~CurrentTimer() {};
	static shared_ptr<CurrentTimer> New()
	{
		return make_shared<CurrentTimer>();
	}

	vtkSmartPointer<vtkTextActor> &getTextActor() { return textActor; }
	void setTextActor(vtkSmartPointer<vtkRenderWindow> &renderWindow) {
		textActor = vtkSmartPointer<vtkTextActor>::New();
		//textActor->SetInput("Hello world");
		//textActor->SetPosition(80, 40);
		//textActor->GetTextProperty()->SetFontSize(24);
		//textActor->GetTextProperty()->SetColor(1.0, 0.0, 0.0);

		////////////////////
		textActor->GetTextProperty()->SetFontSize(16);
		textActor->GetTextProperty()->SetColor(1.0, 1.0, 1.0);
		textActor->SetInput("Current Time = 0");

		int *winsize = renderWindow->GetSize();
		textActor->SetPosition(winsize[0] - 180, winsize[1] - 25);
		
	}
};




#endif
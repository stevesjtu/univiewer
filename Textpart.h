#pragma once
#ifndef TEXTPART_H
#define TEXTPART_H

#include "ParamDefine.h"
// text 
#include "vtkTextActor.h"
#include "vtkTextProperty.h"

using namespace std;

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
	void setTextActor() {
		textActor = vtkSmartPointer<vtkTextActor>::New();
		//textActor->SetInput("Hello world");
		//textActor->SetPosition(80, 40);
		//textActor->GetTextProperty()->SetFontSize(24);
		//textActor->GetTextProperty()->SetColor(1.0, 0.0, 0.0);

		////////////////////
		textActor->GetTextProperty()->SetFontSize(16);
		textActor->GetTextProperty()->SetColor(1.0, 1.0, 1.0);

		textActor->SetPosition(800 - 170, 640 - 25);

	}
};




#endif
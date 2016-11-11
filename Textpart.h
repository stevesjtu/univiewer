#pragma once
#ifndef TEXTPART_H
#define TEXTPART_H
#include "vtkSmartPointer.h"
#include "vtkLabeledDataMapper.h"
#include "vtkXMLUnstructuredGridReader.h"

#include "vtkTextActor.h"
#include "vtkTextProperty.h"

#include<vector>
#include<memory>
using namespace std;

class Labelnode
{
private:
	vtkSmartPointer<vtkActor2D> labelActor;
	vtkSmartPointer<vtkLabeledDataMapper> labelMapper;
public:
	Labelnode() {};
	virtual ~Labelnode() {};

	static shared_ptr<Labelnode> New()
	{
		return make_shared<Labelnode>();
	}

	vtkSmartPointer<vtkActor2D> &getlabelActor() { return labelActor; }
	void setLabelActor(vector<vtkSmartPointer<vtkXMLUnstructuredGridReader>>&ugridReaders) {
		labelMapper = vtkSmartPointer<vtkLabeledDataMapper>::New();

		for (unsigned i = 0; i< (unsigned)ugridReaders.size(); ++i) {
			labelMapper->AddInputConnection(ugridReaders[i]->GetOutputPort());
		}
		labelActor = vtkSmartPointer<vtkActor2D>::New();
		labelActor->SetMapper(labelMapper);
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
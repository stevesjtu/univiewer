#pragma once
#ifndef CONTROLVIEWADDTEXT_H
#define CONTROLVIEWADDTEXT_H

#include"ControlView.h"
#include "vtkVectorText.h"

class ControlView_addText : public ControlView
{
protected:

	vtkSmartPointer<vtkTextActor> addTextActor;
    vector<string> *modelfiles;
    vector<string> *dispfiles;
    
public:
	virtual ~ControlView_addText() {};
	ControlView_addText() {};
	static shared_ptr<ControlView_addText> New() {
		shared_ptr<ControlView_addText> nw = make_shared<ControlView_addText>();
		return nw;
	}

	vtkSmartPointer<vtkTextActor> & getAddTextActor() { return addTextActor; }
    void setfileName(vector<string> *model, vector<string> *disp){
        modelfiles = model;
        dispfiles = disp;
    }
    
	void AddText() {
		addTextActor = vtkSmartPointer<vtkTextActor>::New();
        
        stringstream ss("");
        ss << "Model from: ";
        for(unsigned i=0;i< modelfiles->size()-1; ++i)
            ss << modelfiles->at(i)<<", ";
        ss << modelfiles->at(modelfiles->size()-1)<<"\n";
        
        ss << "Result from: ";
        if(!dispfiles->empty())
            ss<< dispfiles->at(0);
        
		addTextActor->SetInput(ss.str().c_str());
		addTextActor->SetPosition(240, 40);
		addTextActor->GetTextProperty()->SetFontSize(16);
		addTextActor->GetTextProperty()->SetColor(0.0, 0.0, 0.0);
		addTextActor->GetTextProperty()->SetOpacity(0.2);
		renderer->AddActor2D(addTextActor);
	}
    

    
};





#endif //CONTROLVIEWADDTEXT_H
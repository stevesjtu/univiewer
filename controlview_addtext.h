#pragma once
#ifndef CONTROLVIEWADDTEXT_H
#define CONTROLVIEWADDTEXT_H

#include"controlview.h"

namespace univiewer {

class ControlView_addText : public ControlView
{
protected:

	vtkSmartPointer<vtkTextActor> add_text_actor_;
  std::vector<std::string> *modelfiles;
  std::vector<std::string> *dispfiles;
    
public:
	virtual ~ControlView_addText() {};
	ControlView_addText() {};

	vtkSmartPointer<vtkTextActor> & getAddTextActor() { return add_text_actor_; }
    void setfileName(std::vector<std::string> *model, std::vector<std::string> *disp){
        modelfiles = model;
        dispfiles = disp;
    }
    
	void AddText() {
		add_text_actor_ = vtkSmartPointer<vtkTextActor>::New();
        
        std::stringstream ss("");
        ss << "Model from: ";
        for(unsigned i=0;i< modelfiles->size()-1; ++i)
            ss << modelfiles->at(i)<<", ";
        ss << modelfiles->at(modelfiles->size()-1)<<"\n";
        
        ss << "Result from: ";
        if(!dispfiles->empty())
            ss<< dispfiles->at(0);
        
		add_text_actor_->SetInput(ss.str().c_str());
		add_text_actor_->SetPosition(240, 40);
		add_text_actor_->GetTextProperty()->SetFontSize(16);
		add_text_actor_->GetTextProperty()->SetColor(0.0, 0.0, 0.0);
		add_text_actor_->GetTextProperty()->SetOpacity(0.2);
		renderer_->AddActor2D(add_text_actor_);
	}
    

    
};



}

#endif //CONTROLVIEWADDTEXT_H
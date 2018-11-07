#pragma once
#ifndef CONTROLVIEW_H
#define CONTROLVIEW_H

#include "model.h"
#include "axespart.h"
#include "textpart.h"
#include "sliderbarpart.h"
#include "lutpart.h"
#include "plotpart.h"
#include "auxfunc.h"

#include "vtkProperty.h"
// for animation
#include "vtkProgrammableFilter.h"
#include "vtkCallbackCommand.h"
// for windows render and view
#include "vtkRenderWindow.h"
#include "vtkRenderer.h"
#include "vtkRenderWindowInteractor.h"
#include "vtkCamera.h"
#include "vtkInteractorStyleTrackballCamera.h"

#define DATA_RESET 0x0000
#define DATA_MODEL 0x0001
#define DATA_DISPL 0x0002
#define DATA_NODVL 0x0004
#define DATA_CONPR 0x0008

namespace univiewer {

// Timer callback subclass
class CommandSubclass : public vtkCommand
{
public:
	vtkTypeMacro(CommandSubclass, vtkCommand);
	static vtkSmartPointer<CommandSubclass> New()
	{
		return new CommandSubclass;
	}
	virtual void Execute(vtkObject *caller, unsigned long vtkNotUsed(eventId),
		void *vtkNotUsed(callData)) override
	{
		vtkRenderWindowInteractor *iren =
			static_cast<vtkRenderWindowInteractor*>(caller);
		this->ProgrammableFilter->Modified();
		iren->Render();
	}
	vtkSmartPointer<vtkProgrammableFilter> ProgrammableFilter;
};

// real callback function
void TimerCallback(void* arguments);
void WindowModifiedCallback(vtkObject* caller, long unsigned int eventId, void* clientData, void* callData);
void KeypressCallback(vtkObject* caller, long unsigned int eventId, void* clientData, void* callData);

class ControlView
{
protected:
	// for the centering point
	vtkSmartPointer<vtkActor> actor;
	vtkSmartPointer<vtkDataSetMapper> mapper;

	// for windows control
	vtkSmartPointer<vtkRenderWindow> renderWindow;
	vtkSmartPointer<vtkRenderer> renderer;
	vtkSmartPointer<vtkRenderWindowInteractor> renderWindowInteractor;
	vtkSmartPointer<vtkInteractorStyleTrackballCamera> style;

	// for some part of graphics
	shared_ptr<Axesline> axesline;
	shared_ptr<Axesframe> axesframe;
	shared_ptr<Sliderbar> sliderbar;
	shared_ptr<CurrentTimer> currenttimer;
	shared_ptr<LookUpTable> lookuptable;

	shared_ptr<CommandText> cText;
	vector<shared_ptr<CommandText> > cTextBodies;
	// for animation and UI
	vtkSmartPointer<vtkProgrammableFilter> programmableFilter;

	// main part of visulization
	vector<shared_ptr<Model> > pModels;
	vector<shared_ptr<ContactData> > pContacts;
	//////////////////////////////////
	//////////////////////////////////
	int data_type;
	bool play, stepPlay;
	bool ShowMarker;
	bool ShowMesh;
  bool ShowLabel;

	//int nextButtonState, prevButtonState;
  void readSimpleOutModel(ifstream &infile, vector<unsigned int> &modelinfo,
    vector<unsigned int> &elemlist, vector<double> &nodelist);

public:
	virtual ~ControlView() {}
	ControlView(): play(false), stepPlay(false), ShowMarker(true), ShowMesh(true), ShowLabel(false) {};
	static shared_ptr<ControlView> New()
	{
		return make_shared<ControlView>();
	}
	
	inline bool & IsPlay() { return play; }
	inline bool & IsStepPlay() { return stepPlay; }
	inline bool & IsShowMarker() { return ShowMarker; }
	inline bool & IsShowMesh() { return ShowMesh; }
  inline bool & IsShowLabel() { return ShowLabel; }

	vector<shared_ptr<Model> > & getModels() { return pModels; }
	vector<shared_ptr<ContactData> > &getContactData() { return pContacts; }
	vtkSmartPointer<vtkProgrammableFilter> getProgrammableFilter() { return programmableFilter; }

	shared_ptr<CommandText> getCommandText() { return cText; }
	vector<shared_ptr<CommandText> > &getCommandTextBodies() { return cTextBodies; }
	shared_ptr<CurrentTimer> getCurrentTimer() {return currenttimer;}
	shared_ptr<Axesline> getAxesline() { return axesline; }
	shared_ptr<Sliderbar> getSliderbar() { return sliderbar; }
	shared_ptr<LookUpTable> getLookuptable() { return lookuptable; }

	vtkSmartPointer<vtkRenderer> getRenderer(){return renderer;}
	vtkSmartPointer<vtkRenderWindowInteractor> getRenderWindowInteractor(){return renderWindowInteractor;}

  void readDispfile(const vector<string> & filename);

  int readSimpleOutResult(const string& filename);

	int inputModelfiles(std::vector<std::string> &argv);

	virtual void setAnimationMethod( void(*f)(void*)) {
		if (data_type & DATA_DISPL) {
			programmableFilter->SetExecuteMethod(f, this);
			renderWindowInteractor->CreateRepeatingTimer(10);
			vtkSmartPointer<CommandSubclass> timerCallback = vtkSmartPointer<CommandSubclass>::New();
			timerCallback->ProgrammableFilter = programmableFilter;
			renderWindowInteractor->AddObserver(vtkCommand::TimerEvent, timerCallback);
		}
	}

	virtual void setKeyboardMethod(void (*f)(vtkObject* , long unsigned int vtkNotUsed(eventId), void* clientData, void* vtkNotUsed(callData))){
		vtkSmartPointer<vtkCallbackCommand> keypressCallback = vtkSmartPointer<vtkCallbackCommand>::New();
		keypressCallback->SetCallback(f);
		keypressCallback->SetClientData(this);
		renderWindowInteractor->AddObserver(vtkCommand::KeyPressEvent, keypressCallback);
	}

	virtual void setWindowMethod(void (*f)(vtkObject* caller, long unsigned int vtkNotUsed(eventId), void* clientData, void* vtkNotUsed(callData))){
		vtkSmartPointer<vtkCallbackCommand> windowCallback = vtkSmartPointer<vtkCallbackCommand>::New();
		windowCallback->SetCallback(f);
		windowCallback->SetClientData(this);
		renderWindow->AddObserver(vtkCommand::ModifiedEvent, windowCallback);
	}

	void playThisStep(unsigned &step)
	{
		stringstream ss("");
		ss << "Current Time = " << Model::stepCollection[step];
		currenttimer->getTextActor()->SetInput(ss.str().c_str());

		sliderbar->getSliderRep()->SetValue(step);

		for (auto &pmodel : pModels) {
			pmodel->updateDisp(step);
		}

		for (auto &pcontact : pContacts) {
			pcontact->UpdateUGrid(step);
		}
	}

	virtual void Update() {
			
		if (stepPlay) {
			play = false;
			sliderbar->getPlayButton()->setState(0);
			playThisStep(Model::step);
		}
		
		if (play) {
			Model::step = (Model::step == Model::stepNum - 1) ? 0 : Model::step + 1;
			playThisStep(Model::step);
		}
	}

	virtual void setRender() {
		// create an environment actor (center point)
		vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();
		const double p[3] = { 0.0, 0.0, 0.0 };
		// Create the topology of the point (a vertex)
		vtkSmartPointer<vtkCellArray> vertices = vtkSmartPointer<vtkCellArray>::New();
		vtkIdType pid[1];
		pid[0] = points->InsertNextPoint(p);
		vertices->InsertNextCell(1, pid);
		// Create a polydata object
		vtkSmartPointer<vtkPolyData> point = vtkSmartPointer<vtkPolyData>::New();
		// Set the points and vertices we created as the geometry and topology of the polydata
		point->SetPoints(points);
		point->SetVerts(vertices);
		
		programmableFilter = vtkSmartPointer<vtkProgrammableFilter>::New();
		programmableFilter->SetInputData(point);
		
		mapper = vtkSmartPointer<vtkDataSetMapper>::New();
		mapper->SetInputConnection(programmableFilter->GetOutputPort());
		actor = vtkSmartPointer<vtkActor>::New();
		actor->SetMapper(mapper);
		actor->GetProperty()->SetColor(0.0, 0.0, 0.0);
		actor->GetProperty()->SetPointSize(1.0);

		// Create a renderer, render window, and interactor
		renderer = vtkSmartPointer<vtkRenderer>::New();
		renderer->SetNearClippingPlaneTolerance(1e-4);
		renderer->GetActiveCamera()->SetClippingRange(1e-4, 1000);
		renderer->SetAmbient(1.0, 1.0, 1.0);
		renderer->SetLightFollowCamera(1);

		renderWindow = vtkSmartPointer<vtkRenderWindow>::New();
		renderWindow->AddRenderer(renderer);
		renderWindow->SetSize(800, 640);
		renderWindow->SetWindowName("MLV 2.0");
    renderWindow->Render();
		renderWindow->PointSmoothingOn();
		renderWindow->LineSmoothingOn();
		renderWindow->PolygonSmoothingOn();

		renderWindowInteractor = vtkSmartPointer<vtkRenderWindowInteractor>::New();
		renderWindowInteractor->SetRenderWindow(renderWindow);
		style = vtkSmartPointer<vtkInteractorStyleTrackballCamera>::New();
		renderWindowInteractor->SetInteractorStyle(style);

		// Initialize must be called prior to creating timer events.
		renderWindowInteractor->Initialize();

	}

	virtual void Display() {

		// command text
		cText = CommandText::New();
		cText->setRenderWindowInteractor(renderWindowInteractor);
		int *winSize = renderWindowInteractor->GetRenderWindow()->GetSize();
		cText->setCommandTextContent("System View", 1.0, 1.0, 1.0, 0.5, 1);
		//cText->setTextSizePosition(10, winSize[1] - 30, 20);
		vtkSmartPointer<systemReleaseTextCallback> systemRelease = vtkSmartPointer<systemReleaseTextCallback>::New();
		systemRelease->renderWindowInteractor = renderWindowInteractor;
		cText->setTextCallback(systemRelease);

		std::stringstream ss("");
		for (unsigned i = 0; i < pModels.size(); ++i) {
			shared_ptr<CommandText> ct = CommandText::New();
			ct->setRenderWindowInteractor(renderWindowInteractor);
			ss.str("");
			ss.clear();
			ss << "|- Body_" << i;
			ct->setCommandTextContent(ss.str(), 1.0, 1.0, 1.0, 0.5, 1);
			//ct->setTextSizePosition(18, winSize[1] - 30 - 15* (i+1), 15);
			vtkSmartPointer<releaseTextCallback> release = vtkSmartPointer<releaseTextCallback>::New();
			ct->setTextCallback(release);
			cTextBodies.push_back(ct);
			cText->addLeafNode(ct->getTextWidget());
		}


		// axesline part
		axesline = Axesline::New();
		axesline->setAxesActor();
		renderer->AddActor(axesline->getAxesActor());

		// axesframe part
		axesframe = Axesframe::New();
		axesframe->setAxesWidget(renderWindowInteractor);

		// labelnodes part
		for (auto& pmodel : pModels) {
			pmodel->setLabelnode();
			renderer->AddActor2D(pmodel->getLabelactor());
			pmodel->getLabelactor()->VisibilityOff();
		}
		
		if (data_type & DATA_NODVL) {
			// lookuptable part
			lookuptable = LookUpTable::New();
			lookuptable->setScalars(pModels);
			renderer->AddActor2D(lookuptable->getScalarBar());
		}

		// Add the actor to the scene
		renderer->AddActor(actor);
		
		renderer->GradientBackgroundOn();
		renderer->SetBackground2(13.0 / 255.0, 71.0 / 255.0, 161.0 / 255.0);
		renderer->SetBackground(144.0 / 255.0, 202.0 / 255.0, 249.0 / 255.0);
		// Render and interact

		renderWindowInteractor->Start();

		//plotpart p;
		//p.plot(renderWindow, renderer);
	}

};

}

#endif //CONTROLVIEW_H
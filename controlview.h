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
	VTKSubClass(CommandSubclass)

	virtual void Execute(vtkObject *caller, unsigned long vtkNotUsed(eventId),
		void *vtkNotUsed(callData)) override
	{
		vtkRenderWindowInteractor *iren =
			static_cast<vtkRenderWindowInteractor*>(caller);
		this->programmable_filter_->Modified();
		iren->Render();
	}
	vtkSmartPointer<vtkProgrammableFilter> programmable_filter_;
};

// real callback function
void TimerCallback(void* arguments);
void WindowModifiedCallback(vtkObject* caller, long unsigned int eventId, void* clientData, void* callData);
void KeypressCallback(vtkObject* caller, long unsigned int eventId, void* clientData, void* callData);

class ControlView
{
protected:
	// for the centering point
	vtkSmartPointer<vtkActor> actor_;
	vtkSmartPointer<vtkDataSetMapper> mapper_;

	// for windows control
	vtkSmartPointer<vtkRenderWindow> render_window_;
	vtkSmartPointer<vtkRenderer> renderer_;
	vtkSmartPointer<vtkRenderWindowInteractor> render_window_interactor_;
	vtkSmartPointer<vtkInteractorStyleTrackballCamera> interactor_style_;

	// for some part of graphics
	sptr<Axesline> axesline_;
	sptr<Axesframe> axesframe_;
	sptr<SliderBar> sliderbar_;
	sptr<CurrentTimer> currenttimer_;
	sptr<LookUpTable> lookuptable_;

	sptr<CommandText> ctext_;
	std::vector<sptr<CommandText> > ctext_bodies_;
	// for animation and UI
	vtkSmartPointer<vtkProgrammableFilter> programmable_filter_;

	// main part of visulization
	std::vector<sptr<Model> > models_;
	std::vector<sptr<ContactData> > contacts_;
	//////////////////////////////////
	//////////////////////////////////
	int data_type_;
	bool play_, step_play_;
	bool show_marker_;
	bool show_mesh_;
  bool show_label_;

	//int nextButtonState, prevButtonState;
  void ReadSimpleOutModel(std::ifstream &infile, std::vector<unsigned int> &modelinfo,
    std::vector<unsigned int> &elemlist, std::vector<double> &nodelist);

public:
	virtual ~ControlView() {}
	ControlView(): 
    play_(false), 
    step_play_(false), 
    show_marker_(true), 
    show_mesh_(true), 
    show_label_(false),
    data_type_(DATA_RESET) {};
	
  virtual void Reset();

	inline bool & IsPlay() { return play_; }
	inline bool & IsStepPlay() { return step_play_; }
	inline bool & IsShowMarker() { return show_marker_; }
	inline bool & IsShowMesh() { return show_mesh_; }
  inline bool & IsShowLabel() { return show_label_; }

	std::vector<sptr<Model> > & GetModels() { return models_; }
	std::vector<sptr<ContactData> > &GetContactData() { return contacts_; }
	vtkSmartPointer<vtkProgrammableFilter> GetProgrammableFilter() { return programmable_filter_; }

	sptr<CommandText> GetCommandText() { return ctext_; }
	std::vector<sptr<CommandText> > &GetCommandTextBodies() { return ctext_bodies_; }
	sptr<CurrentTimer> GetCurrentTimer() {return currenttimer_;}
	sptr<Axesline> GetAxesline() { return axesline_; }
	sptr<SliderBar> GetSliderbar() { return sliderbar_; }
	sptr<LookUpTable> GetLookuptable() { return lookuptable_; }

	vtkSmartPointer<vtkRenderer> GetRenderer(){return renderer_;}
	vtkSmartPointer<vtkRenderWindowInteractor> GetRenderWindowInteractor(){return render_window_interactor_;}

  void ReadDispFile(const std::vector<std::string> & filename);

  int ReadSimpleOutResult(const std::string& filename);

	int InputModelfiles(std::vector<std::string> &argv);

	virtual void SetAnimationMethod( void(*f)(void*)) {
		if (data_type_ & DATA_DISPL) {
			programmable_filter_->SetExecuteMethod(f, this);
			render_window_interactor_->CreateRepeatingTimer(10);
			vtkSmartPointer<CommandSubclass> timerCallback = vtkSmartPointer<CommandSubclass>::New();
			timerCallback->programmable_filter_ = programmable_filter_;
			render_window_interactor_->AddObserver(vtkCommand::TimerEvent, timerCallback);
		}
	}

	virtual void SetKeyboardMethod(void (*f)(vtkObject* , long unsigned int vtkNotUsed(eventId), void* clientData, void* vtkNotUsed(callData))){
		vtkSmartPointer<vtkCallbackCommand> keypressCallback = vtkSmartPointer<vtkCallbackCommand>::New();
		keypressCallback->SetCallback(f);
		keypressCallback->SetClientData(this);
		render_window_interactor_->AddObserver(vtkCommand::KeyPressEvent, keypressCallback);
	}

	virtual void SetWindowMethod(void (*f)(vtkObject* caller, long unsigned int vtkNotUsed(eventId), void* clientData, void* vtkNotUsed(callData))){
		vtkSmartPointer<vtkCallbackCommand> windowCallback = vtkSmartPointer<vtkCallbackCommand>::New();
		windowCallback->SetCallback(f);
		windowCallback->SetClientData(this);
		render_window_->AddObserver(vtkCommand::ModifiedEvent, windowCallback);
	}

	void PlayThisStep(unsigned &step_)
	{
		std::stringstream ss("");
		ss << "Current Time = " << Model::step_collection_[step_];
		currenttimer_->GetTextActor()->SetInput(ss.str().c_str());

		sliderbar_->GetSliderRep()->SetValue(step_);

		for (auto &pmodel : models_) {
			pmodel->UpdateDisp(step_);
		}

		for (auto &pcontact : contacts_) {
			pcontact->UpdateUGrid(step_);
		}
	}

	virtual void Update() {
			
		if (step_play_) {
			play_ = false;
			sliderbar_->GetPlayButton()->SetState(0);
			PlayThisStep(Model::step_);
		}
		
		if (play_) {
			Model::step_ = (Model::step_ == Model::num_step_ - 1) ? 0 : Model::step_ + 1;
			PlayThisStep(Model::step_);
		}
	}

	virtual void SetRender() {
		// create an environment actor_ (center point)
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
		
		programmable_filter_ = vtkSmartPointer<vtkProgrammableFilter>::New();
		programmable_filter_->SetInputData(point);
		
		mapper_ = vtkSmartPointer<vtkDataSetMapper>::New();
		mapper_->SetInputConnection(programmable_filter_->GetOutputPort());
		actor_ = vtkSmartPointer<vtkActor>::New();
		actor_->SetMapper(mapper_);
		actor_->GetProperty()->SetColor(0.0, 0.0, 0.0);
		actor_->GetProperty()->SetPointSize(1.0);

		// Create a renderer_, render window, and interactor
		renderer_ = vtkSmartPointer<vtkRenderer>::New();
		renderer_->SetNearClippingPlaneTolerance(1e-4);
		renderer_->GetActiveCamera()->SetClippingRange(1e-4, 1000);
		renderer_->SetAmbient(1.0, 1.0, 1.0);
		renderer_->SetLightFollowCamera(1);

		render_window_ = vtkSmartPointer<vtkRenderWindow>::New();
		render_window_->AddRenderer(renderer_);
		render_window_->SetSize(800, 640);
		render_window_->SetWindowName("MLV 2.0");
    render_window_->Render();
		render_window_->PointSmoothingOn();
		render_window_->LineSmoothingOn();
		render_window_->PolygonSmoothingOn();

		render_window_interactor_ = vtkSmartPointer<vtkRenderWindowInteractor>::New();
		render_window_interactor_->SetRenderWindow(render_window_);
		interactor_style_ = vtkSmartPointer<vtkInteractorStyleTrackballCamera>::New();
		render_window_interactor_->SetInteractorStyle(interactor_style_);

		// Initialize must be called prior to creating timer events.
		render_window_interactor_->Initialize();

	}

	virtual void Display() {

		// command text
		ctext_ = CreateOneOf<CommandText>();
		ctext_->SetRenderWindowInteractor(render_window_interactor_);
		int *winSize = render_window_interactor_->GetRenderWindow()->GetSize();
		ctext_->SetCommandTextContent("System View", 1.0, 1.0, 1.0, 0.5, 1);
		//ctext_->SetTextSizePosition(10, winSize[1] - 30, 20);
		vtkSmartPointer<SystemReleaseTextCallback> systemRelease = vtkSmartPointer<SystemReleaseTextCallback>::New();
		systemRelease->render_window_interactor_ = render_window_interactor_;
		ctext_->SetTextCallback(systemRelease);

		std::stringstream ss("");
		for (unsigned i = 0; i < models_.size(); ++i) {
			sptr<CommandText> ct = CreateOneOf<CommandText>();
			ct->SetRenderWindowInteractor(render_window_interactor_);
			ss.str("");
			ss.clear();
			ss << "|- Body_" << i;
			ct->SetCommandTextContent(ss.str(), 1.0, 1.0, 1.0, 0.5, 1);
			//ct->SetTextSizePosition(18, winSize[1] - 30 - 15* (i+1), 15);
			vtkSmartPointer<ReleaseTextCallback> release = vtkSmartPointer<ReleaseTextCallback>::New();
			ct->SetTextCallback(release);
			ctext_bodies_.push_back(ct);
			ctext_->AddLeafNode(ct->GetTextWidget());
		}


		// axesline_ part
		axesline_ = CreateOneOf<Axesline>();
		axesline_->SetAxesActor();
		renderer_->AddActor(axesline_->GetAxesActor());

		// axesframe_ part
		axesframe_ = CreateOneOf<Axesframe>();
		axesframe_->SetAxesWidget(render_window_interactor_);

		// labelnodes part
		for (auto& pmodel : models_) {
			pmodel->SetLabelnode();
			renderer_->AddActor2D(pmodel->GetLabelActor());
			pmodel->GetLabelActor()->VisibilityOff();
		}
		
		if (data_type_ & DATA_NODVL) {
			// lookuptable_ part
			lookuptable_ = CreateOneOf<LookUpTable>();
			lookuptable_->SetScalars(models_);
			renderer_->AddActor2D(lookuptable_->GetScalarBar());
		}

		// Add the actor_ to the scene
		renderer_->AddActor(actor_);
		
		renderer_->GradientBackgroundOn();
		renderer_->SetBackground2(13.0 / 255.0, 71.0 / 255.0, 161.0 / 255.0);
		renderer_->SetBackground(144.0 / 255.0, 202.0 / 255.0, 249.0 / 255.0);
		// Render and interact

		render_window_interactor_->Start();

		//PlotPart p;
		//p.plot(render_window_, renderer_);
	}

};

}

#endif //CONTROLVIEW_H
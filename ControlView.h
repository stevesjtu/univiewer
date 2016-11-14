#pragma once
#ifndef CONTROLVIEW_H
#define CONTROLVIEW_H
#include "vtkXMLUnstructuredGridReader.h"
#include "vtkUnstructuredGrid.h"
#include "vtkAppendFilter.h"

#include "vtkSmartPointer.h"
#include "vtkSphereSource.h"

#include "vtkDataSetMapper.h"
#include "vtkActor.h"
#include "vtkProperty.h"

#include "vtkRenderWindow.h"
#include "vtkRenderer.h"

#include "vtkRenderWindowInteractor.h"
#include "vtkCamera.h"
#include "vtkInteractorStyleTrackballCamera.h"

#include "vtkProgrammableFilter.h"
#include "vtkCallbackCommand.h"

#include "Model.h"
#include "Axespart.h"
#include "Textpart.h"
#include "Sliderbarpart.h"
#include "lutpart.h"

#define AXESLINE_PART		0x0001
#define AXESFRAME_PART		0x0002
#define SLIDEBAR_PART		0x0004
#define LABLENODE_PART		0x0008
#define CURRENTTIMER_PART	0x0010
#define LOOKUPTABLE_PART	0x0020

#define DEFAULT_TIMERCALLBACK TimerCallback
#define DEFAULT_KEYPRESSCALLBACK KeypressCallback
#define DEFAULT_WINDOWCALLBACK WindowModifiedCallback

// Timer callback subclass
class CommandSubclass : public vtkCommand
{
public:
	vtkTypeMacro(CommandSubclass, vtkCommand);
	static CommandSubclass *New()
	{
		return new CommandSubclass;
	}
	void Execute(vtkObject *caller, unsigned long vtkNotUsed(eventId),
		void *vtkNotUsed(callData))
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
	vtkSmartPointer<vtkActor> actor;

    vector<vtkSmartPointer<vtkXMLUnstructuredGridReader> > ugridReaders;
	vtkSmartPointer<vtkDataSetMapper> mapper;

	vtkSmartPointer<vtkRenderWindow> renderWindow;
	vtkSmartPointer<vtkRenderer> renderer;
	vtkSmartPointer<vtkRenderWindowInteractor> renderWindowInteractor;
	vtkSmartPointer<vtkInteractorStyleTrackballCamera> style;

	//
	shared_ptr<Axesline> axesline;

	shared_ptr<Axesframe> axesframe;

	shared_ptr<Labelnode> labelnode;

	shared_ptr<Sliderbar> sliderbar;

	shared_ptr<CurrentTimer> currenttimer;

	shared_ptr<LookUpTable> lookuptable;
	//
	vtkSmartPointer<vtkProgrammableFilter> programmableFilter;
	vtkSmartPointer<CommandSubclass> timerCallback;

	vtkSmartPointer<vtkCallbackCommand> keypressCallback;
	vtkSmartPointer<vtkCallbackCommand> windowCallback;

	//
	shared_ptr<Model> pModel;
	//////////////////////////////////
	//////////////////////////////////
	unsigned stepNum;
	unsigned step;
	bool play;
	bool ShowMarker;
	bool ShowMesh;
    bool ShowLabel;
public:
	virtual ~ControlView() {}
	ControlView(): stepNum(0), step(0),  play(false), ShowMarker(true), ShowMesh(true), ShowLabel(false) {};
	static shared_ptr<ControlView> New()
	{
		shared_ptr<ControlView> nw = make_shared<ControlView>();
		return nw;
	}

	unsigned &getStepNum() { return stepNum; }
	unsigned & getStep() { return step; }
	bool & IsPlay() { return play; }
	bool & IsShowMarker() { return ShowMarker; }
	bool & IsShowMesh() { return ShowMesh; }
    bool & IsShowLabel() { return ShowLabel; }
	
	shared_ptr<Model> & getModel() { return pModel; }
	vtkSmartPointer<vtkProgrammableFilter> & getProgrammableFilter() { return programmableFilter; }
	vtkSmartPointer<vtkActor> &getActor() { return actor; }

	shared_ptr<CurrentTimer> &getCurrentTimer() {return currenttimer;}

	shared_ptr<Axesline> &getAxesline() { return axesline; }

    shared_ptr<Labelnode> &getLabelnode() { return labelnode; }

	shared_ptr<Sliderbar> &getSliderbar() { return sliderbar; }

    vtkSmartPointer<vtkRenderer> &getRenderer(){return renderer;}

    vtkSmartPointer<vtkRenderWindowInteractor> &getRenderWindowInteractor(){return renderWindowInteractor;}
    vtkSmartPointer<CommandSubclass> &getTimerCallback(){return timerCallback;}

	void inputModelfiles(vector<string> &modelFiles, vector<string>&dispFiles,
			const int& argc,  char* argv[]) 
	{
		if (argc == 1) {
			cout << "Type 'Univiewer /h' for more help." << endl;
			exit(0);
		}
		vector<string> *ptr_vector_name = NULL;
		for (int i = 1; i<argc; ++i) {

#ifdef __APPLE__
            if (*argv[i] == '-') {
#else
			if ((*argv[i] == '-') || (*argv[i] == '/')) {
#endif
				switch (*(argv[i] + 1)) {
				case 'm':
					ptr_vector_name = &modelFiles;
					break;
				case 'o':
					ptr_vector_name = &dispFiles;
					break;
				case 'h':
					cout << endl;
					cout << "Usage: Univiewer /m file1.xml file2.xml ... fileN.xml /o disp.dat" << endl;
					cout << endl;
					cout << "The 'fileN.xml' is a model file that is compatible with VTK API, and other programs. It is a text file using *.vtu format, but you can also add your own data." << endl;
					cout << "The 'disp.dat' is a data file that contains the displacements of nodes for all the models(file1.xml, file2.xml and so on). It is a binary file with all the data as the type of double, the detail data sequence is as follows:" << endl << endl;
					cout << "#####################################################################################################" << endl;
					cout << "time1 \n node_1_disp_x node_1_disp_y node_1_disp_z \n node_2_disp_x node_2_disp_y node_2_disp_z \n ... \n node_n_disp_x node_n_disp_y node_n_disp_z" << endl; 
					cout << "time2 \n ......" << endl;
					cout << "#####################################################################################################" << endl << endl;
					cout << "NOTE: Binary files are not same as txt files which contains symbols like \\n \\t, and all the data saved as String. A Binary file contains data saved as its own type without any symbols like \\n \\t." << endl;
					exit(0);
					break;
				default:
					break;
				}
				continue;
			}
			ptr_vector_name->push_back(argv[i]);
		}

		ugridReaders.resize(modelFiles.size());
		auto appendFilter = vtkSmartPointer<vtkAppendFilter>::New();
		for (unsigned i = 0; i< (unsigned)modelFiles.size(); ++i) {
			
			ugridReaders[i] = vtkSmartPointer<vtkXMLUnstructuredGridReader>::New();
			ugridReaders[i]->SetFileName(modelFiles[i].c_str());
			ugridReaders[i]->Update();
			appendFilter->AddInputConnection(ugridReaders[i]->GetOutputPort());
		}
		appendFilter->Update();
		programmableFilter = vtkSmartPointer<vtkProgrammableFilter>::New();
		programmableFilter->SetInputConnection(appendFilter->GetOutputPort());
	}

	void setContent(shared_ptr<Model> &model) {
		pModel = model;
		pModel->setVtkpnt0(programmableFilter->GetUnstructuredGridInput()->GetPoints());
	}

	void Update(){

        if(sliderbar->getButtonRepresentation()->GetState())
            play = true;
        else
            play = false;
		
		stringstream ss("");
		ss << "Current Time = " << pModel->getStep(step);
		currenttimer->getTextActor()->SetInput(ss.str().c_str());

		if (play) {
			step = (step == stepNum - 1) ? 0 : step + 1;
			//pModel->update(step);
			//sliderRep->SetTitleText(ss.str().c_str());
		}
		sliderbar->getSliderRep()->SetValue(step);
		
		unsigned numPts = programmableFilter->GetUnstructuredGridOutput()->GetNumberOfPoints();
		vtkSmartPointer<vtkDoubleArray> &scls = lookuptable->getScalars();
		for (unsigned i = 0; i < numPts; ++i) {
			scls->SetValue(i, static_cast<double>(sin((i + step)*0.1)));
		}

		programmableFilter->GetUnstructuredGridOutput()->GetPointData()->SetScalars(scls);
		programmableFilter->GetUnstructuredGridOutput()->SetPoints(pModel->getvtkPnts(step));	
		
	}

	void setAnimationMethod( void(*f)(void*), vector<string> &dispFiles ) {
		programmableFilter->SetExecuteMethod(f, this);
		pModel->readDispfile(dispFiles);
		stepNum = pModel->getStepNum();
		pModel->initialize();

		renderWindowInteractor->CreateRepeatingTimer(10);
		//timerCallback = vtkSmartPointer<vtkCallbackCommand>::New();
		//timerCallback->SetCallback(TimerCallbackFunction);
		//timerCallback->SetClientData(this);

		timerCallback = vtkSmartPointer<CommandSubclass>::New();
		timerCallback->ProgrammableFilter = programmableFilter;

		renderWindowInteractor->AddObserver(vtkCommand::TimerEvent, timerCallback);

	}

	void setKeyboardMethod(void (*f)(vtkObject* , long unsigned int vtkNotUsed(eventId), void* clientData, void* vtkNotUsed(callData))){
		keypressCallback = vtkSmartPointer<vtkCallbackCommand>::New();
		keypressCallback->SetCallback(f);
		keypressCallback->SetClientData(this);
		renderWindowInteractor->AddObserver(vtkCommand::KeyPressEvent, keypressCallback);
	}

	void setWindowMethod(void (*f)(vtkObject* caller, long unsigned int vtkNotUsed(eventId), void* clientData, void* vtkNotUsed(callData))){
		windowCallback = vtkSmartPointer<vtkCallbackCommand>::New();
		windowCallback->SetCallback(f);
		windowCallback->SetClientData(this);
		renderWindow->AddObserver(vtkCommand::ModifiedEvent, windowCallback);
	}

	void setMainActor() {
		// Create a mapper and actor

		mapper = vtkSmartPointer<vtkDataSetMapper>::New();
		
		mapper->SetInputConnection(programmableFilter->GetOutputPort());

		actor = vtkSmartPointer<vtkActor>::New();
		actor->SetMapper(mapper);
		//actor->GetProperty()->SetColor(0.9, 0.9, 0.9);
		actor->GetProperty()->SetEdgeColor(0.0, 0.0, 0.0);
		actor->GetProperty()->EdgeVisibilityOn();
		//actor->SetScale(0.1);
	}

	void setRender() {
		// Create a renderer, render window, and interactor
		renderer = vtkSmartPointer<vtkRenderer>::New();
		renderer->GetActiveCamera()->SetClippingRange(0.01, 1000);
		renderer->SetAmbient(1.0, 1.0, 1.0);
		renderer->SetLightFollowCamera(1);

		renderWindow = vtkSmartPointer<vtkRenderWindow>::New();
		renderWindow->AddRenderer(renderer);
		renderWindow->SetSize(800, 640);
		renderWindow->SetWindowName("MLV 2.0");
        renderWindow->Render();

		renderWindowInteractor = vtkSmartPointer<vtkRenderWindowInteractor>::New();
		renderWindowInteractor->SetRenderWindow(renderWindow);
		style = vtkSmartPointer<vtkInteractorStyleTrackballCamera>::New();
		renderWindowInteractor->SetInteractorStyle(style);

		// Initialize must be called prior to creating timer events.
		renderWindowInteractor->Initialize();

	}

	void Display(int TobeRegisteredParts) {

		if (TobeRegisteredParts & AXESLINE_PART) {
			axesline = Axesline::New();
			axesline->setAxesActor();
			renderer->AddActor(axesline->getAxesActor());
		}

		if (TobeRegisteredParts & AXESFRAME_PART) {
			axesframe = Axesframe::New();
			axesframe->setAxesWidget(renderWindowInteractor);
		}

		if (TobeRegisteredParts & LABLENODE_PART) {
			labelnode = Labelnode::New();
			labelnode->setLabelActor(ugridReaders);
			renderer->AddActor2D(labelnode->getlabelActor());
			labelnode->getlabelActor()->VisibilityOff();
		}

		if (TobeRegisteredParts & SLIDEBAR_PART) {
			sliderbar = Sliderbar::New();
			sliderbar->setSliderBar(renderer, renderWindowInteractor, step, stepNum);
		}

		if (TobeRegisteredParts & CURRENTTIMER_PART) {
			currenttimer = CurrentTimer::New();
			currenttimer->setTextActor();
			renderer->AddActor2D(currenttimer->getTextActor());
		}
		
		if (TobeRegisteredParts & LOOKUPTABLE_PART) {
			lookuptable = LookUpTable::New();
			lookuptable->setScalars(mapper, pModel->getNodenum());
			renderer->AddActor2D(lookuptable->getScalarBar());
		}
		
		// Add the actor to the scene
		renderer->AddActor(actor);
		
		renderer->GradientBackgroundOn();
		renderer->SetBackground2(13.0 / 255.0, 71.0 / 255.0, 161.0 / 255.0);
		renderer->SetBackground(144.0 / 255.0, 202.0 / 255.0, 249.0 / 255.0);
		// Render and interact

		renderWindowInteractor->Start();
	}

};


//void TimerCallbackFunction(vtkObject* caller, long unsigned int vtkNotUsed(eventId), void* clientData, void* vtkNotUsed(callData))
//{
//	
//	auto cv = static_cast<ControlView*>(clientData);
//	auto programmableFilter = cv->getProgrammableFilter();
//
//	auto *iren = static_cast<vtkRenderWindowInteractor*>(caller);
//	programmableFilter->Modified();
//	iren->Render();
//
//}

void TimerCallback(void* arguments)
{
	ControlView* pCtr = static_cast<ControlView*>(arguments);
	pCtr->Update();
}

void WindowModifiedCallback(vtkObject* caller, long unsigned int vtkNotUsed(eventId), void* clientData, void* vtkNotUsed(callData))
{
	
	vtkRenderWindow* window = static_cast<vtkRenderWindow*>(caller);
	int* windowSize = window->GetSize();
	ControlView* pCtr = static_cast<ControlView*>(clientData);
	
	static int isfirst = 0;   // a trick to prevent read Sliderbar before it is created.
	if (isfirst > 1) {
		if(pCtr->getStepNum()!=0){

			pCtr->getCurrentTimer()->getTextActor()->SetPosition(windowSize[0] - 170, windowSize[1] - 25);

			pCtr->getSliderbar()->getSliderRep()->GetPoint1Coordinate()->SetValue(windowSize[0] - 200 - windowSize[0]/4.0, windowSize[1] - 15);

			pCtr->getSliderbar()->getSliderRep()->GetPoint2Coordinate()->SetValue(windowSize[0] - 200, windowSize[1] - 15);

			vtkSmartPointer<vtkCoordinate> upperRight = vtkSmartPointer<vtkCoordinate>::New();
			upperRight->SetCoordinateSystemToNormalizedDisplay();
			upperRight->SetValue(1.0, 1.0);

			double bds[6];
			double sz = 20.0;

			bds[0] = upperRight->GetComputedDisplayValue(pCtr->getRenderer())[0] - sz - 175;
			bds[1] = bds[0] + sz;
			bds[2] = upperRight->GetComputedDisplayValue(pCtr->getRenderer())[1] - sz - 5;
			bds[3] = bds[2] + sz;
			bds[4] = bds[5] = 0.0;

			// Scale to 1, default is .5
			pCtr->getSliderbar()->getButtonRepresentation()->SetPlaceFactor(1.0);
			pCtr->getSliderbar()->getButtonRepresentation()->PlaceWidget(bds);
		}
	}
	++isfirst;
}

void KeypressCallback(vtkObject* caller, long unsigned int vtkNotUsed(eventId), void* clientData, void* vtkNotUsed(callData))
{

	vtkRenderWindowInteractor *iren =
		static_cast<vtkRenderWindowInteractor*>(caller);

	ControlView* pCtr = static_cast<ControlView*>(clientData);

    if(pCtr->getStepNum()!=0){
	    if (strcmp(iren->GetKeySym(), "space") == 0){
	        pCtr->IsPlay() = !pCtr->IsPlay();
			pCtr->getSliderbar()->getButtonRepresentation()->SetState(!pCtr->getSliderbar()->getButtonRepresentation()->GetState());
	//        if(pCtr->IsPlay()){
	//            pCtr->IsPlay() = !pCtr->IsPlay();
	//            pCtr->getRenderWindowInteractor()->InvokeEvent(vtkCommand::DisableEvent, pCtr->getTimerCallback());
	//        }
	//        else{
	//            pCtr->IsPlay() = !pCtr->IsPlay();
	//            pCtr->getRenderWindowInteractor()->InvokeEvent(vtkCommand::EnableEvent, pCtr->getTimerCallback());
	//        }

	    }


		if (strcmp(iren->GetKeySym(), "b") == 0) {
			pCtr->IsPlay() = false;
			pCtr->getSliderbar()->getButtonRepresentation()->SetState(0);
			pCtr->getStep() = (pCtr->getStep() == pCtr->getStepNum() - 1) ? pCtr->getStepNum() - 1 : pCtr->getStep() + 1;
		}

		if (strcmp(iren->GetKeySym(), "v") == 0) {
			pCtr->IsPlay() = false;
			pCtr->getSliderbar()->getButtonRepresentation()->SetState(0);
			pCtr->getStep() = (pCtr->getStep() == 0) ? 0 : pCtr->getStep() - 1;
		}
    }
	if (strcmp(iren->GetKeySym(), "i") == 0) {
		if (pCtr->IsShowMarker())
            pCtr->getAxesline()->getAxesActor()->VisibilityOff();
		else
			pCtr->getAxesline()->getAxesActor()->VisibilityOn();

		pCtr->IsShowMarker() = !pCtr->IsShowMarker();
	}

	if (strcmp(iren->GetKeySym(), "u") == 0) {
		if (pCtr->IsShowMesh())
			pCtr->getActor()->GetProperty()->EdgeVisibilityOff();
		else
			pCtr->getActor()->GetProperty()->EdgeVisibilityOn();

		pCtr->IsShowMesh() = !pCtr->IsShowMesh();
	}

	if (strcmp(iren->GetKeySym(), "l") == 0) {
        if (pCtr->IsShowLabel())
			pCtr->getLabelnode()->getlabelActor()->VisibilityOff();
        else
			pCtr->getLabelnode()->getlabelActor()->VisibilityOn();

        pCtr->IsShowLabel() = !pCtr->IsShowLabel();
	}

	iren->Render();
}
#endif //CONTROLVIEW_H
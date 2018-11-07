#include"controlview.h"

namespace univiewer {

////////////////////////////////////////////////////////////////////////////
// old reader for dispfile
////////////////////////////////////////////////////////////////////////////
void ControlView::readDispfile(const vector<string> & filename) {
  ifstream infile;
  infile.open(filename[0], ios::in | ios::binary);

  if (infile.is_open()) {
    vector<vector<double> > dispvecCollection;
    unsigned nodedeg = 3;
    if (filename.size() == 2) {
      nodedeg = atoi(filename[1].c_str());
    }
    unsigned dofs = Model::nodeNums * nodedeg;

    vector<double> datavec(dofs);
    double steptime;
    while (1) {
      infile.read((char*)&steptime, sizeof(double));
      infile.read((char*)datavec.data(), sizeof(double) * (dofs));

      if (infile.fail()) break;

      dispvecCollection.push_back(datavec);
      Model::stepCollection.push_back(steptime);
    }

    Model::stepNum = (unsigned)dispvecCollection.size();

    if (filename.size() == 2) {
      unsigned index;
      vector<double> dataPosition;
      for (unsigned s = 0; s < Model::stepNum; ++s) {
        auto &data = dispvecCollection[s];
        index = 0;
        dataPosition.clear();
        for (unsigned i = 0; i < Model::nodeNums; ++i) {
          dataPosition.push_back(data[index]);
          dataPosition.push_back(data[index + 1]);
          dataPosition.push_back(data[index + 2]);
          index += nodedeg;
        }
        dispvecCollection[s].swap(dataPosition);
      }
    }
    infile.close();

    // initialize the nodes position of Model feMesh
    vector<double> node0;
    for (auto& pmodel : pModels) {
      const auto &feMesh = pmodel->getFEMesh();

      node0.resize(pmodel->getNodenum() * 3);
      unsigned index = 0;
      for (unsigned n = 0; n < pmodel->getNodenum(); ++n) {
        double *xyz = feMesh->getUGrid()->GetPoint(n);
        node0[index++] = *xyz;
        node0[index++] = *(xyz + 1);
        node0[index++] = *(xyz + 2);
      }

      feMesh->getpvtkPnts().resize(Model::stepNum);

      double position[3];
      for (unsigned s = 0; s < Model::stepNum; ++s) {
        feMesh->getpvtkPnts(s) = vtkSmartPointer<vtkPoints>::New();
        feMesh->getpvtkPnts(s)->SetNumberOfPoints(pmodel->getNodenum());
        index = pmodel->getOffset();
        for (unsigned n = 0; n < pmodel->getNodenum(); ++n) {
          position[0] = node0[n * 3] + dispvecCollection[s][index];
          position[1] = node0[n * 3 + 1] + dispvecCollection[s][index + 1];
          position[2] = node0[n * 3 + 2] + dispvecCollection[s][index + 2];
          index += 3;
          feMesh->getpvtkPnts(s)->SetPoint(n, position);
        }
      }
    }
  } // if file opened
  else {
    cout << "Can not open dispfiles." << endl;
    exit(0);
  }

}

////////////////////////////////////////////////////////////////////////
// read model mesh
////////////////////////////////////////////////////////////////////////
void ControlView::readSimpleOutModel(ifstream &infile, vector<unsigned int> &modelinfo, vector<unsigned int> &elemlist, vector<double> &nodelist) {
  infile.read((char*)modelinfo.data(), sizeof(unsigned int) * 3);

  elemlist.resize(modelinfo[0] * modelinfo[2]);
  nodelist.resize(modelinfo[1] * 3);

  infile.read((char*)elemlist.data(), sizeof(unsigned int)* elemlist.size());
  infile.read((char*)nodelist.data(), sizeof(double)* nodelist.size());
}

//////////////////////////////////////////////////////////////////////////////
// read simple output results
//////////////////////////////////////////////////////////////////////////////
int ControlView::readSimpleOutResult(const string& filename) {
  ifstream infile(filename, ios::binary);
  if (!infile.is_open()) {
    cout << "Error in opening " << filename << endl;
    exit(1);
  }

  unsigned bodynum;
  infile.read((char*)&bodynum, sizeof(unsigned));

  vector<unsigned int> modelinfo(3);
  vector<unsigned int> elemlist;
  vector<double> nodelist;

  pModels.resize(bodynum);
  for (unsigned i = 0; i < bodynum; ++i) {
    this->readSimpleOutModel(infile, modelinfo, elemlist, nodelist);
    pModels[i] = Model::New();
    pModels[i]->CreateModel(modelinfo, elemlist, nodelist);
  }

  vector<unsigned> nodaldofs(bodynum);
  infile.read((char*)nodaldofs.data(), sizeof(unsigned)* nodaldofs.size());

  if (infile.fail()) {
    infile.close();
    return DATA_MODEL;
  } // means just model imported.

  //////////////////////////////////////////////////////////////
  // if goes here, read disp data 
  //////////////////////////////////////////////////////////////
  unsigned alldofs = 0;
  for(unsigned b=0; b< bodynum; ++b) { 
    alldofs += pModels[b]->getNodenum()* nodaldofs[b];
  }

  vector<vector<double> > dispdata;
  vector<double> disptemp(alldofs);
  
  double temp;
  while(true) { // for each step
    
    infile.read((char*)&temp, sizeof(double));
    if (infile.fail()) break;

    Model::stepCollection.push_back(temp);
    
    unsigned start = 0;
    for(unsigned b=0; b< bodynum; ++b) { // for each body
      pModels[b]->setOffset(start);
      infile.read((char*) (disptemp.data() + start), sizeof(double)* pModels[b]->getNodenum()* nodaldofs[b]);
      start += pModels[b]->getNodenum()* nodaldofs[b];
    }
    dispdata.push_back(disptemp);
  } 

  Model::stepNum = Model::stepCollection.size();
  infile.close();
  // finish reading disp data

  //////////////////////////////////////////////////////////////
  // initialize data base
  //////////////////////////////////////////////////////////////
  vector<double> node0;
  
  for (unsigned b=0; b< bodynum; ++b) {
    const auto &feMesh = pModels[b]->getFEMesh();

    node0.resize(pModels[b]->getNodenum() * 3);
    unsigned index = 0;
    for (unsigned n = 0; n < pModels[b]->getNodenum(); ++n) {
      double *xyz = feMesh->getUGrid()->GetPoint(n);
      node0[index++] = *xyz;
      node0[index++] = *(xyz + 1);
      node0[index++] = *(xyz + 2);
    }

    feMesh->getpvtkPnts().resize(Model::stepNum);
    double position[3];
    for (unsigned s = 0; s < Model::stepNum; ++s) {
      feMesh->getpvtkPnts(s) = vtkSmartPointer<vtkPoints>::New();
      feMesh->getpvtkPnts(s)->SetNumberOfPoints(pModels[b]->getNodenum());
      index = pModels[b]->getOffset();
      for (unsigned n = 0; n < pModels[b]->getNodenum(); ++n) {
        position[0] = node0[n * 3] + dispdata[s][index];
        position[1] = node0[n * 3 + 1] + dispdata[s][index + 1];
        position[2] = node0[n * 3 + 2] + dispdata[s][index + 2];
        index += nodaldofs[b];
        feMesh->getpvtkPnts(s)->SetPoint(n, position);
      }
    }
  }

  return (DATA_MODEL | DATA_DISPL); // means including displacement data for animation

}

int ControlView::inputModelfiles(std::vector<std::string> &argv) {
  vector<string> simple_out_result;
  vector<string> modelFiles;
  vector<string> dispFiles;
  vector<string> contFiles;
  vector<string> nodeDataFiles;
  argParser(argv, simple_out_result, modelFiles, dispFiles, contFiles, nodeDataFiles);
  
  data_type = DATA_RESET;
  if (!simple_out_result.empty()) {
    data_type = this->readSimpleOutResult(simple_out_result[0]);

    for(unsigned i=0; i < pModels.size(); ++i) {
      pModels[i]->getActor()->GetProperty()->SetEdgeColor(0.0, 0.0, 0.0);
      pModels[i]->getActor()->GetProperty()->EdgeVisibilityOn();

      if (pModels[i]->getFEMesh()->getUGrid()->GetCellType(0) == 4)
        pModels[i]->getActor()->GetProperty()->SetLineWidth(2.5);
      else
        pModels[i]->getActor()->GetProperty()->SetLineWidth(1.0);
      
      pModels[i]->getActor()->SetScale(1.0);
      renderer->AddActor(pModels[i]->getActor());
    }

    if (data_type & DATA_DISPL) {
      sliderbar = Sliderbar::New();
      sliderbar->setSliderBar(renderWindowInteractor, Model::step, Model::stepNum, this->play, this->stepPlay);

      currenttimer = CurrentTimer::New();
      currenttimer->setTextActor(renderWindow);
      renderer->AddActor2D(currenttimer->getTextActor());      
    }

  } else { // using simple output results or old fasion type input.
    pModels.resize(modelFiles.size());
    for (unsigned i = 0; i< (unsigned)modelFiles.size(); ++i) {
      pModels[i] = Model::New();
      pModels[i]->setOffset(Model::nodeNums* 3);

      if (modelFiles[i].substr(modelFiles[i].size() - 3).compare(".md") == 0) {
        pModels[i]->ReadTxtModel(modelFiles[i]);
      } else if (modelFiles[i].substr(modelFiles[i].size() - 3).compare("txt") == 0) {
        pModels[i]->ReadTxtModel(modelFiles[i]);
      } else if (modelFiles[i].substr(modelFiles[i].size() - 3).compare("xml") == 0) {
        pModels[i]->ReadXmlModel(modelFiles[i]);
      }
      
      //pModels[i]->getActor()->GetProperty()->SetColor(0.9, 0.9, 0.9);
      //pModels[i]->getActor()->GetProperty()->SetOpacity(1.0);
      pModels[i]->getActor()->GetProperty()->SetEdgeColor(0.0, 0.0, 0.0);
      pModels[i]->getActor()->GetProperty()->EdgeVisibilityOn();

      if (pModels[i]->getFEMesh()->getUGrid()->GetCellType(0) == 4)
        pModels[i]->getActor()->GetProperty()->SetLineWidth(2.5);
      else
        pModels[i]->getActor()->GetProperty()->SetLineWidth(1.0);
      
      pModels[i]->getActor()->SetScale(1.0);

      renderer->AddActor(pModels[i]->getActor());
    }

    if (!dispFiles.empty()) {
      data_type |= DATA_DISPL;

      this->readDispfile(dispFiles);

      sliderbar = Sliderbar::New();
      sliderbar->setSliderBar(renderWindowInteractor, Model::step, Model::stepNum, this->play, this->stepPlay);

      currenttimer = CurrentTimer::New();
      currenttimer->setTextActor(renderWindow);
      renderer->AddActor2D(currenttimer->getTextActor());
    }

    if (!contFiles.empty()) {
      data_type |= DATA_CONPR;

      readContfile(contFiles[0], pContacts, pModels);
      for (auto& pcontact : pContacts) {
        pcontact->InitializeUGrid();

        pcontact->getPActor()->GetProperty()->SetLineWidth(2.0);
        pcontact->getNActor()->GetProperty()->SetLineWidth(2.0);
        pcontact->getPActor()->GetProperty()->SetPointSize(6.0);
        pcontact->getNActor()->GetProperty()->SetPointSize(6.0);
        pcontact->getPActor()->GetProperty()->EdgeVisibilityOff();
        pcontact->getNActor()->GetProperty()->EdgeVisibilityOff();

        renderer->AddActor(pcontact->getPActor());
        renderer->AddActor(pcontact->getNActor());
      }
    }

    if (!nodeDataFiles.empty()) {
      data_type |= DATA_NODVL;
      readNodeDatafile(nodeDataFiles[0], pModels);
    }			
  }
  return 0;
}

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
	
	if(Model::stepNum!=0){
    
		pCtr->getCurrentTimer()->getTextActor()->SetPosition(windowSize[0] - 180, windowSize[1] - 25);

		pCtr->getSliderbar()->getSliderRep()->GetPoint1Coordinate()->SetValue(windowSize[0] - 260 - windowSize[0]/4.0, windowSize[1] - 15);
		pCtr->getSliderbar()->getSliderRep()->GetPoint2Coordinate()->SetValue(windowSize[0] - 260, windowSize[1] - 15);

		int width = 20, height = 20;
		pCtr->getSliderbar()->getPlayButton()->setButtonSizePosition(width, height, windowSize[0] - 235 - width, windowSize[1] - 5 - height);
		pCtr->getSliderbar()->getNextButton()->setButtonSizePosition(width, height, windowSize[0] - 210 - width, windowSize[1] - 5 - height);
		pCtr->getSliderbar()->getPrevButton()->setButtonSizePosition(width, height, windowSize[0] - 185 - width, windowSize[1] - 5 - height);
	}

	pCtr->getCommandText()->setTextSizePosition(10, windowSize[1] - 30, 20);
	for (unsigned i = 0; i < pCtr->getModels().size(); ++i) {
		pCtr->getCommandTextBodies()[i]->setTextSizePosition(10, windowSize[1] - 30 - 15 * (i+1), 15);
	}

	if (pCtr->getLookuptable()) {
		double w = windowSize[0] * 70.0 / 800.0;
		double h = windowSize[1] * 300.0 / 640.0;
		pCtr->getLookuptable()->getScalarBar()->SetPosition(windowSize[0] - w - 10, 10);
		pCtr->getLookuptable()->getScalarBar()->SetPosition2(w, h);
	}

}

void KeypressCallback(vtkObject* caller, long unsigned int vtkNotUsed(eventId), void* clientData, void* vtkNotUsed(callData))
{

	vtkRenderWindowInteractor *iren =
		static_cast<vtkRenderWindowInteractor*>(caller);

	ControlView* pCtr = static_cast<ControlView*>(clientData);

    if(Model::stepNum!=0){
	    if (strcmp(iren->GetKeySym(), "space") == 0){
	        pCtr->IsPlay() = !pCtr->IsPlay();
			if (pCtr->IsPlay()) {
				pCtr->getSliderbar()->getPlayButton()->setState(1);
				pCtr->IsStepPlay() = false;
			}
			else {
				pCtr->getSliderbar()->getPlayButton()->setState(0);
			}
				
	    }

		if (strcmp(iren->GetKeySym(), "b") == 0) {
			pCtr->IsPlay() = false;
			pCtr->getSliderbar()->getPlayButton()->setState(0);
			Model::step = (Model::step == Model::stepNum - 1) ? Model::stepNum - 1 : Model::step + 1;
			pCtr->IsStepPlay() = true;
		}

		if (strcmp(iren->GetKeySym(), "v") == 0) {
			pCtr->IsPlay() = false;
			pCtr->getSliderbar()->getPlayButton()->setState(0);
			Model::step = (Model::step == 0) ? 0 : Model::step - 1;
			pCtr->IsStepPlay() = true;
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
		if (pCtr->IsShowMesh()) {
			for(auto &pmodel: pCtr->getModels())
				pmodel->getActor()->GetProperty()->EdgeVisibilityOff();
		}
		else {
			for (auto &pmodel : pCtr->getModels())
				pmodel->getActor()->GetProperty()->EdgeVisibilityOn();
		}
			
		pCtr->IsShowMesh() = !pCtr->IsShowMesh();
	}

	if (strcmp(iren->GetKeySym(), "l") == 0) {

		if (pCtr->IsShowLabel()) {
			for(auto& pmodel : pCtr->getModels())
				pmodel->getLabelactor()->VisibilityOff();
		}
		else {
			for (auto& pmodel : pCtr->getModels())
				pmodel->getLabelactor()->VisibilityOn();
		}

		pCtr->IsShowLabel() = !pCtr->IsShowLabel();
	}

	iren->Render();
}  

}


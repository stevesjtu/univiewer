#include"controlview.h"
#include"omegafile.h"

namespace univiewer {

///////////////////////////////////////////////////////////////////////////
// Reset data base
///////////////////////////////////////////////////////////////////////////
void ControlView::Reset() {
  
  renderer_->RemoveActor2D(currenttimer_->GetTextActor());
  currenttimer_ = nullptr;

  //sliderbar_ = nullptr;
  //lookuptable_ = nullptr;

  //for (auto & model : models_) {
  //  renderer_->RemoveActor(model->GetActor());
  //  renderer_->RemoveActor2D(model->GetLabelActor());
  //}
  //
  //ctext_bodies_.clear();

  //if (!models_.empty()) models_.clear();
  //if (!contacts_.empty()) contacts_.clear();

  ///////////////////
  data_type_ = DATA_RESET;
  play_ = false;
  step_play_ = false;
  show_marker_ = true;
  show_mesh_ = true;
  show_label_ = false;

}

////////////////////////////////////////////////////////////////////////////
// old reader for dispfile
////////////////////////////////////////////////////////////////////////////
void ControlView::ReadDispFile(const std::vector<std::string> & filename) {
  std::ifstream infile;
  infile.open(filename[0], std::ios::in | std::ios::binary);

  if (infile.is_open()) {
    std::vector<std::vector<double> > dispvecCollection;
    Uint nodedeg = 3;
    if (filename.size() == 2) {
      nodedeg = atoi(filename[1].c_str());
    }
    Uint dofs = Model::num_nodes_ * nodedeg;

    std::vector<double> datavec(dofs);
    double steptime;
    while (1) {
      infile.read((char*)&steptime, sizeof(double));
      infile.read((char*)datavec.data(), sizeof(double) * (dofs));

      if (infile.fail()) break;

      dispvecCollection.push_back(datavec);
      Model::step_collection_.push_back(steptime);
    }

    Model::num_step_ = (Uint)dispvecCollection.size();

    if (filename.size() == 2) {
      Uint index;
      std::vector<double> dataPosition;
      for (Uint s = 0; s < Model::num_step_; ++s) {
        auto &data = dispvecCollection[s];
        index = 0;
        dataPosition.clear();
        for (Uint i = 0; i < Model::num_nodes_; ++i) {
          dataPosition.push_back(data[index]);
          dataPosition.push_back(data[index + 1]);
          dataPosition.push_back(data[index + 2]);
          index += nodedeg;
        }
        dispvecCollection[s].swap(dataPosition);
      }
    }
    infile.close();

    // initialize the nodes position of Model fe_mesh_
    std::vector<double> node0;
    for (auto& pmodel : models_) {
      const auto &fe_mesh_ = pmodel->GetFEMesh();

      node0.resize(pmodel->GetNumOfNode() * 3);
      Uint index = 0;
      for (Uint n = 0; n < pmodel->GetNumOfNode(); ++n) {
        double *xyz = fe_mesh_->GetUGrid()->GetPoint(n);
        node0[index++] = *xyz;
        node0[index++] = *(xyz + 1);
        node0[index++] = *(xyz + 2);
      }

      fe_mesh_->GetPvtkpnts().resize(Model::num_step_);

      double position[3];
      for (Uint s = 0; s < Model::num_step_; ++s) {
        fe_mesh_->GetPvtkpnts(s) = vtkSmartPointer<vtkPoints>::New();
        fe_mesh_->GetPvtkpnts(s)->SetNumberOfPoints(pmodel->GetNumOfNode());
        index = pmodel->GetOffset();
        for (Uint n = 0; n < pmodel->GetNumOfNode(); ++n) {
          position[0] = node0[n * 3] + dispvecCollection[s][index];
          position[1] = node0[n * 3 + 1] + dispvecCollection[s][index + 1];
          position[2] = node0[n * 3 + 2] + dispvecCollection[s][index + 2];
          index += 3;
          fe_mesh_->GetPvtkpnts(s)->SetPoint(n, position);
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
void ControlView::ReadSimpleOutModel(std::ifstream &infile, std::vector<Uint> &modelinfo, std::vector<Uint> &elemlist, std::vector<double> &nodelist) {
  infile.read((char*)modelinfo.data(), sizeof(Uint) * 3);

  elemlist.resize(modelinfo[0] * modelinfo[2]);
  nodelist.resize(modelinfo[1] * 3);

  infile.read((char*)elemlist.data(), sizeof(Uint)* elemlist.size());
  infile.read((char*)nodelist.data(), sizeof(double)* nodelist.size());
}

//////////////////////////////////////////////////////////////////////////////
// read simple output results
//////////////////////////////////////////////////////////////////////////////
int ControlView::ReadSimpleOutResult(const std::string& filename) {
  std::ifstream infile(filename, std::ios::binary);
  if (!infile.is_open()) {
    std::cout << "Error in opening " << filename << std::endl;
    exit(1);
  }

  Uint bodynum;
  infile.read((char*)&bodynum, sizeof(Uint));

  std::vector<Uint> modelinfo(3);
  std::vector<Uint> elemlist;
  std::vector<double> nodelist;

  models_.resize(bodynum);
  for (Uint i = 0; i < bodynum; ++i) {
    this->ReadSimpleOutModel(infile, modelinfo, elemlist, nodelist);
    models_[i] = CreateOneOf<Model>();
    models_[i]->CreateModel(modelinfo, elemlist, nodelist);
  }

  std::vector<Uint> bodydofs(bodynum);
  infile.read((char*)bodydofs.data(), sizeof(Uint)* bodydofs.size());

  if (infile.fail()) {
    infile.close();
    return DATA_MODEL;
  } // means just model imported.

  //////////////////////////////////////////////////////////////
  // if goes here, read disp data 
  //////////////////////////////////////////////////////////////
  Uint alldofs = 0;
  for(Uint b=0; b< bodynum; ++b) { 
    alldofs += bodydofs[b];
  }

  std::vector<std::vector<double> > dispdata;
  std::vector<double> disptemp(alldofs);
  
  double temp;
  while(true) { // for each step_
    
    infile.read((char*)&temp, sizeof(double));
    if (infile.fail()) break;

    Model::step_collection_.push_back(temp);
    
    Uint start = 0;
    for(Uint b=0; b< bodynum; ++b) { // for each body
      models_[b]->SetOffset(start);
      infile.read((char*) (disptemp.data() + start), sizeof(double)* bodydofs[b]);
      start += bodydofs[b];
    }
    dispdata.push_back(disptemp);
  } 

  Model::num_step_ = (Uint)Model::step_collection_.size();
  infile.close();
  // finish reading disp data

  //////////////////////////////////////////////////////////////
  // initialize data base
  //////////////////////////////////////////////////////////////
  std::vector<double> node0;
  
  for (Uint b=0; b< bodynum; ++b) {
    const auto &fe_mesh_ = models_[b]->GetFEMesh();

    node0.resize(models_[b]->GetNumOfNode() * 3);
    Uint index = 0;
    for (Uint n = 0; n < models_[b]->GetNumOfNode(); ++n) {
      double *xyz = fe_mesh_->GetUGrid()->GetPoint(n);
      node0[index++] = *xyz;
      node0[index++] = *(xyz + 1);
      node0[index++] = *(xyz + 2);
    }

    fe_mesh_->GetPvtkpnts().resize(Model::num_step_);
    double position[3];
    for (Uint s = 0; s < Model::num_step_; ++s) {
      fe_mesh_->GetPvtkpnts(s) = vtkSmartPointer<vtkPoints>::New();
      fe_mesh_->GetPvtkpnts(s)->SetNumberOfPoints(models_[b]->GetNumOfNode());
      index = models_[b]->GetOffset();
      for (Uint n = 0; n < models_[b]->GetNumOfNode(); ++n) {
        position[0] = node0[n * 3] + dispdata[s][index];
        position[1] = node0[n * 3 + 1] + dispdata[s][index + 1];
        position[2] = node0[n * 3 + 2] + dispdata[s][index + 2];
        index += 3;
        fe_mesh_->GetPvtkpnts(s)->SetPoint(n, position);

        //std::cout << "(" << position[0] << ", " << position[1] << ", " << position[2] << ")" << std::endl;

      }
    }
    models_[b]->UpdateDisp(0);
  }

  return (DATA_MODEL | DATA_DISPL); // means including displacement data for animation

}

int ControlView::ReadHDF5Result(const std::string &filename) {
  OmegaFile oft(filename);
  oft.Initialize();

  std::string meshtype;
  std::vector<Uint> meshinfo, elem;
  std::vector<double> node;
  models_.resize(oft.GetBodyName().size());

  for (Uint i = 0, isize = oft.GetBodyName().size(); i < isize; ++i) {
    meshtype.clear();
    meshinfo.clear();
    elem.clear();
    node.clear();
    oft.GetMeshByBodyName(oft.GetBodyName(i), meshtype, meshinfo, elem, node);
    models_[i] = CreateOneOf<Model>();
    models_[i]->CreateModel(meshinfo, elem, node);
  }

  oft.GetTimeSeries(Model::step_collection_);
  Model::num_step_ = (Uint)Model::step_collection_.size();

  if (Model::num_step_ == 0) return DATA_MODEL;

  // read nodal position at each step
  std::vector<std::vector<double>> nodepos;
  double position[3];

  for (Uint i = 0, isize = oft.GetBodyName().size(); i < isize; ++i) {
    nodepos.clear();
    oft.GetNodalPositionByBodyName(oft.GetBodyName(i), nodepos);

    const auto &fe_mesh_ = models_[i]->GetFEMesh();
    fe_mesh_->GetPvtkpnts().resize(Model::num_step_);

    for (Uint s = 0; s < Model::num_step_; ++s) {
      fe_mesh_->GetPvtkpnts(s) = vtkSmartPointer<vtkPoints>::New();
      fe_mesh_->GetPvtkpnts(s)->SetNumberOfPoints(models_[i]->GetNumOfNode());


      for (Uint n = 0; n < models_[i]->GetNumOfNode(); ++n) {
        position[0] = nodepos[s][n * 3];
        position[1] = nodepos[s][n * 3 + 1];
        position[2] = nodepos[s][n * 3 + 2];

        // std::cout << "(" << position[0] << ", " << position[1] << ", " << position[2] << ")" << std::endl;
        fe_mesh_->GetPvtkpnts(s)->SetPoint(n, position);
        
      }

    }

    models_[i]->UpdateDisp(0);

  }

  return (DATA_MODEL | DATA_DISPL);
}

int ControlView::InputModelfiles(std::vector<std::string> &argv) {
  std::vector<std::string> hdf5files;
  std::vector<std::string> simple_out_result;
  std::vector<std::string> modelFiles;
  std::vector<std::string> dispFiles;
  std::vector<std::string> contFiles;
  std::vector<std::string> nodeDataFiles;
  ArgParser(argv, hdf5files, simple_out_result, modelFiles, dispFiles, contFiles, nodeDataFiles);
  
  data_type_ = DATA_RESET;
  if (!hdf5files.empty()) {
    data_type_ = this->ReadHDF5Result(hdf5files[0]);

    for(Uint i=0; i < models_.size(); ++i) {
      models_[i]->GetActor()->GetProperty()->SetEdgeColor(0.0, 0.0, 0.0);
      models_[i]->GetActor()->GetProperty()->EdgeVisibilityOn();

      if (models_[i]->GetFEMesh()->GetUGrid()->GetCellType(0) == 4)
        models_[i]->GetActor()->GetProperty()->SetLineWidth(2.5);
      else
        models_[i]->GetActor()->GetProperty()->SetLineWidth(1.0);
      
      models_[i]->GetActor()->SetScale(1.0);
      renderer_->AddActor(models_[i]->GetActor());
    }

    if (data_type_ & DATA_DISPL) {
      sliderbar_ = CreateOneOf<SliderBar>();
      sliderbar_->SetSliderBar(render_window_interactor_, Model::step_, Model::num_step_, this->play_, this->step_play_);

      currenttimer_ = CreateOneOf<CurrentTimer>();
      currenttimer_->SetTextActor(render_window_);
      renderer_->AddActor2D(currenttimer_->GetTextActor());
    }

  } else if (!simple_out_result.empty()) {
    data_type_ = this->ReadSimpleOutResult(simple_out_result[0]);

    for(Uint i=0; i < models_.size(); ++i) {
      models_[i]->GetActor()->GetProperty()->SetEdgeColor(0.0, 0.0, 0.0);
      models_[i]->GetActor()->GetProperty()->EdgeVisibilityOn();

      if (models_[i]->GetFEMesh()->GetUGrid()->GetCellType(0) == 4)
        models_[i]->GetActor()->GetProperty()->SetLineWidth(2.5);
      else
        models_[i]->GetActor()->GetProperty()->SetLineWidth(1.0);
      
      models_[i]->GetActor()->SetScale(1.0);
      renderer_->AddActor(models_[i]->GetActor());
    }

    if (data_type_ & DATA_DISPL) {
      sliderbar_ = CreateOneOf<SliderBar>();
      sliderbar_->SetSliderBar(render_window_interactor_, Model::step_, Model::num_step_, this->play_, this->step_play_);

      currenttimer_ = CreateOneOf<CurrentTimer>();
      currenttimer_->SetTextActor(render_window_);
      renderer_->AddActor2D(currenttimer_->GetTextActor());      
    }

  } else { // using simple output results or old fasion type input.
    models_.resize(modelFiles.size());
    for (Uint i = 0; i< (Uint)modelFiles.size(); ++i) {
      models_[i] = CreateOneOf<Model>();
      models_[i]->SetOffset(Model::num_nodes_* 3);

      if (modelFiles[i].substr(modelFiles[i].size() - 3).compare(".md") == 0) {
        models_[i]->ReadTxtModel(modelFiles[i]);
      } else if (modelFiles[i].substr(modelFiles[i].size() - 3).compare("txt") == 0) {
        models_[i]->ReadTxtModel(modelFiles[i]);
      } else if (modelFiles[i].substr(modelFiles[i].size() - 3).compare("xml") == 0) {
        models_[i]->ReadXmlModel(modelFiles[i]);
      }
      
      //models_[i]->GetActor()->GetProperty()->SetColor(0.9, 0.9, 0.9);
      //models_[i]->GetActor()->GetProperty()->SetOpacity(1.0);
      models_[i]->GetActor()->GetProperty()->SetEdgeColor(0.0, 0.0, 0.0);
      models_[i]->GetActor()->GetProperty()->EdgeVisibilityOn();

      if (models_[i]->GetFEMesh()->GetUGrid()->GetCellType(0) == 4)
        models_[i]->GetActor()->GetProperty()->SetLineWidth(2.5);
      else
        models_[i]->GetActor()->GetProperty()->SetLineWidth(1.0);
      
      models_[i]->GetActor()->SetScale(1.0);

      renderer_->AddActor(models_[i]->GetActor());
    }

    if (!dispFiles.empty()) {
      data_type_ |= DATA_DISPL;

      this->ReadDispFile(dispFiles);

      sliderbar_ = CreateOneOf<SliderBar>();
      sliderbar_->SetSliderBar(render_window_interactor_, Model::step_, Model::num_step_, this->play_, this->step_play_);

      currenttimer_ = CreateOneOf<CurrentTimer>();
      currenttimer_->SetTextActor(render_window_);
      renderer_->AddActor2D(currenttimer_->GetTextActor());
    }

    if (!contFiles.empty()) {
      data_type_ |= DATA_CONPR;

      ReadContactFile(contFiles[0], contacts_, models_);
      for (auto& pcontact : contacts_) {
        pcontact->InitializeUGrid();

        pcontact->GetPrevActor()->GetProperty()->SetLineWidth(2.0);
        pcontact->GetNextActor()->GetProperty()->SetLineWidth(2.0);
        pcontact->GetPrevActor()->GetProperty()->SetPointSize(6.0);
        pcontact->GetNextActor()->GetProperty()->SetPointSize(6.0);
        pcontact->GetPrevActor()->GetProperty()->EdgeVisibilityOff();
        pcontact->GetNextActor()->GetProperty()->EdgeVisibilityOff();

        renderer_->AddActor(pcontact->GetPrevActor());
        renderer_->AddActor(pcontact->GetNextActor());
      }
    }

    if (!nodeDataFiles.empty()) {
      data_type_ |= DATA_NODVL;
      ReadNodalDataFile(nodeDataFiles[0], models_);
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
	
	if(Model::num_step_!=0){
    
		pCtr->GetCurrentTimer()->GetTextActor()->SetPosition(windowSize[0] - 180, windowSize[1] - 25);

		pCtr->GetSliderbar()->GetSliderRep()->GetPoint1Coordinate()->SetValue(windowSize[0] - 260 - windowSize[0]/4.0, windowSize[1] - 15);
		pCtr->GetSliderbar()->GetSliderRep()->GetPoint2Coordinate()->SetValue(windowSize[0] - 260, windowSize[1] - 15);

		int width = 20, height = 20;
		pCtr->GetSliderbar()->GetPlayButton()->SetButtonSizePosition(width, height, windowSize[0] - 235 - width, windowSize[1] - 5 - height);
		pCtr->GetSliderbar()->GetNextButton()->SetButtonSizePosition(width, height, windowSize[0] - 210 - width, windowSize[1] - 5 - height);
		pCtr->GetSliderbar()->GetPrevButton()->SetButtonSizePosition(width, height, windowSize[0] - 185 - width, windowSize[1] - 5 - height);
	}

	pCtr->GetCommandText()->SetTextSizePosition(10, windowSize[1] - 30, 20);
	for (Uint i = 0; i < pCtr->GetModels().size(); ++i) {
		pCtr->GetCommandTextBodies()[i]->SetTextSizePosition(10, windowSize[1] - 30 - 15 * (i+1), 15);
	}

	if (pCtr->GetLookuptable()) {
		double w = windowSize[0] * 70.0 / 800.0;
		double h = windowSize[1] * 300.0 / 640.0;
		pCtr->GetLookuptable()->GetScalarBar()->SetPosition(windowSize[0] - w - 10, 10);
		pCtr->GetLookuptable()->GetScalarBar()->SetPosition2(w, h);
	}

}

void KeypressCallback(vtkObject* caller, long unsigned int vtkNotUsed(eventId), void* clientData, void* vtkNotUsed(callData))
{

	vtkRenderWindowInteractor *iren =
		static_cast<vtkRenderWindowInteractor*>(caller);

	ControlView* pCtr = static_cast<ControlView*>(clientData);

    if(Model::num_step_!=0){
	    if (strcmp(iren->GetKeySym(), "space") == 0){
	        pCtr->IsPlay() = !pCtr->IsPlay();
			if (pCtr->IsPlay()) {
				pCtr->GetSliderbar()->GetPlayButton()->SetState(1);
				pCtr->IsStepPlay() = false;
			}
			else {
				pCtr->GetSliderbar()->GetPlayButton()->SetState(0);
			}
				
	    }

		if (strcmp(iren->GetKeySym(), "b") == 0) {
			pCtr->IsPlay() = false;
			pCtr->GetSliderbar()->GetPlayButton()->SetState(0);
			Model::step_ = (Model::step_ == Model::num_step_ - 1) ? Model::num_step_ - 1 : Model::step_ + 1;
			pCtr->IsStepPlay() = true;
		}

		if (strcmp(iren->GetKeySym(), "v") == 0) {
			pCtr->IsPlay() = false;
			pCtr->GetSliderbar()->GetPlayButton()->SetState(0);
			Model::step_ = (Model::step_ == 0) ? 0 : Model::step_ - 1;
			pCtr->IsStepPlay() = true;
		}
    }
	if (strcmp(iren->GetKeySym(), "i") == 0) {
		if (pCtr->IsShowMarker())
            pCtr->GetAxesline()->GetAxesActor()->VisibilityOff();
		else
			pCtr->GetAxesline()->GetAxesActor()->VisibilityOn();

		pCtr->IsShowMarker() = !pCtr->IsShowMarker();
	}

	if (strcmp(iren->GetKeySym(), "u") == 0) {
		if (pCtr->IsShowMesh()) {
			for(auto &pmodel: pCtr->GetModels())
				pmodel->GetActor()->GetProperty()->EdgeVisibilityOff();
		}
		else {
			for (auto &pmodel : pCtr->GetModels())
				pmodel->GetActor()->GetProperty()->EdgeVisibilityOn();
		}
			
		pCtr->IsShowMesh() = !pCtr->IsShowMesh();
	}

	if (strcmp(iren->GetKeySym(), "l") == 0) {

		if (pCtr->IsShowLabel()) {
			for(auto& pmodel : pCtr->GetModels())
				pmodel->GetLabelActor()->VisibilityOff();
		}
		else {
			for (auto& pmodel : pCtr->GetModels())
				pmodel->GetLabelActor()->VisibilityOn();
		}

		pCtr->IsShowLabel() = !pCtr->IsShowLabel();
	}

  if (strcmp(iren->GetKeySym(), "K") == 0) {
    pCtr->Reset();
  }

	iren->Render();
}  

}


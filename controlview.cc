#include"controlview.h"

void ControlView::readDispfile(const vector<string> & filename) {
  ifstream infile;
  infile.open(filename[0], ios::in | ios::binary);

  if (infile.is_open()) {
    vector<vector<double>> dispvecCollection;
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

void ControlView::readSimpleOutModel(ifstream &infile, vector<int> &modelinfo, vector<int> &elemlist, vector<double> &nodelist) {
  infile.read((char*)modelinfo.data(), sizeof(int) * 3);

  elemlist.resize(modelinfo[0] * modelinfo[2]);
  nodelist.resize(modelinfo[1] * 3);

  infile.read((char*)elemlist.data(), sizeof(int)* elemlist.size());
  infile.read((char*)nodelist.data(), sizeof(double)* nodelist.size());
}


void ControlView::readSimpleOutResult(const string& filename) {
  ifstream infile(filename);
  if (!infile.is_open()) {
    cout << "Error in opening " << filename << endl;
    exit(1);
  }

  unsigned bodynum;
  infile.read((char*)bodynum, sizeof(bodynum));

  vector<int> modelinfo(3);
  vector<int> elemlist;
  vector<double> nodelist;

  for (unsigned i = 0; i < bodynum; ++i) {
    this->readSimpleOutModel(infile, modelinfo, elemlist, nodelist);
    pModels[i]->CreateModel(modelinfo, elemlist, nodelist);
  }

  vector<>
  infile.read()


}
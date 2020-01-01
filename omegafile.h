#pragma once
#ifndef OMEGA_FILE_H
#define OMEGA_FILE_H

#include"h5cpp.h"
#include<string>
#include<vector>
#include<iostream>


class OmegaFile {
public:
  OmegaFile() {};
  ~OmegaFile() {
    _bodygp.close();
    _consgp.close();
    _solution.close();
    _h5file.close();
  };
  OmegaFile(const std::string &filename) {
    _h5file.openFile(filename, H5F_ACC_RDWR);
    _solution = _h5file.openGroup("solution");
    _bodygp = _solution.openGroup("body");
    _consgp = _solution.openGroup("constraint");
  }

  const unsigned GetBodyNum () const;
  const unsigned GetConstraintNum() const;
  void Initialize();
  
  const std::string GetBodyName(unsigned i) const {
    return _bodyname_list[i];
  }

  const std::vector<std::string> GetBodyName() const {
    return _bodyname_list;
  }

  const std::string GetConsName(unsigned i) const {
    return _consname_list[i];
  }

  const std::vector<std::string> GetConsName() const {
    return _consname_list;
  }

  // basic operation
  void GetDoubleData(H5::DataSet &dset, std::vector<double> &val_out,
                      unsigned row_offset = 0, unsigned col_offset = 0,
                      unsigned row_count = 0, unsigned col_count = 0);

  void GetUnsignedData(H5::DataSet &dset, std::vector<unsigned> &val_out,
                        unsigned row_offset = 0, unsigned col_offset = 0,
                        unsigned row_count = 0, unsigned col_count = 0);

  void GetDataDimension(H5::DataSet &dset, unsigned &rows, unsigned &cols);

  void GetStringAttribute(H5::H5Object &h5obj, const std::string &attrname, std::string &attr_out);
  void GetDoubleAttribute(H5::H5Object &h5obj, const std::string &attrname, double &attr_out);
  void GetUnsignedAttribute(H5::H5Object &h5obj, const std::string &attrname, unsigned &attr_out);

  // special operation
  void GetMeshByBodyName(const std::string &bodyname,
                          std::string &meshtype,
                          std::vector<unsigned> &meshinfo,
                          std::vector<unsigned> &elem,
                          std::vector<double> &node);
  void GetNodalPositionByBodyName(const std::string &bodyname, std::vector<std::vector<double>> &nodepos);
  void GetTimeSeries(std::vector<double> &times);

protected:
  H5::H5File _h5file;
  H5::Group _solution;
  H5::Group _bodygp;
  H5::Group _consgp;

  std::vector<std::string> _bodyname_list;
  std::vector<std::string> _consname_list;


};

#endif


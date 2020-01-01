#include "omegafile.h"

template<typename T>
std::ostream & operator<<(std::ostream &ss, std::vector<T>c) {
  for (unsigned i = 0, isize = c.size(); i < isize; ++i)
    ss << c[i] << " ";
  return ss;
}

herr_t Op(hid_t group, const char *name, void *data) {
  H5G_stat_t statbuf;
  H5Gget_objinfo(group, name, false, &statbuf); //name是每一轮遍历的名称，statbuf存储文件信息，内部的H5G_obj_t可以鉴别是数据集还是group类型。data是外部传入的数据，可以把遍历文件的信息导出去
  std::vector<std::string> *opdata = (std::vector<std::string> *)data;
  opdata->push_back(name);
  return 0;
}

void GenerateRmat(double p[4], double R[3][3]) {
  R[0][0] = 2.0 * (p[0] * p[0] + p[1] * p[1]) - 1.0;
  R[0][1] = 2.0 * (p[1] * p[2] - p[0] * p[3]);
  R[0][2] = 2.0 * (p[1] * p[3] + p[0] * p[2]);

  R[1][0] = 2.0 *(p[1] * p[2] + p[0] * p[3]);
  R[1][1] = 2.0 * (p[0] * p[0] + p[2] * p[2]) - 1.0;
  R[1][2] = 2.0 *(p[2] * p[3] - p[0] * p[1]);

  R[2][0] = 2.0*(p[1] * p[3] - p[0] * p[2]);
  R[2][1] = 2.0*(p[2] * p[3] + p[0] * p[1]);
  R[2][2] = 2.0*(p[0] * p[0] + p[3] * p[3]) - 1.0;
}

void RotateVector(double R[3][3], double vin[3], double vout[3]) {
  vout[0] = R[0][0] * vin[0] + R[0][1] * vin[1] + R[0][2] * vin[2];
  vout[1] = R[1][0] * vin[0] + R[1][1] * vin[1] + R[1][2] * vin[2];
  vout[2] = R[2][0] * vin[0] + R[2][1] * vin[1] + R[2][2] * vin[2];
}

void AddVector(double r[3], double rin[3]) {
  r[0] += rin[0];
  r[1] += rin[1];
  r[2] += rin[2];
}

//############################################################################
// member functions
//############################################################################
const unsigned OmegaFile::GetBodyNum() const {
  return (unsigned)_bodygp.getNumObjs();
}

const unsigned OmegaFile::GetConstraintNum() const {
  return (unsigned)_consgp.getNumObjs();
}

void OmegaFile::Initialize() {
  unsigned int num = GetBodyNum();
  std::vector<int> idbody(num);
  _bodygp.iterateElems("./", idbody.data(), Op, &_bodyname_list);

  num = GetConstraintNum();
  std::vector<int> idcons(num);
  _consgp.iterateElems("./", idcons.data(), Op, &_consname_list);
}

void OmegaFile::GetUnsignedData(H5::DataSet &dset, std::vector<unsigned> &val_out, 
                                unsigned row_offset, unsigned col_offset,
                                unsigned row_count, unsigned col_count) {
  H5::DataSpace val_dspace = dset.getSpace();
  hsize_t dims[2];
  val_dspace.getSimpleExtentDims(dims, NULL);
  hsize_t offset[2] = { row_offset, col_offset };
  if (row_count != 0) dims[0] = row_count;
  if (col_count != 0) dims[1] = col_count;

  H5::DataSpace mspace(2, dims, NULL);
  val_dspace.selectHyperslab(H5S_SELECT_SET, dims, offset);
  val_out.resize(dims[0] * dims[1]);
  dset.read(val_out.data(), H5::PredType::NATIVE_UINT, mspace, val_dspace);
}

void OmegaFile::GetDoubleData(H5::DataSet &dset, std::vector<double> &val_out,
                              unsigned row_offset, unsigned col_offset,
                              unsigned row_count, unsigned col_count) {
  H5::DataSpace val_dspace = dset.getSpace();
  hsize_t dims[2];
  val_dspace.getSimpleExtentDims(dims, NULL);
  hsize_t offset[2] = { row_offset, col_offset };
  if (row_count != 0) dims[0] = row_count;
  if (col_count != 0) dims[1] = col_count;

  H5::DataSpace mspace(2, dims, NULL);
  val_dspace.selectHyperslab(H5S_SELECT_SET, dims, offset);
  val_out.resize(dims[0] * dims[1]);
  dset.read(val_out.data(), H5::PredType::NATIVE_DOUBLE, mspace, val_dspace);
}

void OmegaFile::GetDataDimension(H5::DataSet &dset, unsigned &rows, unsigned &cols) {
  H5::DataSpace val_dspace = dset.getSpace();
  hsize_t dims[2];
  val_dspace.getSimpleExtentDims(dims, NULL);
  rows = dims[0];
  cols = dims[1];
}

void OmegaFile::GetStringAttribute(H5::H5Object & h5obj, const std::string & attrname, std::string & attr_out) {
  H5::Attribute attr = h5obj.openAttribute(attrname);
  H5::DataType type = attr.getDataType();
  attr.read(attr.getDataType(), attr_out);
}

void OmegaFile::GetDoubleAttribute(H5::H5Object & h5obj, const std::string & attrname, double & attr_out) {
  H5::Attribute attr = h5obj.openAttribute(attrname);
  attr.read(attr.getDataType(), &attr_out);
}

void OmegaFile::GetUnsignedAttribute(H5::H5Object & h5obj, const std::string & attrname, unsigned & attr_out) {
  H5::Attribute attr = h5obj.openAttribute(attrname);
  attr.read(attr.getDataType(), &attr_out);
}

void OmegaFile::GetMeshByBodyName(const std::string & bodyname, 
                                  std::string &meshtype,
                                  std::vector<unsigned> &meshinfo, 
                                  std::vector<unsigned> &elem, 
                                  std::vector<double> &node) {
  H5::Group body = _bodygp.openGroup(bodyname);
  H5::Group mesh = body.openGroup("mesh");

  // read mesh attribute
  GetStringAttribute(mesh, "mesh type", meshtype);
  unsigned elemnum, nodenum, nofe;
  GetUnsignedAttribute(mesh, "num elem", elemnum);
  GetUnsignedAttribute(mesh, "num node", nodenum);
  GetUnsignedAttribute(mesh, "num elemental node", nofe);
  meshinfo.push_back(elemnum);
  meshinfo.push_back(nodenum);
  meshinfo.push_back(nofe);

  // read element and node list
  GetUnsignedData(mesh.openDataSet("element"), elem);

  std::vector<double> buffnode;
  GetDoubleData(mesh.openDataSet("node"), buffnode);
  node.resize(buffnode.size());

  unsigned buffidx = 0, buffidy = nodenum, buffidz = 2* nodenum;
  for (unsigned i = 0; i < nodenum; ++i) {
    node[i * 3] = buffnode[buffidx++];
    node[i * 3 + 1] = buffnode[buffidy++];
    node[i * 3 + 2] = buffnode[buffidz++];
  }

}

void OmegaFile::GetNodalPositionByBodyName(const std::string & bodyname, std::vector<std::vector<double>> &nodepos) {
  H5::Group body = _bodygp.openGroup(bodyname);
  std::string bodytype;
  GetStringAttribute(body, "type", bodytype);
  
  H5::DataSet val_dset = body.openDataSet("val");

  std::vector<double> nodepos_buff;
  GetDoubleData(val_dset, nodepos_buff);

  unsigned stepnum, dofs;
  GetDataDimension(val_dset, stepnum, dofs);
  
  std::string meshtype;
  std::vector<unsigned> meshinfo, elem;
  std::vector<double> node;

  nodepos.resize(stepnum);

  if (bodytype.compare("RigidBody") == 0) {
    GetMeshByBodyName(bodyname, meshtype, meshinfo, elem, node);
    double r0[3], r[3], rc[3], p[4], R[3][3];
    for (unsigned s = 0; s < stepnum; ++s) {
      nodepos[s].resize(meshinfo[1]* 3);
      
      for (unsigned i = 0; i < meshinfo[1]; ++i) {
        r0[0] = node[i * 3];
        r0[1] = node[i * 3 + 1];
        r0[2] = node[i * 3 + 2];

        p[0] = nodepos_buff[s* dofs + 3];
        p[1] = nodepos_buff[s* dofs + 4];
        p[2] = nodepos_buff[s* dofs + 5];
        p[3] = nodepos_buff[s* dofs + 6];

        GenerateRmat(p, R);
        RotateVector(R, r0, r);

        rc[0] = nodepos_buff[s* dofs];
        rc[1] = nodepos_buff[s* dofs + 1];
        rc[2] = nodepos_buff[s* dofs + 2];

        AddVector(r, rc);

        nodepos[s][i* 3] = r[0];
        nodepos[s][i* 3 + 1] = r[1];
        nodepos[s][i* 3 + 2] = r[2];

      }
    }
  }

}

void OmegaFile::GetTimeSeries(std::vector<double> &times) {
  H5::DataSet time_series = _solution.openDataSet("time");
  GetDoubleData(time_series, times);
}
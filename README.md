# univiewer
a lightweight FEM viewer API using VTK

## what is it
this code was implemented by VTK and Eigen, it can be a tool to display FEM model.

## who may be interested
FEM or Multibody major researcher

## how to use it
>**syntax**   
> `shared_ptr<Univiewer> uv = Univiewer::New();`   
> `uv->plotModel(elem, node);`   
or   
> `uv->plotModel(meshes);`   
> where `elem = MatrixXu(3, elem_num); node = MatrixXd(3, node_num);`   
> or `meshes = vector<TriangleMesh>(elem, node);`   

Because of the dependency of the code, you should have,  
1. VTK toolkit  
2. Eigen  


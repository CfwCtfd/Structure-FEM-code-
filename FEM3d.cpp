//****************************************************************************80
//
//  Purpose:
//
//    This is an FEM solver on Tetrahedral Mesh with 10 nodes
//
//  Discussion:
//
//    Solves the Elasticity equation for bending problem which arises in aero-elasticity
//    although other load condition can be easily implemented. This code works with the 
//    pre-processor MeshPreProcessor.cpp which in turn uses the c++ code of John Burkardt
//    tet_mesh_l2q.cpp to generate the tet10 elements. The basic mesh is genarated using 
//    pointwise and exported for SU2 mesh format. The mesh is further processed using 
//    to generate three input files for the current program namely:
//    1. element.dat which has the element node information as well as info on 
//                     the faces, normals, areas, volumes
//
//    2. nodes.dat contains the nodes coordinate
//
//    3. BoundaryOrdering.dat which has boundary node information.
//
//    References: The Quadratic Tetrahedron:lecture notes by Carlos Fellipa.
//
//   Modified:
//
//    14 June 2017
//
//   Author: Kaushik K N 
//
//    N.A.L. 


     #include<cstdlib>
     #include<iostream>
     #include<iomanip>
     #include<fstream>
     #include<ctime>
     #include<cmath>
     #include<string>
     #include<sstream>
     #include<vector>
     #include<complex>
     #include<sys/types.h>
     #include<sys/stat.h>
     #include<unistd.h>
     #include<algorithm>
     #include<functional>
     #include <iterator>
     using namespace std;
     #include<Eigen/Core>
     #include<Eigen/Dense>
     #include<Eigen/Eigenvalues> 
     #include<Eigen/SVD>
     #include<Eigen/Sparse>
     #include<Eigen/SparseCholesky>
     #include<Eigen/SparseLU>
    using namespace Eigen;
    typedef Eigen::SparseMatrix<double> SpMat;
    struct NodeInfo
     {
      double x,y,z;
     };
    struct ElementInfo
    {
      int node[10];
      double volume, area[4],nx[4],ny[4],nz[4];
    };
   struct constants
   {
    double YoungsModulus, PoissonRatio, Density, gravity;
   };
   struct Boundarytype
   {
    int Btype;
    int Nbtype;
    std::vector<std::pair<int,int> > ElemFace;
   };

   std::istream& operator>>(std::istream& is, ElementInfo& element)
   {
    is >> element.node[0] >> element.node[1] >> element.node[2] >>element.node[3] 
   >> element.node[4] >> element.node[5] >> element.node[6] >> element.node[7]
   >> element.node[8] >> element.node[9] >> element.volume
   >> element.nx[0] >> element.ny[0] >>element.nz[0] >> element.area[0] 
   >> element.nx[1] >> element.ny[1] >>element.nz[1] >> element.area[1]
   >> element.nx[2] >> element.ny[2] >>element.nz[2] >> element.area[2]
   >> element.nx[3] >> element.ny[3] >>element.nz[3] >> element.area[3];
  return is;
  }
   std::istream& operator>>(std::istream& is, NodeInfo& node)
  {
    is >> node.x >>node.y >> node.z ;
    return is;
  }

  void ComputeElasticityMat(Eigen::MatrixXd& El, constants& M);
  void ComputeElementStiffnessBodyForce(ElementInfo& ithEle, vector<NodeInfo>& node, Eigen::MatrixXd& K,
  Eigen::MatrixXd& E, Eigen::Matrix4d& eta, double weight,constants& M,Eigen::VectorXd& BodyForce);
  void ComputeSurfaceForce(ElementInfo& ithEle,int LocalFace,Eigen::VectorXd& SurfaceForce, 
                            double ElementLoad, int (*FaceNodes)[6]);
  void AssembleElement(ElementInfo& ithEle, Eigen::MatrixXd& KelemMat, Eigen::VectorXd& BodyForce, Eigen::VectorXd& SurfaceForce,
                     Eigen::MatrixXd& KGl, Eigen::VectorXd&  BFGl, int npe, int ndf);
  void PhysicalBoundaryConditions(Eigen::MatrixXd& KGl, Eigen::VectorXd& BFGl, vector<Boundarytype>&  BoundaryOrdering, 
                    vector<ElementInfo> Elem, int boundaryNum, int TotalDOF, int ndf,std::vector<NodeInfo>& node, int (*FaceNodes)[6]);
  void SparseSolve(Eigen::MatrixXd& KGl, Eigen::VectorXd& BFGl, Eigen::VectorXd& wsol);
  void WriteSolution(char* filename, std::vector<ElementInfo>& Elem, std::vector<NodeInfo>& node, Eigen::VectorXd& wsol );
int main()
 {
   char filename[] = "elems_tet10.dat";
   char filename1[] ="nodes10.dat";
   char OutputFile[] = "FEMSolution.dat";
   int nelem, nnodes, TotalDOF;
   int ndf=3;//number of degrees of freedom at each node
   int npe=10;//number of points per element
   int nbtypes, LocalFace;
// Gauss Integration points
   double alfa=(5.0+3.0*sqrt(5.0))/20;
   double beta =(5.0-sqrt(5.0))/20;
   double weight=1.0/4.0;
   Matrix4d GaussPoints;
   GaussPoints<<alfa, beta, beta, beta,
                beta, alfa, beta, beta,
                beta, beta, alfa, beta,
                beta, beta, beta, alfa;
   int FaceNodes[4][6]={{1,2,3,5,8,9},{0,2,3,6,7,9},{0,1,3,4,7,8},{0,1,2,4,5,6}};
// To define the element properties
   constants MaterialProperties={.YoungsModulus=70e9,.PoissonRatio=1.0/3.0, .Density=2830, .gravity=9.81}; 
// Read the data file for the elements and the nodes
   std::vector<ElementInfo> Elem;
   std::vector<NodeInfo> Node;
   std::ifstream ifs(filename), ifs1(filename1);
   ifs>>nelem;
   cout<<"The total number of elements is="<<nelem<<endl;
   if (ifs) {
     std::copy(std::istream_iterator<ElementInfo>(ifs), 
            std::istream_iterator<ElementInfo>(),
            std::back_inserter(Elem));
   }
   else {
     std::cerr << "Couldn't open " << filename << " for reading\n";
   }
     cout<<"Element file read sucessfuly  "<<endl;

   if (ifs1) {
    std::copy(std::istream_iterator<NodeInfo>(ifs1), 
            std::istream_iterator<NodeInfo>(),
            std::back_inserter(Node));
   }
   else {
    std::cerr << "Couldn't open " << filename1 << " for reading\n";
   }
    cout<<"Node file read sucessfuly  "<<endl;

     std::vector<Boundarytype> BoundaryOrdering;
     std::ifstream bound("BoundaryOrdered.dat");
     std::ofstream ofs;

     int nface,type;
     bound>>nbtypes;
     cout<<"The number of boundary types="<<nbtypes<<endl;
     BoundaryOrdering.resize(nbtypes);
     for(int j =0 ; j < nbtypes ; j++ )
     {
	     bound>>type>>nface;
             if(bound.eof()){break;}
	     cout<<nface<<"    "<<type<<endl;
            
	     BoundaryOrdering[type].Nbtype=nface;
             BoundaryOrdering[type].Btype=type;
             BoundaryOrdering[type].ElemFace.resize(nface);
	     for(std::vector<std::pair<int,int> >::iterator it= BoundaryOrdering[type].ElemFace.begin(); it!= BoundaryOrdering[type].ElemFace.end();it++)
		     bound>>it->first>>it->second;
        
     }	     
     
  ofs.close();
  bound.close();   


  cout<<"The total number of Nodes is="<<Node.size() <<endl;
  nnodes=Node.size();

  //To form the Elasticity Matrix
   TotalDOF=nnodes*ndf;
  cout<<"The Total number of Degrees of freedom=:\n"<<TotalDOF<<endl;


   MatrixXd EelemMat(6,6), KelemMat(30,30),  KGl(TotalDOF, TotalDOF);
   VectorXd BodyForce(30), BFGl(TotalDOF), SurfaceForce(30);
   EelemMat.setZero(6,6);

   ComputeElasticityMat(EelemMat,MaterialProperties);
   cout<<"Elasticity Matrix=\n"<<EelemMat<<endl;
   cout<<"Elem size="<<Elem.size()<<endl;
   //To find the Stiffness Matrix
    if(nelem!=Elem.size()){
    cout<<"Something wrong in element file...please check!"<<endl;
    exit(0);
   }
  double ElementLoad = -1750; 
  VectorXd wsol(TotalDOF);

   
//0:Back 1:Bottom 2:Front 3:Left 4:Right 5:Top

 for(int i=0;i<nelem;i++)
   {
     for(int k=0;k<BodyForce.size();k++){
        BodyForce(k)=0.0;
        SurfaceForce(k)=0.0;}
        KelemMat.setZero(30,30);
        ComputeElementStiffnessBodyForce(Elem[i],Node,KelemMat,
        EelemMat,GaussPoints,weight,MaterialProperties,BodyForce);
  
             for(std::vector<std::pair<int,int> >::iterator it= BoundaryOrdering[5].ElemFace.begin(); it!= BoundaryOrdering[5].ElemFace.end();it++)
                {	
                  if(it->first==i){
                   LocalFace=it->second;
                   ComputeSurfaceForce(Elem[i],LocalFace, SurfaceForce,ElementLoad,FaceNodes);
                  } 
                }
       for(int k=0;k<BodyForce.size();k++) BodyForce(k)=0.0;
        AssembleElement(Elem[i],KelemMat,BodyForce,SurfaceForce, KGl,BFGl,npe,ndf);  
    }
 
   cout<<"Finished Assembling the element"<<endl;
   PhysicalBoundaryConditions(KGl,BFGl,BoundaryOrdering,Elem,3,TotalDOF,ndf, Node, FaceNodes);
   cout<<"Finished applying the physical boundary conditions"<<endl;
  
   cout<<"Solving the system...."<<endl;
   SparseSolve(KGl,BFGl,wsol);
   cout<<"Finished solving the force system"<<endl;

   cout<<"Writing the solution...."<<endl;
   WriteSolution(OutputFile, Elem, Node, wsol);
   cout<<"Finished writing the solution"<<endl;
   return 0;    
 } 

 void ComputeSurfaceForce(ElementInfo& ithEle,int LocalFace,Eigen::VectorXd& SurfaceForce, double ElementLoad, int (*FaceNodes)[6])
  {
   double weight=1.0/3.0;
   double SurfaceArea;
   Matrix3d GaussPoints;
   Vector3d b,xi;
   MatrixXd N(3,30);
   VectorXd shape(6);
   GaussPoints << 0.5, 0.5, 0.0,
                  0.0, 0.5, 0.5,
                  0.5, 0.0, 0.5;
   Vector3d temp;
   int stride=3;
   double Normal;

   b(0)=ElementLoad*ithEle.nx[LocalFace];
   b(1)=ElementLoad*ithEle.ny[LocalFace];
   b(2)=ElementLoad*ithEle.nz[LocalFace];
   Normal=sqrt(pow(ithEle.nx[LocalFace],2)+pow(ithEle.ny[LocalFace],2)+pow(ithEle.nz[LocalFace],2));
   for(int intPoints=0;intPoints<3;intPoints++)
    {
     for(int j=0;j<3;j++)
       xi(j)=GaussPoints(intPoints,j);
        N.setZero(3,30);
         for(int k=0;k<3;k++)
          shape(k)=xi(k)*(2.0*xi(k)-1.0); 
          shape(3)=4.0*xi(0)*xi(1);
          shape(4)=4.0*xi(1)*xi(2);
          shape(5)=4.0*xi(0)*xi(2);

          for(int m=0;m<6;m++){
             temp<<shape(m),shape(m),shape(m);
             N.middleCols(FaceNodes[LocalFace][m]*stride,stride)=temp.asDiagonal();       
           }
                 
           SurfaceForce += weight*(N.transpose()*b)*Normal*ithEle.area[LocalFace];
    }
          return;
  }
  

 void ComputeElasticityMat(Eigen::MatrixXd& El, constants& M)
  {
   double K=M.YoungsModulus/((1.0+M.PoissonRatio)*  (1.0-2.0*M.PoissonRatio));
    El(0,0)=1.0-M.PoissonRatio;
    El(0,1)=M.PoissonRatio;
    El(0,2)=M.PoissonRatio;

    El(1,0)=M.PoissonRatio;
    El(1,1)=1.0-M.PoissonRatio;
    El(1,2)=M.PoissonRatio;

    El(2,0)=M.PoissonRatio;
    El(2,1)=M.PoissonRatio;
    El(2,2)=1.0-M.PoissonRatio;

    El(3,3)=0.5-M.PoissonRatio;
    El(4,4)=El(3,3);
    El(5,5)=El(3,3);

    El=El*K;

    return;
  }

 void ComputeElementStiffnessBodyForce(ElementInfo& ithEle, vector<NodeInfo>& node, Eigen::MatrixXd& K,
  Eigen::MatrixXd& E, Eigen::Matrix4d& eta, double weight,constants& M,Eigen::VectorXd& BodyForce)
{
  double jx1,jx2,jx3,jx4,jy1,jy2,jy3,jy4,jz1,jz2,jz3,jz4,
       a1,a2,a3,a4,b1,b2,b3,b4,c1,c2,c3,c4,Jdet;
  
  VectorXd Nfx(10),Nfy(10),Nfz(10),xi(4), shape(10);
  MatrixXd B(6,30), Bt; 
  MatrixXd J(4,4), P, Iaug(4,3),Jinv(4,4);
  double v01,v02,v03,v04,V;
  MatrixXd N(3,30);
  Vector3d b;
  double test_int=0.0;
  b<<0,0,-M.Density*M.gravity;
  for(int intPoints=0;intPoints<4;intPoints++)
  {
    N.setZero(3,30);
     for(int j=0;j<4;j++)
       xi(j) = eta(intPoints,j);

       jx1=4.0*(node[ithEle.node[0]].x*(xi(0)-0.25)+node[ithEle.node[4]].x*xi(1)+node[ithEle.node[6]].x*xi(2)+node[ithEle.node[7]].x*xi(3));
       jy1=4.0*(node[ithEle.node[0]].y*(xi(0)-0.25)+node[ithEle.node[4]].y*xi(1)+node[ithEle.node[6]].y*xi(2)+node[ithEle.node[7]].y*xi(3));
       jz1=4.0*(node[ithEle.node[0]].z*(xi(0)-0.25)+node[ithEle.node[4]].z*xi(1)+node[ithEle.node[6]].z*xi(2)+node[ithEle.node[7]].z*xi(3));

       jx2=4.0*(node[ithEle.node[4]].x*xi(0)+node[ithEle.node[1]].x*(xi(1)-0.25)+node[ithEle.node[5]].x*xi(2)+node[ithEle.node[8]].x*xi(3));
       jy2=4.0*(node[ithEle.node[4]].y*xi(0)+node[ithEle.node[1]].y*(xi(1)-0.25)+node[ithEle.node[5]].y*xi(2)+node[ithEle.node[8]].y*xi(3));
       jz2=4.0*(node[ithEle.node[4]].z*xi(0)+node[ithEle.node[1]].z*(xi(1)-0.25)+node[ithEle.node[5]].z*xi(2)+node[ithEle.node[8]].z*xi(3));

       jx3=4.0*(node[ithEle.node[6]].x*xi(0)+node[ithEle.node[5]].x*xi(1)+node[ithEle.node[2]].x*(xi(2)-0.25)+node[ithEle.node[9]].x*xi(3));
       jy3=4.0*(node[ithEle.node[6]].y*xi(0)+node[ithEle.node[5]].y*xi(1)+node[ithEle.node[2]].y*(xi(2)-0.25)+node[ithEle.node[9]].y*xi(3));
       jz3=4.0*(node[ithEle.node[6]].z*xi(0)+node[ithEle.node[5]].z*xi(1)+node[ithEle.node[2]].z*(xi(2)-0.25)+node[ithEle.node[9]].z*xi(3));

       jx4=4.0*(node[ithEle.node[7]].x*xi(0)+node[ithEle.node[8]].x*xi(1)+node[ithEle.node[9]].x*xi(2)+node[ithEle.node[3]].x*(xi(3)-0.25));
       jy4=4.0*(node[ithEle.node[7]].y*xi(0)+node[ithEle.node[8]].y*xi(1)+node[ithEle.node[9]].y*xi(2)+node[ithEle.node[3]].y*(xi(3)-0.25));
       jz4=4.0*(node[ithEle.node[7]].z*xi(0)+node[ithEle.node[8]].z*xi(1)+node[ithEle.node[9]].z*xi(2)+node[ithEle.node[3]].z*(xi(3)-0.25));

       J.row(0)<<1,1,1,1;
       J.row(1)<<jx1,jx2,jx3,jx4;
       J.row(2)<<jy1,jy2,jy3,jy4;
       J.row(3)<<jz1,jz2,jz3,jz4;

       Jdet=J.determinant();
       Jinv=J.inverse();

       Iaug<<0,0,0,
             1,0,0,
             0,1,0,
             0,0,1;

       P=Jinv*Iaug;

       a1=P(0,0);
       a2=P(1,0);
       a3=P(2,0);
       a4=P(3,0);

       b1=P(0,1);
       b2=P(1,1);
       b3=P(2,1);
       b4=P(3,1);

       c1=P(0,2);
       c2=P(1,2);
       c3=P(2,2);
       c4=P(3,2);

       Nfx << (4.0*xi(0)-1)*a1, (4.0*xi(1)-1)*a2, (4.0*xi(2)-1)*a3, (4.0*xi(3)-1)*a4,
              4.0*(a1*xi(1)+a2*xi(0)), 4.0*(a2*xi(2)+a3*xi(1)), 4.0*(a1*xi(2)+a3*xi(0)),
              4.0*(a1*xi(3)+a4*xi(0)), 4.0*(a2*xi(3)+a4*xi(1)), 4.0*(a3*xi(3)+a4*xi(2));

       Nfy << (4.0*xi(0)-1)*b1, (4.0*xi(1)-1)*b2, (4.0*xi(2)-1)*b3, (4.0*xi(3)-1)*b4,
              4.0*(b1*xi(1)+b2*xi(0)), 4.0*(b2*xi(2)+b3*xi(1)), 4.0*(b1*xi(2)+b3*xi(0)),
              4.0*(b1*xi(3)+b4*xi(0)), 4.0*(b2*xi(3)+b4*xi(1)), 4.0*(b3*xi(3)+b4*xi(2));

       Nfz << (4.0*xi(0)-1)*c1, (4.0*xi(1)-1)*c2, (4.0*xi(2)-1)*c3, (4.0*xi(3)-1)*c4,
              4.0*(c1*xi(1)+c2*xi(0)), 4.0*(c2*xi(2)+c3*xi(1)), 4.0*(c1*xi(2)+c3*xi(0)),
              4.0*(c1*xi(3)+c4*xi(0)), 4.0*(c2*xi(3)+c4*xi(1)), 4.0*(c3*xi(3)+c4*xi(2));          

       B.row(0) << Nfx(0),0,0, Nfx(1),0,0, Nfx(2),0,0, Nfx(3),0,0, Nfx(4),0,0, Nfx(5),0,0, Nfx(6),0,0, Nfx(7),0,0, Nfx(8),0,0, Nfx(9),0,0;
       B.row(1) << 0,Nfy(0),0, 0,Nfy(1),0, 0,Nfy(2),0, 0,Nfy(3),0, 0,Nfy(4),0, 0,Nfy(5),0, 0,Nfy(6),0, 0,Nfy(7),0, 0,Nfy(8),0, 0,Nfy(9),0;
       B.row(2) << 0,0,Nfz(0), 0,0,Nfz(1), 0,0,Nfz(2), 0,0,Nfz(3), 0,0,Nfz(4), 0,0,Nfz(5), 0,0,Nfz(6), 0,0,Nfz(7), 0,0,Nfz(8), 0,0,Nfz(9);
       B.row(3) << Nfy(0),Nfx(0),0,  Nfy(1),Nfx(1),0, Nfy(2),Nfx(2),0, Nfy(3),Nfx(3),0, Nfy(4),Nfx(4),0, Nfy(5),Nfx(5),0, Nfy(6),Nfx(6),0,
                   Nfy(7),Nfx(7),0,Nfy(8),Nfx(8),0,  Nfy(9),Nfx(9),0;
       B.row(4) << 0,Nfz(0),Nfy(0), 0,Nfz(1),Nfy(1), 0,Nfz(2),Nfy(2), 0,Nfz(3),Nfy(3), 0,Nfz(4),Nfy(4), 0,Nfz(5),Nfy(5), 0,Nfz(6),Nfy(6),
                   0,Nfz(7),Nfy(7), 0,Nfz(8),Nfy(8), 0,Nfz(9),Nfy(9);
       B.row(5) << Nfz(0),0,Nfx(0), Nfz(1),0,Nfx(1), Nfz(2),0,Nfx(2), Nfz(3),0,Nfx(3), Nfz(4),0,Nfx(4), Nfz(5),0,Nfx(5), Nfz(6),0,Nfx(6),
                   Nfz(7),0,Nfx(7), Nfz(8),0,Nfx(8), Nfz(9),0,Nfx(9);

       Bt = B.transpose();
       K += weight*(Bt*E*B*Jdet/6.0); 
       for(int i=0;i<4;i++)
         shape(i) = xi(i)*(2.0*xi(i)-1.0);
         shape(4) = 4.0*xi(0)*xi(1);
         shape(5) = 4.0*xi(1)*xi(2);
         shape(6) = 4.0*xi(0)*xi(2);
         shape(7) = 4.0*xi(0)*xi(3);
         shape(8) = 4.0*xi(1)*xi(3);
         shape(9) = 4.0*xi(2)*xi(3);  

 
       N.row(0) << shape(0),0,0,shape(1),0,0,shape(2),0,0,shape(3),0,0,shape(4),0,0,shape(5),0,0,shape(6),0,0,shape(7),0,0,shape(8),0,0,shape(9),0,0;
       N.row(1) << 0,shape(0),0,0,shape(1),0,0,shape(2),0,0,shape(3),0,0,shape(4),0,0,shape(5),0,0,shape(6),0,0,shape(7),0,0,shape(8),0,0,shape(9),0;
       N.row(2) << 0,0,shape(0),0,0,shape(1),0,0,shape(2),0,0,shape(3),0,0,shape(4),0,0,shape(5),0,0,shape(6),0,0,shape(7),0,0,shape(8),0,0,shape(9); 


       BodyForce+=weight*N.transpose()*b*Jdet/6.0;
       test_int+=weight*Jdet/6.0;

 }
   return;
}

void AssembleElement(ElementInfo& ithEle, Eigen::MatrixXd& KelemMat, Eigen::VectorXd& BodyForce, Eigen::VectorXd& SurfaceForce, 
                     Eigen::MatrixXd& KGl, Eigen::VectorXd&  BFGl, int npe, int ndf)
{
 for(int m=0;m<npe;m++)
 {
    int nr=(ithEle.node[m])*ndf;
    for(int ii=0;ii<ndf;ii++)
    {
     nr+=1;
     int l =m*ndf+ii;
       BFGl(nr-1) += BodyForce(l)+SurfaceForce(l);
 
        for(int p=0;p<npe;p++)
         {
         int ncl=(ithEle.node[p])*ndf;
          for(int kk=0;kk<ndf;kk++)
            {
               ncl+=1;
              int q=p*ndf+kk;
                KGl(nr-1,ncl-1)=KGl(nr-1,ncl-1)+KelemMat(l,q);
               

            }

          }

      }
 
   }

return;
}

void PhysicalBoundaryConditions(Eigen::MatrixXd& KGl, Eigen::VectorXd& BFGl, vector<Boundarytype>&  BoundaryOrdering, 
                    vector<ElementInfo> Elem, int boundaryNum, int TotalDOF, int ndf,std::vector<NodeInfo>& node, int (*FaceNodes)[6])
{
	int igl,node_num,m;
            
     	for(std::vector<std::pair<int,int> >::iterator it= BoundaryOrdering[boundaryNum].ElemFace.begin(); it!= BoundaryOrdering[boundaryNum].ElemFace.end();it++)
        {	
               
		for(int i=0;i<6;i++) //6 is the total number of face nodes
		{
	               node_num=Elem[it->first].node[FaceNodes[it->second][i]];
                       igl=node_num*ndf;
                        for(int j=0;j<ndf;j++)
			{
				m=igl+j;
				KGl(m,m)=1.0;
				BFGl(m)=0.0;
			
			 for(int k=0;k<TotalDOF;k++)
                          {
			     if(k!=m){
				KGl(k,m)=0.0;
				KGl(m,k)=0.0;
				}
			  }
			    
		       }
                   
                    
		}
     
	}
	
 return;
}

void SparseSolve(Eigen::MatrixXd& KGl, Eigen::VectorXd& BFGl, Eigen::VectorXd& wsol)
{
 SpMat SparseKGl;
 SparseKGl = KGl.sparseView();
 BiCGSTAB<SparseMatrix<double>, IncompleteLUT<double> > bicg; 
 wsol=bicg.compute(SparseKGl).solve(BFGl);	
 std::cout << "#iterations:     " << bicg.iterations() << std::endl;
 std::cout << "estimated error: " << bicg.error()      << std::endl;
// update b, and solve again
// wsol = bicg.solve(BFGl);
 cout<<"Size of wsol is="<<wsol.size()<<endl;
 cout<<"The maximum value of solution="<<wsol.maxCoeff()<<endl;
return;
}

void WriteSolution(char* filename, std::vector<ElementInfo>& Elem, std::vector<NodeInfo>& node, Eigen::VectorXd& wsol )
{
std::ofstream ofs;
ofs.open(filename);
     ofs<<"Title=Solution from the FEM Solver"<<endl;
     ofs<<"VARIABLES = \"X\", \"Y\", \"Z\", \"u\",\"v\",\"w\""<<endl;
     ofs<<"ZONE"<<"  "<<"N="<<node.size()<<" "<<"E="<<Elem.size()<<"  "<<"Datapacking=Point"<<"  "<<"ZONETYPE=FETETRAHEDRON"
           <<endl;
   int count=0;
   for(int i=0;i<node.size();i++){
   ofs<< std::setprecision(3)<<node[i].x<<"\t"<<node[i].y<<"\t"<<node[i].z<<"\t"<<wsol(count)<<"\t"<<wsol(count+1)<<"\t"<<wsol(count+2)<<endl;
    count+=3;
    }
   ofs<<"\n";
   for(int i=0;i<Elem.size();i++)
   ofs<<std::setprecision(3)<<Elem[i].node[0]+1<<"\t"<<Elem[i].node[1]+1<<"\t"<<Elem[i].node[2]+1<<"\t"<<Elem[i].node[3]+1<<endl;
   ofs.close();
return;
}

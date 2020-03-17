//****************************************************************************80
//
//  Purpose:
//
//    This is the pre-processor for the fem solver on Tetrahedral Mesh with 10 nodes
//
//  Discussion:
//
//    The basic mesh is genarated using pointwise and exported for SU2 mesh format. 
//
//    10 noded tet is generated using code of John Burkardt: tet_mesh_l2q.cpp
//
//    The code generates three files for the structural solver namely
//
//    1. element.dat which has the element node information as well as info on 
//                     the faces, normals, areas, volumes
//
//    2. nodes.dat contains the nodes coordinate
//
//    3. BoundaryOrdering.dat which has boundary node information.
//
//    Modified:
//
//    14 June 2017
//
//   Author:
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
#include <fstream>

using namespace std;
struct Node
{
  double x,y,z;
};
struct Element
{
   int node[4];
};
struct BoundaryNodes
{
    int nodeB[3];
};
	
struct Boundary
{
	int nfaces;
	std::vector<BoundaryNodes> Face;
};
struct ElementInfo
{
     
    double volume, area[4],nx[4],ny[4],nz[4];
};
struct Boundarytype
{
int Btype;
int Nbtype;
std::vector<std::pair<int,int> > ElemFace;
};
struct ElementInfoTet10
{
    int node[10];
    double volume, area[4],nx[4],ny[4],nz[4];
};
struct ElementsTet10
{
   int node[10];
};
struct NodeInfo
{
 double x,y,z;
};
std::istream& operator>>(std::istream& is, NodeInfo& node)
{
    is >> node.x >>node.y >> node.z ;
    return is;
}

 void ComputeElementMetric(int, Element*, Node*, ElementInfo*);
int main()
{
	int ndim,nelem,nbtypes,nnodes,nface;
        int eleType, eleNum,junk;
	char filename[] = "beam3d.su2";
     	std::string input; 
	std::ifstream ifs(filename);
        std::istringstream line(input);
        std::vector<Boundary> boundary;
      	std::getline(ifs, input, '=');
        cout<<input<<endl;
        ifs>>ndim;
        cout<<"Dimension="<<ndim<<endl;
	std::getline(ifs, input, '=');
        cout<<input<<endl;
        ifs>>nelem;
        cout<<"Number of elements="<<nelem<<endl;
	Element *element;
	element= new Element [nelem];
	for(int i=0;i<nelem;i++)
	ifs>>eleType>>element[i].node[0] >> element[i].node[1]>>element[i].node[2]>>element[i].node[3]>>eleNum;
        
         std::getline(ifs, input, '=');
        ifs>>nnodes;
        cout<<"Number of Nodes="<<nnodes<<endl;
        Node *node;
        node=new Node[nnodes];
        for(int i=0;i<nnodes;i++)	
	ifs>>node[i].x>>node[i].y>>node[i].z>>junk;
     
        
	
        std::getline(ifs, input, '=');
        ifs>>nbtypes;
        cout<<"Number of boundary types= "<<nbtypes<<endl;
        boundary.resize(nbtypes);
        for(int j=0;j<nbtypes;j++)
          {
       	    std::getline(ifs, input, '=');
            cout<<input<<endl;
            std::getline(ifs, input, '=');
            cout<<input<<endl;
    	    ifs>>nface;
            boundary[j].nfaces=nface;
            cout<<j<<"th Boundary nfaces="<<nface<<endl;
            boundary[j].Face.resize(nface);
            for(std::vector<BoundaryNodes>::iterator it= boundary[j].Face.begin(); it!= boundary[j].Face.end();it++)
                    ifs>>junk>>it->nodeB[0]>>it->nodeB[1]>>it->nodeB[2];
          }	
         ElementInfo *elementInfo;
         elementInfo = new ElementInfo[nelem];        
           for(int i=0;i<nelem;i++)
             ComputeElementMetric(i, element,  node, elementInfo);
         
         std::ofstream ofs;
         ofs.open("elems.dat");
         ofs<<nelem<<endl;    
          for(int i=0;i<nelem;i++)
          {
            ofs << std::setprecision(10)<<element[i].node[0] <<"   "<<element[i].node[1] <<"   "<< element[i].node[2] << "   "<<element[i].node[3] << "   "<<elementInfo[i].volume <<"\n"
            << elementInfo[i].nx[0] <<"   "<< elementInfo[i].ny[0] <<"   "<<elementInfo[i].nz[0] <<"   "<< elementInfo[i].area[0] <<"\n"
            << elementInfo[i].nx[1] <<"   "<< elementInfo[i].ny[1] <<"   "<<elementInfo[i].nz[1] <<"   "<< elementInfo[i].area[1] <<"\n"
            << elementInfo[i].nx[2] <<"   "<< elementInfo[i].ny[2] <<"   "<<elementInfo[i].nz[2] <<"   "<< elementInfo[i].area[2] <<"\n"
            << elementInfo[i].nx[3] <<"   "<< elementInfo[i].ny[3] <<"   "<<elementInfo[i].nz[3] <<"   "<< elementInfo[i].area[3]<<endl;
          }
         ofs.close();


        ofs.open("nodes.dat");
          ofs<<nnodes<<endl;
          for(int i=0;i<nnodes;i++)
           ofs<<std::setprecision(20)<<node[i].x<<"\t"<< node[i].y<<"\t"<<node[i].z<<endl;
        ofs.close();

  //To search the element corresponding to a boundary face

      ofs.open("elem_bdryTest.dat");
      for(int j=0;j<nbtypes;j++)
      {
      ofs<<j<<"\t"<<boundary[j].nfaces<<endl;

         for(std::vector<BoundaryNodes>::iterator it= boundary[j].Face.begin(); it!= boundary[j].Face.end();it++)
		     ofs<<it->nodeB[0]<<"\t"<<it->nodeB[1]<<"\t"<<it->nodeB[2]<<endl;		     
         
      
      }
      ofs.close();

     ofstream outfile,outfile1;
     int matchingNodes,localFace, matchingElem;
     outfile.open("BoundaryInfo.dat");
   
     for(int j =0 ; j < nbtypes ; j++ )
     {
         outfile<<j<<"\t"<<boundary[j].nfaces<<endl;
   	for(std::vector<BoundaryNodes>::iterator it= boundary[j].Face.begin(); it!= boundary[j].Face.end();it++)  
	{
	
		for (int itelem=0; itelem < nelem; ++itelem)
		{

             
			matchingNodes=0;
	           for(int i=0;i<3;i++)
		           {
		           if(it->nodeB[i]==element[itelem].node[0] || it->nodeB[i]==element[itelem].node[1]  
                             || it->nodeB[i]==element[itelem].node[2] || it->nodeB[i]==element[itelem].node[3]  )

				matchingNodes+=1;
                      
   		           }
           
                      if(matchingNodes==3)
                       {
			     matchingElem = itelem;
                           for(int k=0;k<4;k++)
			    {
		               if(it->nodeB[0]!=element[itelem].node[k] && it->nodeB[1]!=element[itelem].node[k]  && it->nodeB[2]!=element[itelem].node[k])
			        localFace=k;
                                                      
			    }
			     
		       }
              
	       
           }
      

              outfile<<matchingElem<<"\t"<<localFace<<endl;
     
     }
     cout<<"Completed the matching for the "<<j<<"Boundary"<<endl;

}

     outfile.close();
     std::vector<Boundarytype> BoundaryOrdering;
     std::ifstream bound("BoundaryInfo.dat");
     ofs.open("BoundaryOrdered.dat");
     BoundaryOrdering.resize(nbtypes);
     int type;
     ofs<<nbtypes<<endl;
     for(int j =0 ; j < nbtypes ; j++ )
       {
	     bound>>type>>nface;
             if(bound.eof()){break;}
	     cout<<nface<<"    "<<type<<endl;
            
	     BoundaryOrdering[type].Nbtype=nface;
             BoundaryOrdering[type].Btype=type;
             BoundaryOrdering[type].ElemFace.resize(nface);
	     for(std::vector<std::pair<int,int> >::iterator it= BoundaryOrdering[type].ElemFace.begin(); it!= BoundaryOrdering[type].ElemFace.end();it++){
		     bound>>it->first>>it->second;
               }
               std::sort(BoundaryOrdering[j].ElemFace.begin(), BoundaryOrdering[j].ElemFace.end());
               
               ofs<<type<<"    "<<nface<<endl;
               for(std::vector<std::pair<int,int> >::iterator it= BoundaryOrdering[type].ElemFace.begin(); it!= BoundaryOrdering[type].ElemFace.end();it++){
               ofs<<it->first<<"\t"<<it->second<<endl;

             }
            

       }	

        ofs.close();

        ofs.open("mesh_elements.txt");
        for(int i=0;i<nelem;i++)
        ofs << std::setprecision(15)<<element[i].node[0] <<"      "<<element[i].node[1] <<"      "<< element[i].node[2] << "      "<<element[i].node[3] << endl;
        ofs.close();         

        
        ofs.open("mesh_nodes.txt");
        for(int i=0;i<nnodes;i++)
        ofs<<std::scientific<<node[i].x<<"     "<< node[i].y<<"      "<<node[i].z<<endl;
        ofs.close();


        if (system(NULL)) puts ("Ok");
        else exit (EXIT_FAILURE);
        cout<<"Executing tet10 command...\n"<<endl;
        int  error=system ("./tet10 <input");
        cout<<"The value returned was"<<"   "<<error<<endl;
        system("cp mesh_l2q_elements.txt mesh_l2q_elements.dat");
        system("cp mesh_l2q_nodes.txt mesh_l2q_nodes.dat");
        system("cp mesh_l2q_nodes.txt nodes10.dat");
        ElementsTet10 *elements10;
        elements10=new ElementsTet10[nelem];
        std::ifstream ifst;
        ifst.open("mesh_l2q_elements.dat");
        for(int i=0;i<nelem;i++)
        {
          for(int j=0;j<10;j++ )
          ifst>>elements10[i].node[j];
        }
        ifst.close();

        ElementInfoTet10 *elementsTet10;
        elementsTet10 = new ElementInfoTet10[nelem];
        
        ofs.open("elems_tet10.dat");
        ofs<<nelem<<endl;
        for(int i=0;i<nelem;i++)
        {
         ofs<<std::setprecision(4)<<elements10[i].node[0]<<"   "<<elements10[i].node[1]<<"   "<<elements10[i].node[2]<<"   "<<elements10[i].node[3]<<"   "<<
            elements10[i].node[4]<<"   "<<elements10[i].node[7]<<"   "<<elements10[i].node[5]<<"   "<<elements10[i].node[6]<<"   "<<elements10[i].node[8]
             <<"   "<<elements10[i].node[9]<<"\n"<<elementInfo[i].volume<<"   "<<elementInfo[i].nx[0] <<
               "   "<< elementInfo[i].ny[0] <<"   "<<elementInfo[i].nz[0] <<"   "<< elementInfo[i].area[0] <<"\n"
            << elementInfo[i].nx[1] <<"   "<< elementInfo[i].ny[1] <<"   "<<elementInfo[i].nz[1] <<"   "<< elementInfo[i].area[1] <<"\n"
            << elementInfo[i].nx[2] <<"   "<< elementInfo[i].ny[2] <<"   "<<elementInfo[i].nz[2] <<"   "<< elementInfo[i].area[2] <<"\n"
            << elementInfo[i].nx[3] <<"   "<< elementInfo[i].ny[3] <<"   "<<elementInfo[i].nz[3] <<"   "<< elementInfo[i].area[3]<<endl;
        }
        ofs.close();

      
        std::vector<NodeInfo> node10;
        ifst.open("mesh_l2q_nodes.dat");
        if (ifst) {
        std::copy(std::istream_iterator<NodeInfo>(ifst), 
                std::istream_iterator<NodeInfo>(),
                std::back_inserter(node10));
         }
            else {
          std::cerr << "Couldn't open the file for reading\n";
          }
        ifst.close();
 
        cout<<"The number of nodes in the new mesh="<<node10.size()<<endl;
        return 0;
}
 void ComputeElementMetric(int i, Element* element, Node* node, ElementInfo* elementInfo)
  {
        double x12,x21,x13,x31,x14,x41,x23,x32,x24,x42,x34,x43,
               y12,y21,y13,y31,y14,y41,y23,y32,y24,y42,y34,y43, 
               z12,z21,z13,z31,z14,z41,z23,z32,z24,z42,z34,z43;
  
        double Jdet,a[4],b[4],c[4];
         
        x12=node[element[i].node[0]].x-node[element[i].node[1]].x;
        x21=-x12;
        x13=node[element[i].node[0]].x-node[element[i].node[2]].x;
        x31=-x13;
        x14=node[element[i].node[0]].x-node[element[i].node[3]].x;
        x41=-x14;
        x23=node[element[i].node[1]].x-node[element[i].node[2]].x;  
        x32=-x23;
        x24=node[element[i].node[1]].x-node[element[i].node[3]].x; 
        x42=-x24;
        x34=node[element[i].node[2]].x-node[element[i].node[3]].x;
        x43=-x34;

        y12=node[element[i].node[0]].y-node[element[i].node[1]].y;
        y21=-y12;
        y13=node[element[i].node[0]].y-node[element[i].node[2]].y;
        y31=-y13;
        y14=node[element[i].node[0]].y-node[element[i].node[3]].y;
        y41=-y14;
        y23=node[element[i].node[1]].y-node[element[i].node[2]].y;  
        y32=-y23;
        y24=node[element[i].node[1]].y-node[element[i].node[3]].y; 
        y42=-y24;
        y34=node[element[i].node[2]].y-node[element[i].node[3]].y;
        y43=-y34;

        z12=node[element[i].node[0]].z-node[element[i].node[1]].z;
        z21=-z12;
        z13=node[element[i].node[0]].z-node[element[i].node[2]].z;
        z31=-z13;
        z14=node[element[i].node[0]].z-node[element[i].node[3]].z;
        z41=-z14;
        z23=node[element[i].node[1]].z-node[element[i].node[2]].z;  
        z32=-z23;
        z24=node[element[i].node[1]].z-node[element[i].node[3]].z; 
        z42=-z24;
        z34=node[element[i].node[2]].z-node[element[i].node[3]].z;
        z43=-z34;
  
        Jdet=x21*(y23*z34-y34*z23)+x32*(y34*z12-y12*z34)+x43*(y12*z23-y23*z12);
        elementInfo[i].volume= Jdet/6.0;

        a[0]=y42*z32-y32*z42;
        a[1]=y31*z43-y34*z13;
        a[2]=y24*z14-y14*z24;
        a[3]=y13*z21-y12*z31;

        b[0]=x32*z42-x42*z32;
        b[1]=x43*z31-x13*z34;
        b[2]=x14*z24-x24*z14;
        b[3]=x21*z13-x31*z12;
  
        c[0]=x42*y32-x32*y42;
        c[1]=x31*y43-x34*y13;
        c[2]=x24*y14-x14*y24;
        c[3]=x13*y21-x12*y31;
        double Si, area;
        for(int j=0;j<4;j++)
        {
               
            Si = sqrt(a[j]*a[j]+b[j]*b[j]+c[j]*c[j]);
            area=0.5*Si;
            elementInfo[i].area[j]=area;
            elementInfo[i].nx[j]= -a[j]/Si;
            elementInfo[i].ny[j]= -b[j]/Si;
            elementInfo[i].nz[j]= -c[j]/Si;

         }

   return;
  }



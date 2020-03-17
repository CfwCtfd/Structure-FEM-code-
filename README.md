# Structure-FEM-code-


  Purpose:

    This is an FEM solver on Tetrahedral Mesh with 10 nodes

  Discussion:

    Solves the Elasticity equation for bending problem which arises in aero-elasticity
    although other load condition can be easily implemented. This code works with the 
    pre-processor MeshPreProcessor.cpp which in turn uses the c++ code of John Burkardt
    tet_mesh_l2q.cpp to generate the tet10 elements. The basic mesh is genarated using 
    pointwise and exported for SU2 mesh format. The mesh is further processed using 
    to generate three input files for the current program namely:
    1. element.dat which has the element node information as well as info on 
                     the faces, normals, areas, volumes

    2. nodes.dat contains the nodes coordinate

    3. BoundaryOrdering.dat which has boundary node information.

    References: The Quadratic Tetrahedron:lecture notes by Carlos Fellipa.


   Author: Kaushik K N 

 

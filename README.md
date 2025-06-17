# CrackSystems
2D linear elasticity plane strain FEM codes for the computation of stresses in systems of en echelon and randomly distributed cracks

README

In this repository, you may find the codes to calculate stress tensor components in 2D section of a body with different crack systems undergoing linear elastic deformation.
The codes are written in the MATLAB environment, and they are based on the plane strain approximation.
The main files are those named Mesh_4 and Mesh_5. 
Mesh_4 is for the systems of en echelon cracks, Mesh_5 for the randomly distributed cracks. In the main file, all the geometric properties can be edited, the cracks length, aspect ratio, separation, orientation, the Randomness Degree, the size oh the whole domain, the orientation of the cracks' alignment, the boundary conditions (Neumann and Dirichlet) and the material properties. Files MeshM4 and MeshM5 are the correspondent function files which generate the mesh. Note that the meshes are produced through triangle mesh generator, which is an open source unstructured mesh generator that can be downloaded. For reference, see Shewchuk (1996) and Shewchuk (2002). The file k_el_hoop is the function responsible for the computation of the elementary stiffnes matrices. Stress_computation calculates all the stress components, element by element. Nodebynodetri6 evaluates the stress components node by node (for 6-node triangular elements).
Finally, the files Plots_Mesh4 and Plots_Mesh5redblue can be use to plot the stress components.


References:

Shewchuk, J. R. (1996). Triangle: Engineering a 2D quality mesh generator and Delaunay triangulator. In Workshop on Applied Computational Geometry (pp. 203-222). Berlin, Heidelberg: Springer Berlin Heidelberg.
Shewchuk, J. R. (2002). Delaunay refinement algorithms for triangular mesh generation. Computational geometry, 22(1-3), 21-74.



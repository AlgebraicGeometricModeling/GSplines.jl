var documenterSearchIndex = {"docs":
[{"location":"iga/#IsoGeometric-Analysis-and-Point-Cloud-Fitting","page":"IsoGeometric Analysis and Point Cloud Fitting","title":"IsoGeometric Analysis and Point Cloud Fitting","text":"","category":"section"},{"location":"iga/","page":"IsoGeometric Analysis and Point Cloud Fitting","title":"IsoGeometric Analysis and Point Cloud Fitting","text":"Pages = [\"iga.md\"]","category":"page"},{"location":"iga/","page":"IsoGeometric Analysis and Point Cloud Fitting","title":"IsoGeometric Analysis and Point Cloud Fitting","text":"g1dimension\ng1hom_matrix","category":"page"},{"location":"iga/#GSplines.g1dimension","page":"IsoGeometric Analysis and Point Cloud Fitting","title":"GSplines.g1dimension","text":"g1dimension(m::{HMesh, Mesh})\n\nCount the dimension of the space of biquintic g1 splines on the mesh, with symmetric quadratic glueing data, using a combinatorial formula depending on the type of edges and vertices in the mesh.\n\n\n\n\n\ng1dimension(m::{HMesh, Mesh}, kn::Vector)\n\nComputes the dimension of the space of G1 splines associated to the knot sequence kn on the mesh m, with symmetric quadratic glueing data.\n\n\n\n\n\n","category":"function"},{"location":"iga/#GSplines.g1hom_matrix","page":"IsoGeometric Analysis and Point Cloud Fitting","title":"GSplines.g1hom_matrix","text":"Compute the matrix defining the G1 splines on the mesh hm with the knot distribution kn using quadratic glueing data.\n\nIt assumes that boundary vertices are regular.\n\n\n\n\n\n","category":"function"},{"location":"iga/","page":"IsoGeometric Analysis and Point Cloud Fitting","title":"IsoGeometric Analysis and Point Cloud Fitting","text":"g1basis","category":"page"},{"location":"iga/#GSplines.g1basis","page":"IsoGeometric Analysis and Point Cloud Fitting","title":"GSplines.g1basis","text":"This function takes in input a quad mesh and returns a sparse matrix containing the coefficients defining a set of G1 biquintic basis functions on the input mesh. The construction holds for both planar and nonplanar quad meshes. The input mesh is assumed to be with isolated EVs and with no boundary EVs.\n\nncols of the sparse matrix gives the dimension of the spline space,\nnrows is the total number of control points in the mesh i.e. nfaces*36\n\nExample\n\nusing GSplines\nm = offdata(\"triangle_planar.off\")\nbasis = g1basis(m)\n\n\n\n\n\n","category":"function"},{"location":"iga/","page":"IsoGeometric Analysis and Point Cloud Fitting","title":"IsoGeometric Analysis and Point Cloud Fitting","text":"ToGismo","category":"page"},{"location":"iga/#GSplines.ToGismo","page":"IsoGeometric Analysis and Point Cloud Fitting","title":"GSplines.ToGismo","text":"This function produces the .xml file containing both the spline geometry and associated biquintic basis functions (stored as sparse matrix) that can be loaded in G+Smo scripts to be used in IsoGeometric Analysis simulations or in point cloud fitting problems. \n\nIt may also creates .xml files containing only the geometry or the set of basis functions.\n\nExample\n\nusing GSplines\nm = offdata(\"triangle_planar.off\")\ns1 = g1surface(m)\nbasis = g1basis(m)\nToGismo(s1,basis,s1.knots,'filename')\n\n\n\n\n\n","category":"function"},{"location":"package/#About-MomentTools.jl","page":"About the package","title":"About MomentTools.jl","text":"","category":"section"},{"location":"package/#Installation","page":"About the package","title":"Installation","text":"","category":"section"},{"location":"package/","page":"About the package","title":"About the package","text":"The package can be installed from Julia as follows:","category":"page"},{"location":"package/","page":"About the package","title":"About the package","text":"] add https://github.com/AlgebraicGeometricModeling/GSplines.jl\n","category":"page"},{"location":"package/","page":"About the package","title":"About the package","text":"A version of Julia >= 1.3 should be used.","category":"page"},{"location":"package/","page":"About the package","title":"About the package","text":"It depends on the following packages, which are installed automatically:","category":"page"},{"location":"package/","page":"About the package","title":"About the package","text":"SemiAlgebraicTypes","category":"page"},{"location":"package/#Development","page":"About the package","title":"Development","text":"","category":"section"},{"location":"package/","page":"About the package","title":"About the package","text":"The git project of the package is     https://github.com/AlgebraicGeometricModeling/GSplines.jl.","category":"page"},{"location":"package/","page":"About the package","title":"About the package","text":"The main developpers are (so far)","category":"page"},{"location":"package/","page":"About the package","title":"About the package","text":"Michelangelo Marsala\nBernard Mourrain","category":"page"},{"location":"converters/#Readers,-writers,-converters","page":"Readers, writers, converters","title":"Readers, writers, converters","text":"","category":"section"},{"location":"converters/","page":"Readers, writers, converters","title":"Readers, writers, converters","text":"Pages = [\"converters.md\"]","category":"page"},{"location":"converters/","page":"Readers, writers, converters","title":"Readers, writers, converters","text":"set_boundary!","category":"page"},{"location":"converters/#GSplines.set_boundary!","page":"Readers, writers, converters","title":"GSplines.set_boundary!","text":"This function modifies the half edge data structure of a mesh re-defining some edges to be boundaries.      It is done giving in input a list of indices identifying the edges (in the half edge structure) to be set as boundaries.\n\nExample\n\n  hm = hmesh(offdata(\"cube.off\"))\n  set_boundary!(hm, [1,2])\n\n\n\n\n\n","category":"function"},{"location":"converters/","page":"Readers, writers, converters","title":"Readers, writers, converters","text":"Here is an example of how to set an edge as boundary and visualize it with 'Axl':","category":"page"},{"location":"converters/","page":"Readers, writers, converters","title":"Readers, writers, converters","text":"using GSplines, SemiAlgebraicTypes, Axl\n  \n  hm = hmesh(offdata(\"cube.off\")) # read the mesh and convert it in half edge structure\n  set_boundary!(hm,[1 2])        # define the edge with indicized vertices [1 15] as boundary\n  s1 = g1surface(hm)              # create the G1 surface on the modifies half mesh\n  @axlview s1                     # view the surface using Axl","category":"page"},{"location":"converters/","page":"Readers, writers, converters","title":"Readers, writers, converters","text":"(Image: scube)","category":"page"},{"location":"converters/","page":"Readers, writers, converters","title":"Readers, writers, converters","text":"bspline","category":"page"},{"location":"converters/#GSplines.bspline","page":"Readers, writers, converters","title":"GSplines.bspline","text":"Extract the bspline representation of face f\n\n\n\n\n\nExtract the bspline representations of the faces as a vector of bspline  functions\n\n\n\n\n\n","category":"function"},{"location":"converters/","page":"Readers, writers, converters","title":"Readers, writers, converters","text":"offread\noffdata","category":"page"},{"location":"converters/#GSplines.offread","page":"Readers, writers, converters","title":"GSplines.offread","text":"Read an off file and ouput a mesh.\n\noffread(\"file.off\")\n\n\n\n\n\n","category":"function"},{"location":"converters/#GSplines.offdata","page":"Readers, writers, converters","title":"GSplines.offdata","text":"Read the off data in the 'file' of the data repository.\n\nm = offdata(\"cube.off\")\n\n\n\n\n\n","category":"function"},{"location":"converters/","page":"Readers, writers, converters","title":"Readers, writers, converters","text":"axldata","category":"page"},{"location":"converters/#GSplines.axldata","page":"Readers, writers, converters","title":"GSplines.axldata","text":"Read the axl data in the 'file' of the data repository.\n\nm = axldata(\"y1m1.axl\")\n\n\n\n\n\n","category":"function"},{"location":"gismo/#Project-using-GSmo","page":"Using G+Smo","title":"Project using G+Smo","text":"","category":"section"},{"location":"gismo/#How-to-compile-GSmo","page":"Using G+Smo","title":"How to compile G+Smo","text":"","category":"section"},{"location":"gismo/","page":"Using G+Smo","title":"Using G+Smo","text":"We can point our new project directly to the build folder where G+Smo is compiled:","category":"page"},{"location":"gismo/","page":"Using G+Smo","title":"Using G+Smo","text":"cd /path/to/myproject\ncmake . -Dgismo_DIR=/path/to/build","category":"page"},{"location":"gismo/","page":"Using G+Smo","title":"Using G+Smo","text":"Alternatively, if we would like to point to our install location the path would change:","category":"page"},{"location":"gismo/","page":"Using G+Smo","title":"Using G+Smo","text":"cmake . -Dgismo_DIR=/path/to/install","category":"page"},{"location":"gismo/","page":"Using G+Smo","title":"Using G+Smo","text":"This will create a \"Makefile\" (therefore will overwrite any existing makefiles).  For this solution the installation of G+Smo is not mandatory.","category":"page"},{"location":"gismo/","page":"Using G+Smo","title":"Using G+Smo","text":"It is important to see in the output of CMake the message:","category":"page"},{"location":"gismo/","page":"Using G+Smo","title":"Using G+Smo","text":"\"G+Smo shared library found in /some/path\"","category":"page"},{"location":"gismo/","page":"Using G+Smo","title":"Using G+Smo","text":"This reveals the actual dynamic or static library file that will be used for linking. If the path is not the one instructed in gsimoDIR, then we might need to remove the CMakeCache.txt file and try again. Another complication might be that CMake usually \"remembers\" a specific build location and uses that one, eg. when we do not explicitly set gismoDIR or the path we set is not correct. In that case we should delete the cache file (CMakeCache.txt) and also remove the following folder:","category":"page"},{"location":"gismo/","page":"Using G+Smo","title":"Using G+Smo","text":"rm -rf ~/.cmake/packages/gismo","category":"page"},{"location":"gismo/","page":"Using G+Smo","title":"Using G+Smo","text":"and restart the process.","category":"page"},{"location":"gismo/","page":"Using G+Smo","title":"Using G+Smo","text":"If the configuration is successfull, then typing \"make\" should create the binary file \"main\", which can be executed as:","category":"page"},{"location":"gismo/","page":"Using G+Smo","title":"Using G+Smo","text":"make\n./main","category":"page"},{"location":"gismo/","page":"Using G+Smo","title":"Using G+Smo","text":"Hello G+Smo.","category":"page"},{"location":"gismo/","page":"Using G+Smo","title":"Using G+Smo","text":"Note A: CMake can produce project files for Visual Studio, Ninja, CodeBlocks, and many other build systems, see here.","category":"page"},{"location":"gismo/","page":"Using G+Smo","title":"Using G+Smo","text":"Note B: On Windows (but also linux), we can use the graphical tool cmake-gui to configure and compile. On linux the interactive tool ccmake is also an option.","category":"page"},{"location":"gismo/#How-to-use-GSmo-in-your-application","page":"Using G+Smo","title":"How to use G+Smo in your application","text":"","category":"section"},{"location":"gismo/","page":"Using G+Smo","title":"Using G+Smo","text":"Assuming that you have the source tree at /path/to/gismo, here is how to use G+Smo in your application:","category":"page"},{"location":"gismo/","page":"Using G+Smo","title":"Using G+Smo","text":"Choose a folder for building (eg. folder /path/to/build for demonstation):","category":"page"},{"location":"gismo/","page":"Using G+Smo","title":"Using G+Smo","text":"mkdir /path/to/build\ncd /path/to/build","category":"page"},{"location":"gismo/","page":"Using G+Smo","title":"Using G+Smo","text":"Configure, build and install on a predefined location (eg. `./path/to/install for demonstration):","category":"page"},{"location":"gismo/","page":"Using G+Smo","title":"Using G+Smo","text":"The configure step:","category":"page"},{"location":"gismo/","page":"Using G+Smo","title":"Using G+Smo","text":"cmake /path/to/gismo -DCMAKE_INSTALL_PREFIX=/path/to/install -DGISMO_EXAMPLES=OFF","category":"page"},{"location":"gismo/","page":"Using G+Smo","title":"Using G+Smo","text":"The build step:","category":"page"},{"location":"gismo/","page":"Using G+Smo","title":"Using G+Smo","text":"make -j2","category":"page"},{"location":"gismo/","page":"Using G+Smo","title":"Using G+Smo","text":"Note that here we assume that Makefiles were used (autotools).  We could explicitly choose a builder during configuration or CMake will resort to an available builder in the system.  Another popular builder is ninja, in which case we would always issue \"ninja\" instead of \"make\".","category":"page"},{"location":"gismo/","page":"Using G+Smo","title":"Using G+Smo","text":"The install step (optional), requires that the build step above succeeded:","category":"page"},{"location":"gismo/","page":"Using G+Smo","title":"Using G+Smo","text":"make -j2 install","category":"page"},{"location":"gismo/","page":"Using G+Smo","title":"Using G+Smo","text":"The installation will create the folder /path/to/install, containing:","category":"page"},{"location":"gismo/","page":"Using G+Smo","title":"Using G+Smo","text":"./include/gismo   Include path needed by the library (eg. gismo.h)\n./lib             Shared library (eg. libgismo.so)\n./bin             Compiled examples (only if GISMO_EXAMPLES=ON)","category":"page"},{"location":"gismo/","page":"Using G+Smo","title":"Using G+Smo","text":"During development that might also involve frequent re-building of G+Smo the install step would have to be executed after each build. Therefore one could also skip this step and use the build folder to deploy the library to an external project.","category":"page"},{"location":"gismo/","page":"Using G+Smo","title":"Using G+Smo","text":"There are several solutions for using the library in your project:","category":"page"},{"location":"#GSplines","page":"Home","title":"GSplines","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"The package provide tools for the construction and use of Geometrically Smooth (G1) Splines .","category":"page"},{"location":"","page":"Home","title":"Home","text":"Here is an example of construction of a G1 surface, visualized with 'Axl':","category":"page"},{"location":"","page":"Home","title":"Home","text":"using GSplines, SemiAlgebraicTypes, Axl\n\nm = axldata(\"y1m1.axl\")[1]  # read mesh from data file \"y1m1.axl\"\nm[:color] = Axl.red         # set its color to red\n@axlview m                  # view the mesh using Axl","category":"page"},{"location":"","page":"Home","title":"Home","text":"(Image: y1m1)","category":"page"},{"location":"","page":"Home","title":"Home","text":"s = g1surface(m)            # compute a G1 smooth surface from the mesh m \n@axlview s                  # view the G1 spline surface using Axl","category":"page"},{"location":"","page":"Home","title":"Home","text":"(Image: y1g1)","category":"page"},{"location":"","page":"Home","title":"Home","text":"Here is another example with an off file.","category":"page"},{"location":"","page":"Home","title":"Home","text":"m1 = offdata(\"cube.off\")    # read mesh in off format from data file \"cube.off\" \nm1[:color]= Axl.blue        # set its color to blue\n@axlview m1                 # view the mesh using Axl","category":"page"},{"location":"","page":"Home","title":"Home","text":"(Image: m1cube)","category":"page"},{"location":"","page":"Home","title":"Home","text":"s1 = g1surface(m1)          # compute a G1 smooth surface from the mesh m1 \n@axlview s1                 # view the G1 spline using Axl ","category":"page"},{"location":"","page":"Home","title":"Home","text":"(Image: g1cube)","category":"page"},{"location":"","page":"Home","title":"Home","text":"Notice that the mesh has been subdivided using one step of Catmull-Clark subdivision, so that no extraordinary vertices (of valence neq 4) are connected by an edge.","category":"page"},{"location":"geometric_model/#Geometric-models","page":"Geometric models","title":"Geometric models","text":"","category":"section"},{"location":"geometric_model/","page":"Geometric models","title":"Geometric models","text":"Pages = [\"geometric_models.md\"]","category":"page"},{"location":"geometric_model/","page":"Geometric models","title":"Geometric models","text":"g0surface","category":"page"},{"location":"geometric_model/#GSplines.g0surface","page":"Geometric models","title":"GSplines.g0surface","text":"g0surface(hm::{HMesh, Mesh}, d=3)\n\nThis function takes in input a mesh (or an half egde data structure of a mesh) and returns a (bicubic) g0 multipatch surface obtained using the Approximate Catmull-Clark scheme (ACC3).     If the input degree is greater than 3, the resulting surface is degree elevated to the desired degree.\n\nExample\n\nusing G1Splines\nm = offdata(\"cube.off\")\ng0 = g0surface(m)\n\n\n\n\n\n","category":"function"},{"location":"geometric_model/","page":"Geometric models","title":"Geometric models","text":"g1surface","category":"page"},{"location":"geometric_model/#GSplines.g1surface","page":"Geometric models","title":"GSplines.g1surface","text":"g1surface(hm::HMesh,S::String=\"CS-S\")\n\nThis function takes in input a mesh in the half egde data structure and a string with the G1 solving strategy and will return an array composed by the patches constituing the surface.\n\nBy default, if two EVs are connected, the mesh is split. If the option check_ev = false, this does not happen, but the construction may be wrong.\n\nThe input string contains the solving strategy for the construction of the G1 surface to be selected from the     following four: \"CS-S\",\"CS-AS\",\"NCS-S\",\"NCS-AS\". The default strategy is \"CS-S\".\n\nExample\n\nusing G1Splines\nm = offdata(\"cube.off\")\ng1 = g1surface(m)\n\n\n\n\n\n","category":"function"}]
}

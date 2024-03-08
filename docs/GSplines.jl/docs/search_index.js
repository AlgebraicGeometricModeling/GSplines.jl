var documenterSearchIndex = {"docs":
[{"location":"iga/#IsoGeometric-Analysis","page":"IsoGeometric Analysis","title":"IsoGeometric Analysis","text":"","category":"section"},{"location":"iga/","page":"IsoGeometric Analysis","title":"IsoGeometric Analysis","text":"Pages = [\"iga.md\"]","category":"page"},{"location":"iga/","page":"IsoGeometric Analysis","title":"IsoGeometric Analysis","text":"g1basis","category":"page"},{"location":"iga/#GSplines.g1basis","page":"IsoGeometric Analysis","title":"GSplines.g1basis","text":"This function takes in input a quad mesh and returns a sparse matrix containing the coefficients defining a set of G1 biquintic basis functions on the input mesh.\n\nncols of the sparse matrix gives the dimension of the spline space,\nnrows is the total number of control points in the mesh i.e. nfaces*36\n\n\n\n\n\n","category":"function"},{"location":"package/#About-MomentTools.jl","page":"About the package","title":"About MomentTools.jl","text":"","category":"section"},{"location":"package/#Installation","page":"About the package","title":"Installation","text":"","category":"section"},{"location":"package/","page":"About the package","title":"About the package","text":"The package can be installed from Julia as follows:","category":"page"},{"location":"package/","page":"About the package","title":"About the package","text":"] add https://github.com/AlgebraicGeometricModeling/GSplines.jl\n","category":"page"},{"location":"package/","page":"About the package","title":"About the package","text":"A version of Julia >= 1.3 should be used.","category":"page"},{"location":"package/","page":"About the package","title":"About the package","text":"It depends on the following packages, which are installed automatically:","category":"page"},{"location":"package/","page":"About the package","title":"About the package","text":"SemiAlgebraicTypes","category":"page"},{"location":"package/#Development","page":"About the package","title":"Development","text":"","category":"section"},{"location":"package/","page":"About the package","title":"About the package","text":"The git project of the package is     https://github.com/AlgebraicGeometricModeling/GSplines.jl.","category":"page"},{"location":"package/","page":"About the package","title":"About the package","text":"The main developpers are (so far)","category":"page"},{"location":"package/","page":"About the package","title":"About the package","text":"Michelangelo Marsala\nBernard Mourrain","category":"page"},{"location":"converters/#Readers,-writers,-converters","page":"Readers, writers, converters","title":"Readers, writers, converters","text":"","category":"section"},{"location":"converters/","page":"Readers, writers, converters","title":"Readers, writers, converters","text":"Pages = [\"converters.md\"]","category":"page"},{"location":"converters/","page":"Readers, writers, converters","title":"Readers, writers, converters","text":"bspline","category":"page"},{"location":"converters/#GSplines.bspline","page":"Readers, writers, converters","title":"GSplines.bspline","text":"Extract the bspline representation of face f\n\n\n\n\n\nExtract the bspline representations of the faces as a vector of bspline  functions\n\n\n\n\n\n","category":"function"},{"location":"converters/","page":"Readers, writers, converters","title":"Readers, writers, converters","text":"offread\noffdata","category":"page"},{"location":"converters/#GSplines.offread","page":"Readers, writers, converters","title":"GSplines.offread","text":"Read an off file and ouput a mesh.\n\noffread(\"file.off\")\n\n\n\n\n\n","category":"function"},{"location":"converters/#GSplines.offdata","page":"Readers, writers, converters","title":"GSplines.offdata","text":"Read the off data in the 'file' of the data repository.\n\nm = offdata(\"cube.off\")\n\n\n\n\n\n","category":"function"},{"location":"converters/","page":"Readers, writers, converters","title":"Readers, writers, converters","text":"axldata","category":"page"},{"location":"converters/#GSplines.axldata","page":"Readers, writers, converters","title":"GSplines.axldata","text":"Read the off data in the 'file' of the data repository.\n\nm = axldata(\"y1m1.axl\")\n\n\n\n\n\n","category":"function"},{"location":"gismo/#Project-using-GSmo","page":"Using G+Smo","title":"Project using G+Smo","text":"","category":"section"},{"location":"gismo/#How-to-compile-GSmo","page":"Using G+Smo","title":"How to compile G+Smo","text":"","category":"section"},{"location":"gismo/","page":"Using G+Smo","title":"Using G+Smo","text":"We can point our new project directly to the build folder where G+Smo is compiled:","category":"page"},{"location":"gismo/","page":"Using G+Smo","title":"Using G+Smo","text":"cd /path/to/myproject\ncmake . -Dgismo_DIR=/path/to/build","category":"page"},{"location":"gismo/","page":"Using G+Smo","title":"Using G+Smo","text":"Alternatively, if we would like to point to our install location the path would change:","category":"page"},{"location":"gismo/","page":"Using G+Smo","title":"Using G+Smo","text":"cmake . -Dgismo_DIR=/path/to/install","category":"page"},{"location":"gismo/","page":"Using G+Smo","title":"Using G+Smo","text":"This will create a \"Makefile\" (therefore will overwrite any existing makefiles).  For this solution the installation of G+Smo is not mandatory.","category":"page"},{"location":"gismo/","page":"Using G+Smo","title":"Using G+Smo","text":"It is important to see in the output of CMake the message:","category":"page"},{"location":"gismo/","page":"Using G+Smo","title":"Using G+Smo","text":"\"G+Smo shared library found in /some/path\"","category":"page"},{"location":"gismo/","page":"Using G+Smo","title":"Using G+Smo","text":"This reveals the actual dynamic or static library file that will be used for linking. If the path is not the one instructed in gsimoDIR, then we might need to remove the CMakeCache.txt file and try again. Another complication might be that CMake usually \"remembers\" a specific build location and uses that one, eg. when we do not explicitly set gismoDIR or the path we set is not correct. In that case we should delete the cache file (CMakeCache.txt) and also remove the following folder:","category":"page"},{"location":"gismo/","page":"Using G+Smo","title":"Using G+Smo","text":"rm -rf ~/.cmake/packages/gismo","category":"page"},{"location":"gismo/","page":"Using G+Smo","title":"Using G+Smo","text":"and restart the process.","category":"page"},{"location":"gismo/","page":"Using G+Smo","title":"Using G+Smo","text":"If the configuration is successfull, then typing \"make\" should create the binary file \"main\", which can be executed as:","category":"page"},{"location":"gismo/","page":"Using G+Smo","title":"Using G+Smo","text":"make\n./main","category":"page"},{"location":"gismo/","page":"Using G+Smo","title":"Using G+Smo","text":"Hello G+Smo.","category":"page"},{"location":"gismo/","page":"Using G+Smo","title":"Using G+Smo","text":"Note A: CMake can produce project files for Visual Studio, Ninja, CodeBlocks, and many other build systems, see here.","category":"page"},{"location":"gismo/","page":"Using G+Smo","title":"Using G+Smo","text":"Note B: On Windows (but also linux), we can use the graphical tool cmake-gui to configure and compile. On linux the interactive tool ccmake is also an option.","category":"page"},{"location":"gismo/#How-to-use-GSmo-in-your-application","page":"Using G+Smo","title":"How to use G+Smo in your application","text":"","category":"section"},{"location":"gismo/","page":"Using G+Smo","title":"Using G+Smo","text":"Assuming that you have the source tree at /path/to/gismo, here is how to use G+Smo in your application:","category":"page"},{"location":"gismo/","page":"Using G+Smo","title":"Using G+Smo","text":"Choose a folder for building (eg. folder /path/to/build for demonstation):","category":"page"},{"location":"gismo/","page":"Using G+Smo","title":"Using G+Smo","text":"mkdir /path/to/build\ncd /path/to/build","category":"page"},{"location":"gismo/","page":"Using G+Smo","title":"Using G+Smo","text":"Configure, build and install on a predefined location (eg. `./path/to/install for demonstration):","category":"page"},{"location":"gismo/","page":"Using G+Smo","title":"Using G+Smo","text":"The configure step:","category":"page"},{"location":"gismo/","page":"Using G+Smo","title":"Using G+Smo","text":"cmake /path/to/gismo -DCMAKE_INSTALL_PREFIX=/path/to/install -DGISMO_EXAMPLES=OFF","category":"page"},{"location":"gismo/","page":"Using G+Smo","title":"Using G+Smo","text":"The build step:","category":"page"},{"location":"gismo/","page":"Using G+Smo","title":"Using G+Smo","text":"make -j2","category":"page"},{"location":"gismo/","page":"Using G+Smo","title":"Using G+Smo","text":"Note that here we assume that Makefiles were used (autotools).  We could explicitly choose a builder during configuration or CMake will resort to an available builder in the system.  Another popular builder is ninja, in which case we would always issue \"ninja\" instead of \"make\".","category":"page"},{"location":"gismo/","page":"Using G+Smo","title":"Using G+Smo","text":"The install step (optional), requires that the build step above succeeded:","category":"page"},{"location":"gismo/","page":"Using G+Smo","title":"Using G+Smo","text":"make -j2 install","category":"page"},{"location":"gismo/","page":"Using G+Smo","title":"Using G+Smo","text":"The installation will create the folder /path/to/install, containing:","category":"page"},{"location":"gismo/","page":"Using G+Smo","title":"Using G+Smo","text":"./include/gismo   Include path needed by the library (eg. gismo.h)\n./lib             Shared library (eg. libgismo.so)\n./bin             Compiled examples (only if GISMO_EXAMPLES=ON)","category":"page"},{"location":"gismo/","page":"Using G+Smo","title":"Using G+Smo","text":"During development that might also involve frequent re-building of G+Smo the install step would have to be executed after each build. Therefore one could also skip this step and use the build folder to deploy the library to an external project.","category":"page"},{"location":"gismo/","page":"Using G+Smo","title":"Using G+Smo","text":"There are several solutions for using the library in your project:","category":"page"},{"location":"#GSplines","page":"Home","title":"GSplines","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"The package provide tools for the construction and use of Geometrically Smooth (G1) Splines .","category":"page"},{"location":"","page":"Home","title":"Home","text":"Here is an example of construction of a G1 surface, visualized with 'Axl':","category":"page"},{"location":"","page":"Home","title":"Home","text":"using GSplines, SemiAlgebraicTypes, Axl\n\nm = axldata(\"y1m1.axl\")[1]  #read mesh from data file \"y1m1.axl\"\nm[:color] = Axl.red","category":"page"},{"location":"","page":"Home","title":"Home","text":"(Image: y1m1)","category":"page"},{"location":"","page":"Home","title":"Home","text":"s = g1surface(m)\n@axlview s","category":"page"},{"location":"","page":"Home","title":"Home","text":"(Image: y1g1)","category":"page"},{"location":"geometric_model/#Geometric-models","page":"Geometric models","title":"Geometric models","text":"","category":"section"},{"location":"geometric_model/","page":"Geometric models","title":"Geometric models","text":"Pages = [\"geometric_models.md\"]","category":"page"},{"location":"geometric_model/","page":"Geometric models","title":"Geometric models","text":"g0surface","category":"page"},{"location":"geometric_model/#GSplines.g0surface","page":"Geometric models","title":"GSplines.g0surface","text":"g0surface(hm::HMesh,d=3), g0surface(m::Mesh,d=3)\n\nThis function takes in input a mesh (or an half egde data structure of a mesh) and returns a g0 multipatch surface obtained using the Approximate Catmull-Clark scheme (ACC3)\n\nExample\n\nusing G1Splines\nm = offdata(\"cube.off\")\ng0 = g0surface(m)\n\n\n\n\n\n","category":"function"},{"location":"geometric_model/","page":"Geometric models","title":"Geometric models","text":"g1surface","category":"page"},{"location":"geometric_model/#GSplines.g1surface","page":"Geometric models","title":"GSplines.g1surface","text":"g1surface(hm::HMesh,S::String=\"CS-S\")\n\nThis function takes in input a mesh in the half egde data structure and a string with the G1 solving strategy and will return an array composed by the patches constituing the surface.\n\nBy default, if two EVs are connected, the mesh is split. If the option check_ev = false, this does not happen, but the construction may be wrong.\n\nExample\n\nusing G1Splines\nm = offdata(\"cube.off\")\ng1 = g1surface(m)\n\n\n\n\n\n","category":"function"},{"location":"geometric_model/","page":"Geometric models","title":"Geometric models","text":"GSplines.g1matrix\nGSplines.G1matrix\nGSplines.g1matrix_edge\nGSplines.fake_edges","category":"page"},{"location":"geometric_model/#GSplines.g1matrix","page":"Geometric models","title":"GSplines.g1matrix","text":"Compute the matrix defining the G1 splines on the mesh hm with the knot distribution kn using quadratic glueing data.\n\nIt assumes that boundary vertices are regular.\n\n\n\n\n\n","category":"function"},{"location":"geometric_model/#GSplines.G1matrix","page":"Geometric models","title":"GSplines.G1matrix","text":"This function creates the smoothing matrix to create a G1 smooth patch around extraordinary vertices of any valence\n\nN–> Valence of the EV\nS–> Solving strategy NCS-S (non circulant system+symm), NCS-AS(non circulant system+asymm), CS-AS(circulant system+asymm), CS-S(circulant system+symm)\n\nT–>Type of mesh face The output matrix B has smoothing masks ordered as follow\n\nControl points labelling for inner EVs       Biquintic patch\n\n    28 31 34 | 21 20 19\n    29 32 35 | 24 23 22\n    30 33 36 | 25 26 27\n    -------------------       B=[1|2|3|...|34|35|36] columns OF B\n    7  8  9  | 18 15 12\n    4  5  6  | 17 14 11\n    1  2  3  | 16 13 10\n\nFor regular vertices N=4 it will return ACC3 smoothing masks\n\n\n\n\n\n","category":"function"},{"location":"geometric_model/#GSplines.g1matrix_edge","page":"Geometric models","title":"GSplines.g1matrix_edge","text":"Compute the matrix defining the edge spline basis, in a compressed form.\n\n\n\n\n\n","category":"function"},{"location":"geometric_model/#GSplines.fake_edges","page":"Geometric models","title":"GSplines.fake_edges","text":"This function fakes the half edge datas tructure of a mesh re-defining some edges to be boundaries from a given list of edges in input\n\n\n\n\n\n","category":"function"}]
}

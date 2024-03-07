
# Project using G+Smo

## How to compile G+Smo

We can point our new project directly to the build folder where G+Smo is compiled:

```
cd /path/to/myproject
cmake . -Dgismo_DIR=/path/to/build
```

Alternatively, if we would like to point to our install location the
path would change:

```
cmake . -Dgismo_DIR=/path/to/install
```

This will create a "Makefile" (therefore will overwrite any existing
makefiles).  For this solution the installation of G+Smo is not
mandatory.

It is important to see in the output of CMake the message:

"G+Smo shared library found in /some/path"

This reveals the actual dynamic or static library file that will be
used for linking. If the path is not the one instructed in gsimo_DIR,
then we might need to remove the CMakeCache.txt file and try again.
Another complication might be that CMake usually "remembers" a
specific build location and uses that one, eg. when we do not
explicitly set gismo_DIR or the path we set is not correct. In that
case we should delete the cache file (CMakeCache.txt) and also remove
the following folder:

```
rm -rf ~/.cmake/packages/gismo
```

and restart the process.

If the configuration is successfull, then typing "make" should create
the binary file "main", which can be executed as:

```
make
./main
```

```Hello G+Smo.```

**Note A:** CMake can produce project files for Visual Studio, Ninja,
CodeBlocks, and many other build systems,
see [here](https://cmake.org/cmake/help/v3.0/manual/cmake-generators.7.html).

**Note B:** On Windows (but also linux), we can use the graphical tool
cmake-gui to configure and compile. On linux the interactive tool
ccmake is also an option.


## How to use G+Smo in your application

Assuming that you have the source tree at /path/to/gismo, here is how
to use G+Smo in your application:


* Choose a folder for building (eg. folder /path/to/build for demonstation):

```
mkdir /path/to/build
cd /path/to/build
```

* Configure, build and install on a predefined location
  (eg. `./path/to/install for demonstration):

The configure step:
```
cmake /path/to/gismo -DCMAKE_INSTALL_PREFIX=/path/to/install -DGISMO_EXAMPLES=OFF
```

The build step:

```
make -j2
```

Note that here we assume that Makefiles were used (autotools).  We
could explicitly choose a builder during configuration or CMake will
resort to an available builder in the system.  Another popular builder
is ninja, in which case we would always issue "ninja" instead of "make".

The install step (optional), requires that the build step above succeeded:

```
make -j2 install
```

The installation will create the folder /path/to/install, containing:

- ./include/gismo   Include path needed by the library (eg. gismo.h)
- ./lib             Shared library (eg. libgismo.so)
- ./bin             Compiled examples (only if GISMO_EXAMPLES=ON)

During development that might also involve frequent re-building of
G+Smo the install step would have to be executed after each build.
Therefore one could also skip this step and use the build folder to
deploy the library to an external project.


There are several solutions for using the library in your project:

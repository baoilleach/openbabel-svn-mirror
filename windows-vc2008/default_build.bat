@rem The cmake/bin directory should be in PATH
@rem The top level directory of Eigen should be in the the environment variable
@rem EIGEN2_INCLUDE_DIR. 
@rem If Eigen is not found then qeq.cpp and qtpie.cpp are not built.
@rem To build the GUI, the top level directory of wxWidgets should be in the
@rem environment variable WXWIN and -DBUILD_GUI=ON specified 
@if NOT EXIST build mkdir build
@cd build
cmake.exe -DCMAKE_INSTALL_PREFIX=..\install -G "Visual Studio 9 2008" -DLIBXML2_LIBRARIES="%CD%\..\libs\i386\libxml2.lib" -DMINIMAL_BUILD=OFF -DLIBXML2_INCLUDE_DIR=. -DZLIB_LIBRARY="%CD%\..\libs\i386\zlib1.lib" -DZLIB_INCLUDE_DIR=. -DINCHI_LIBRARY="%CD%\..\libs\i386\libinchi.lib" -DINCHI_INCLUDE_DIR=. -DXDR_LIBRARY="%CD%\..\libs\i386\xdr.lib" -DEIGEN2_INCLUDE_DIR="%EIGEN2_INCLUDE_DIR%" -DRUN_SWIG=ON -DJAVA_BINDINGS=ON -DENABLE_TESTS=OFF -DBUILD_GUI=ON -DwxWidgets_ROOT_DIR="%WXWIN%" -DwxWidgets_LIB_DIR="%WXWIN%/lib/vc_lib" -DwxWidgets_CONFIGURATION=msw %1 %2 %3 %4 %5 %6 %7 %8 ..\..
@cd ..

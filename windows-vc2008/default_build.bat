@rem -DMINIMAL_BUILD=OFF -DBUILD_GUI=ON -DENABLE_TESTS=ON

@if NOT EXIST build mkdir build
@cd build
"C:\Program Files\CMake 2.6\bin\cmake.exe" -DCMAKE_INSTALL_PREFIX=..\install -G "Visual Studio 9 2008" -DLIBXML2_LIBRARIES=..\..\windows-vc2008\libxml2.lib -DMINIMAL_BUILD=OFF -DBUILD_GUI=ON -DLIBXML2_INCLUDE_DIR=. -DZLIB_LIBRARY=..\..\windows-vc2008\zlib1.lib -DZLIB_INCLUDE_DIR=. -DENABLE_TESTS=OFF %1 %2 %3 %4 %5 %6 %7 %8 ..\..
@cd ..

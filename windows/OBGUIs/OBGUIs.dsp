# Microsoft Developer Studio Project File - Name="OBGUIs" - Package Owner=<4>
# Microsoft Developer Studio Generated Build File, Format Version 6.00
# ** DO NOT EDIT **

# TARGTYPE "Win32 (x86) Application" 0x0101

CFG=OBGUIs - Win32 Debug
!MESSAGE This is not a valid makefile. To build this project using NMAKE,
!MESSAGE use the Export Makefile command and run
!MESSAGE 
!MESSAGE NMAKE /f "OBGUIs.mak".
!MESSAGE 
!MESSAGE You can specify a configuration when running NMAKE
!MESSAGE by defining the macro CFG on the command line. For example:
!MESSAGE 
!MESSAGE NMAKE /f "OBGUIs.mak" CFG="OBGUIs - Win32 Debug"
!MESSAGE 
!MESSAGE Possible choices for configuration are:
!MESSAGE 
!MESSAGE "OBGUIs - Win32 Release" (based on "Win32 (x86) Application")
!MESSAGE "OBGUIs - Win32 Debug" (based on "Win32 (x86) Application")
!MESSAGE 

# Begin Project
# PROP AllowPerConfigDependencies 0
# PROP Scc_ProjName ""
# PROP Scc_LocalPath ""
CPP=cl.exe
MTL=midl.exe
RSC=rc.exe

!IF  "$(CFG)" == "OBGUIs - Win32 Release"

# PROP BASE Use_MFC 0
# PROP BASE Use_Debug_Libraries 0
# PROP BASE Output_Dir "Release"
# PROP BASE Intermediate_Dir "Release"
# PROP BASE Target_Dir ""
# PROP Use_MFC 1
# PROP Use_Debug_Libraries 0
# PROP Output_Dir "Release"
# PROP Intermediate_Dir "Release"
# PROP Ignore_Export_Lib 0
# PROP Target_Dir ""
# ADD BASE CPP /nologo /W3 /GX /O2 /D "WIN32" /D "NDEBUG" /D "_WINDOWS" /D "_MBCS" /YX /FD /c
# ADD CPP /nologo /MT /W3 /GR /GX /I "..\obgui" /I "..\..\src" /D "NDEBUG" /D "WIN32" /D "_WINDOWS" /D "_MBCS" /D "GUI" /FD /c
# SUBTRACT CPP /YX
# ADD BASE MTL /nologo /D "NDEBUG" /mktyplib203 /win32
# ADD MTL /nologo /D "NDEBUG" /mktyplib203 /win32
# ADD BASE RSC /l 0x809 /d "NDEBUG"
# ADD RSC /l 0x809 /d "NDEBUG"
BSC32=bscmake.exe
# ADD BASE BSC32 /nologo
# ADD BSC32 /nologo
LINK32=link.exe
# ADD BASE LINK32 kernel32.lib user32.lib gdi32.lib winspool.lib comdlg32.lib advapi32.lib shell32.lib ole32.lib oleaut32.lib uuid.lib odbc32.lib odbccp32.lib /nologo /subsystem:windows /machine:I386
# ADD LINK32 nafxcw.lib libcmt.lib /nologo /subsystem:windows /machine:I386 /nodefaultlib:"nafxcw.lib libcmt.lib" /out:"OBGUIs.exe"

!ELSEIF  "$(CFG)" == "OBGUIs - Win32 Debug"

# PROP BASE Use_MFC 0
# PROP BASE Use_Debug_Libraries 1
# PROP BASE Output_Dir "Debug"
# PROP BASE Intermediate_Dir "Debug"
# PROP BASE Target_Dir ""
# PROP Use_MFC 0
# PROP Use_Debug_Libraries 1
# PROP Output_Dir "Debug"
# PROP Intermediate_Dir "Debug"
# PROP Ignore_Export_Lib 0
# PROP Target_Dir ""
# ADD BASE CPP /nologo /W3 /Gm /GX /ZI /Od /D "WIN32" /D "_DEBUG" /D "_WINDOWS" /D "_MBCS" /YX /FD /GZ /c
# ADD CPP /nologo /MTd /W3 /Gm /GR /GX /ZI /Od /I "..\obgui" /I "..\..\src" /D "_DEBUG" /D "WIN32" /D "_WINDOWS" /D "_MBCS" /D "GUI" /FR /FD /GZ /c
# SUBTRACT CPP /YX
# ADD BASE MTL /nologo /D "_DEBUG" /mktyplib203 /win32
# ADD MTL /nologo /D "_DEBUG" /mktyplib203 /win32
# ADD BASE RSC /l 0x809 /d "_DEBUG"
# ADD RSC /l 0x809 /d "_DEBUG"
BSC32=bscmake.exe
# ADD BASE BSC32 /nologo
# ADD BSC32 /nologo
LINK32=link.exe
# ADD BASE LINK32 kernel32.lib user32.lib gdi32.lib winspool.lib comdlg32.lib advapi32.lib shell32.lib ole32.lib oleaut32.lib uuid.lib odbc32.lib odbccp32.lib /nologo /subsystem:windows /debug /machine:I386 /pdbtype:sept
# ADD LINK32 kernel32.lib user32.lib gdi32.lib winspool.lib comdlg32.lib advapi32.lib shell32.lib ole32.lib oleaut32.lib uuid.lib odbc32.lib odbccp32.lib nafxcwd.lib libcmtd.lib /nologo /subsystem:windows /debug /machine:I386 /nodefaultlib:"nafxcwd.lib libcmtd.lib" /pdbtype:sept

!ENDIF 

# Begin Target

# Name "OBGUIs - Win32 Release"
# Name "OBGUIs - Win32 Debug"
# Begin Group "Source Files"

# PROP Default_Filter "cpp;c;cxx;rc;def;r;odl;idl;hpj;bat"
# Begin Source File

SOURCE=..\..\src\atom.cpp
# End Source File
# Begin Source File

SOURCE=..\..\src\base.cpp
# End Source File
# Begin Source File

SOURCE=..\..\src\bgf.cpp
# End Source File
# Begin Source File

SOURCE=..\..\src\binary.cpp
# End Source File
# Begin Source File

SOURCE=..\..\src\bitvec.cpp
# End Source File
# Begin Source File

SOURCE=..\..\src\bond.cpp
# End Source File
# Begin Source File

SOURCE=..\..\src\box.cpp
# End Source File
# Begin Source File

SOURCE=..\..\src\c3d.cpp
# End Source File
# Begin Source File

SOURCE=..\..\src\cache.cpp
# End Source File
# Begin Source File

SOURCE=..\..\src\car.cpp
# End Source File
# Begin Source File

SOURCE=..\..\src\ccc.cpp
# End Source File
# Begin Source File

SOURCE=..\..\src\chains.cpp
# End Source File
# Begin Source File

SOURCE=..\..\src\chdrw.cpp
# End Source File
# Begin Source File

SOURCE=..\..\src\chiral.cpp
# End Source File
# Begin Source File

SOURCE=..\..\src\formats\cml.cpp
# End Source File
# Begin Source File

SOURCE=..\..\src\formats\cmlformat.cpp
# End Source File
# Begin Source File

SOURCE=..\..\src\csr.cpp
# End Source File
# Begin Source File

SOURCE=..\..\src\cssr.cpp
# End Source File
# Begin Source File

SOURCE=..\..\src\data.cpp
# End Source File
# Begin Source File

SOURCE=..\..\src\dlhandler_win32..cpp
# End Source File
# Begin Source File

SOURCE=..\..\src\dmol.cpp
# End Source File
# Begin Source File

SOURCE=..\OBGUI\DynamicOptions.cpp
# End Source File
# Begin Source File

SOURCE=..\..\src\generic.cpp
# End Source File
# Begin Source File

SOURCE=..\..\src\grid.cpp
# End Source File
# Begin Source File

SOURCE=..\..\src\hin.cpp
# End Source File
# Begin Source File

SOURCE=..\..\src\matrix.cpp
# End Source File
# Begin Source File

SOURCE=..\..\src\math\matrix3x3.cpp
# End Source File
# Begin Source File

SOURCE=..\..\src\formats\mdlformat.cpp
# End Source File
# Begin Source File

SOURCE=..\..\src\mol.cpp
# End Source File
# Begin Source File

SOURCE=..\..\src\formats\mol2format.cpp
# End Source File
# Begin Source File

SOURCE=..\..\src\molchrg.cpp
# End Source File
# Begin Source File

SOURCE=..\..\src\obconversion.cpp
# End Source File
# Begin Source File

SOURCE=..\..\src\oberror.cpp
# End Source File
# Begin Source File

SOURCE=..\OBGUI\OBGUI.cpp
# End Source File
# Begin Source File

SOURCE=..\OBGUI\OBGUI.rc
# End Source File
# Begin Source File

SOURCE=..\OBGUI\OBGUIDlg.cpp
# End Source File
# Begin Source File

SOURCE=..\..\src\obutil.cpp
# End Source File
# Begin Source File

SOURCE=..\..\src\parsmart.cpp
# End Source File
# Begin Source File

SOURCE=..\..\src\patty.cpp
# End Source File
# Begin Source File

SOURCE=..\..\src\phmodel.cpp
# End Source File
# Begin Source File

SOURCE=..\..\src\rand.cpp
# End Source File
# Begin Source File

SOURCE=..\..\src\report.cpp
# End Source File
# Begin Source File

SOURCE=..\..\src\residue.cpp
# End Source File
# Begin Source File

SOURCE=..\..\src\ring.cpp
# End Source File
# Begin Source File

SOURCE=..\..\src\rotor.cpp
# End Source File
# Begin Source File

SOURCE=..\..\src\formats\rxnformat.cpp
# End Source File
# Begin Source File

SOURCE=..\..\src\formats\smilesformat.cpp
# End Source File
# Begin Source File

SOURCE=..\OBGUI\StdAfx.cpp
# End Source File
# Begin Source File

SOURCE=..\..\src\tinker.cpp
# End Source File
# Begin Source File

SOURCE=..\..\src\tokenst.cpp
# End Source File
# Begin Source File

SOURCE=..\..\src\transform.cpp
# End Source File
# Begin Source File

SOURCE=..\..\src\formats\turbomoleformat.cpp
# End Source File
# Begin Source File

SOURCE=..\..\src\typer.cpp
# End Source File
# Begin Source File

SOURCE=..\..\src\math\vector3.cpp
# End Source File
# End Group
# Begin Group "Header Files"

# PROP Default_Filter "h;hpp;hxx;hm;inl"
# Begin Source File

SOURCE=..\..\src\aromatic.h
# End Source File
# Begin Source File

SOURCE=..\..\src\atomtyp.h
# End Source File
# Begin Source File

SOURCE=..\..\src\babelconfig.h
# End Source File
# Begin Source File

SOURCE=..\..\src\base.h
# End Source File
# Begin Source File

SOURCE=..\..\src\binary.h
# End Source File
# Begin Source File

SOURCE=..\..\src\bitvec.h
# End Source File
# Begin Source File

SOURCE=..\..\src\bondtyp.h
# End Source File
# Begin Source File

SOURCE=..\..\src\bondtyper.h
# End Source File
# Begin Source File

SOURCE=..\..\src\chains.h
# End Source File
# Begin Source File

SOURCE=..\..\src\chiral.h
# End Source File
# Begin Source File

SOURCE=..\..\src\crk.h
# End Source File
# Begin Source File

SOURCE=..\..\src\data.h
# End Source File
# Begin Source File

SOURCE=..\OBGUI\DynamicOptions.h
# End Source File
# Begin Source File

SOURCE=..\..\src\element.h
# End Source File
# Begin Source File

SOURCE=..\..\src\extable.h
# End Source File
# Begin Source File

SOURCE=..\..\src\generic.h
# End Source File
# Begin Source File

SOURCE=..\..\src\grid.h
# End Source File
# Begin Source File

SOURCE=..\..\src\isotope.h
# End Source File
# Begin Source File

SOURCE=..\..\src\matrix.h
# End Source File
# Begin Source File

SOURCE=..\..\src\mol.h
# End Source File
# Begin Source File

SOURCE=..\..\src\molchrg.h
# End Source File
# Begin Source File

SOURCE=..\..\src\obconversion.h
# End Source File
# Begin Source File

SOURCE=..\..\src\oberror.h
# End Source File
# Begin Source File

SOURCE=..\OBGUI\OBGUI.h
# End Source File
# Begin Source File

SOURCE=..\OBGUI\OBGUIDlg.h
# End Source File
# Begin Source File

SOURCE=..\..\src\obifstream.h
# End Source File
# Begin Source File

SOURCE=..\..\src\obutil.h
# End Source File
# Begin Source File

SOURCE=..\..\src\parsmart.h
# End Source File
# Begin Source File

SOURCE=..\..\src\patty.h
# End Source File
# Begin Source File

SOURCE=..\..\src\phmodel.h
# End Source File
# Begin Source File

SOURCE=..\..\src\phmodeldata.h
# End Source File
# Begin Source File

SOURCE=..\..\src\resdata.h
# End Source File
# Begin Source File

SOURCE=..\OBGUI\Resource.h
# End Source File
# Begin Source File

SOURCE=..\..\src\ring.h
# End Source File
# Begin Source File

SOURCE=..\..\src\rotor.h
# End Source File
# Begin Source File

SOURCE=..\..\src\snprintf.h
# End Source File
# Begin Source File

SOURCE=..\OBGUI\StdAfx.h
# End Source File
# Begin Source File

SOURCE=..\..\src\typer.h
# End Source File
# Begin Source File

SOURCE=..\..\src\types.h
# End Source File
# End Group
# Begin Group "Resource Files"

# PROP Default_Filter "ico;cur;bmp;dlg;rc2;rct;bin;rgs;gif;jpg;jpeg;jpe"
# End Group
# End Target
# End Project

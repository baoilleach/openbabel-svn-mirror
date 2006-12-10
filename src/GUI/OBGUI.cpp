/**********************************************************************
OBGUI.cpp -  Implementation of a cross-platform Graphical User Interface

Copyright (C) 2006 by Chris Morley

This file is part of the Open Babel project.
For more information, see <http://openbabel.sourceforge.net/>

This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation version 2 of the License.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.
***********************************************************************/
#include <openbabel/babelconfig.h>
#include <openbabel/stdwx.h>
#include <sstream>
#include <openbabel/obconversion.h>
#include <openbabel/dlhandler.h>
#include <openbabel/selformats.h>
#include <openbabel/OBGUI.h>

/*
#ifdef _DEBUG
#undef THIS_FILE
static char THIS_FILE[]=__FILE__;
#define new DEBUG_NEW
#endif
*/
using namespace OpenBabel;

// the event tables connect the wxWidgets events with the functions (event
// handlers) which process them. 
BEGIN_EVENT_TABLE(OBGUIFrame, wxFrame)
	EVT_MENU(ID_CONVERT,		 OBGUIFrame::OnConvert)
	EVT_MENU(wxID_EXIT,      OBGUIFrame::OnQuit)
	EVT_MENU(wxID_SAVE,      OBGUIFrame::OnSaveInputText)
	EVT_MENU(ID_SELFORMATS,  OBGUIFrame::OnSelectFormats)
	EVT_MENU(ID_RESTRICTFORMATS,  OBGUIFrame::OnRestrictFormats)
	EVT_MENU_RANGE(ID_SHOWCONVOPTIONS,ID_SHOWOUTOPTIONS, OBGUIFrame::OnChangeFormat)
	EVT_MENU(wxID_ABOUT, OBGUIFrame::OnAbout)
	EVT_MENU(wxID_HELP, OBGUIFrame::OnHelp)
	EVT_BUTTON(ID_INGETFILES, OBGUIFrame::OnGetInputFile)
	EVT_BUTTON(ID_OUTGETFILES, OBGUIFrame::OnGetOutputFile)
	EVT_BUTTON(ID_ININFO, OBGUIFrame::OnInFormatInfo)
	EVT_BUTTON(ID_OUTINFO, OBGUIFrame::OnOutFormatInfo)
	EVT_UPDATE_UI_RANGE(ID_OUTFILENAME,ID_OUTGETFILES, OBGUIFrame::OnOutFileNameUpdate)
	EVT_UPDATE_UI_RANGE(ID_INFILENAME, ID_INGETFILES, OBGUIFrame::OnInFileNameUpdate)
	EVT_BUTTON(ID_CONVERT, OBGUIFrame::OnConvert)
	EVT_CHECKBOX(ID_INPUTHERE,OBGUIFrame::OnChangeInputHere) 
	EVT_CHOICE(ID_INFORMAT,OBGUIFrame::OnChangeFormat)
	EVT_CHOICE(ID_OUTFORMAT,OBGUIFrame::OnChangeFormat)
	EVT_MOUSEWHEEL(OBGUIFrame::OnMouseWheel)
	EVT_CLOSE(OBGUIFrame::OnClose)
	EVT_DROP_FILES(OBGUIFrame::OnDropFiles)
 END_EVENT_TABLE()

IMPLEMENT_APP(OBGUIApp)

// 'Main program' equivalent: the program execution "starts" here
bool OBGUIApp::OnInit()
{
	//Read in the stored window sizes and extensions previously used (or use the defaults)
	wxConfig config("OpenBabelGUI");
	wxSize size;
	size.SetWidth(config.Read("Width",1020));
	size.SetHeight(config.Read("Height",570));
	wxPoint position;
	position.x = config.Read("Left",10);
	position.y = config.Read("Top",30);

	//Save the full path of the help file - assumed to be in the starting working directory
	wxFileName help("OpenBabelGUI.html");
	help.MakeAbsolute();
	HelpFile = help.GetFullPath();

	// create the main application window
  OBGUIFrame *frame = new OBGUIFrame(_T("OpenBabelGUI"), position, size);

  // and show it (the frames, unlike simple controls, are not shown when
  // created initially)
  frame->Show(true);
  frame->SetInitialFocus(); //doesn't work!
  return true;
}

// ----------------------------------------------------------------------------
// main frame
// ----------------------------------------------------------------------------

OBGUIFrame::OBGUIFrame(const wxString& title, wxPoint position, wxSize size)
       : wxFrame(NULL,wxID_ANY,title, position,size), m_pGenOptsPanel(NULL)
{
	// set the frame icon
	SetIcon(wxICON(sample));
	DragAcceptFiles(true);

	wxConfig config("OpenBabelGUI");

	// create a menu bar
	fileMenu = new wxMenu;
	viewMenu = new wxMenu;
	helpMenu = new wxMenu;
	viewMenu->AppendCheckItem(ID_RESTRICTFORMATS, _T("Use &restricted set of formats"));
	viewMenu->Append(ID_SELFORMATS, _T("&Select set of formats"));
	viewMenu->AppendSeparator();
	viewMenu->AppendCheckItem(ID_SHOWAPIOPTIONS, _T("&API Options"),_T("e.g. errorlevel"));
	viewMenu->AppendCheckItem(ID_SHOWCONVOPTIONS, _T("&Conversion Options"));
	viewMenu->AppendCheckItem(ID_SHOWOBJOPTIONS1, _T("Chemical Object Options(&single char)"),
		_T("relating to molecules or reactions\ne.g. -h Add hydrogens"));
	viewMenu->AppendCheckItem(ID_SHOWOBJOPTIONS2, _T("Chemical Object Options(&multichar)"),
		_T("e.g. --addtotitle"));
	viewMenu->AppendCheckItem(ID_SHOWINOPTIONS, _T("&Input Format Options"));
	viewMenu->AppendCheckItem(ID_SHOWOUTOPTIONS, _T("&Output Format Options"));
	viewMenu->AppendSeparator();
	viewMenu->AppendCheckItem(ID_INWRAPPED, _T("Wrap input text (on restart)"));
	viewMenu->AppendCheckItem(ID_OUTWRAPPED, _T("Wrap output text (on restart)"));
	fileMenu->Append(ID_CONVERT, _T("&Convert"));
	fileMenu->Append(wxID_SAVE, _T("&Save Input Text As...\tCtrl+S"),
		_T("Brings up Save dialog"));
	fileMenu->Append(wxID_EXIT, _T("E&xit\tAlt-X"), _T("Quit OpenBabelGUI"));
	helpMenu->Append(wxID_HELP, _T("&Instructions"),
		_T("Opens default browser to view help file"));
	helpMenu->Append(wxID_ABOUT, _T("&About..."), _T("Show about dialog"));

	bool chk;
	config.Read("ShowConvOptions",&chk,true);
	viewMenu->Check(ID_SHOWCONVOPTIONS,chk);
	config.Read("ShowAPIOptions",&chk,true);
	viewMenu->Check(ID_SHOWAPIOPTIONS,chk);
	config.Read("ShowObjOptions1",&chk,true);
	viewMenu->Check(ID_SHOWOBJOPTIONS1,chk);
	config.Read("ShowObjOptions2",&chk,true);
	viewMenu->Check(ID_SHOWOBJOPTIONS2,chk);
	config.Read("ShowInOptions",&chk,true);
	viewMenu->Check(ID_SHOWINOPTIONS,chk);
	config.Read("ShowOutOptions",&chk,true);
	viewMenu->Check(ID_SHOWOUTOPTIONS,chk);
	config.Read("InWrapped",&chk,true);
	viewMenu->Check(ID_INWRAPPED,chk);
	config.Read("OutWrapped",&chk,true);
	viewMenu->Check(ID_OUTWRAPPED,chk);

	chk = m_ActiveFormats.ReadConfig(config);
	viewMenu->Check(ID_RESTRICTFORMATS,chk);

	// now append the freshly created menu to the menu bar...
	wxMenuBar *menuBar = new wxMenuBar();
	menuBar->Append(fileMenu, _T("&File"));
	menuBar->Append(viewMenu, _T("&View"));
	menuBar->Append(helpMenu, _T("&Help"));

	// ... and attach this menu bar to the frame
	SetMenuBar(menuBar);

//******************************************************
//**************** Controls (in tab order)**************

	wxPanel* panel = new wxPanel(this, wxID_ANY);

	m_pInFormat    = new wxChoice(panel,ID_INFORMAT,	wxDefaultPosition,wxDefaultSize,
		0, NULL);
	m_pInInfo      = new wxButton  (panel, ID_ININFO, wxT("Format Info"),
				wxDefaultPosition);
	m_pForceInFormat  = new wxCheckBox(panel,ID_INFORCEFORMAT,
				wxT("Use this format for all &input files (ignore file extensions)"));
	m_pInPath      = new wxStaticText(panel,wxID_STATIC,wxT(""),
				wxDefaultPosition,wxDefaultSize,wxST_NO_AUTORESIZE );
	m_pInFilename  = new CFilenames(panel, ID_INFILENAME, wxEmptyString,
				wxDefaultPosition,wxSize(150,20));
	m_pInFiles     = new wxButton  (panel, ID_INGETFILES, wxT("..."),
				wxDefaultPosition,wxSize(35,20));
	m_pInputHere   = new wxCheckBox(panel,ID_INPUTHERE,
				wxT("Input below (ignore input file)"));
	long notwrapped = viewMenu->IsChecked(ID_INWRAPPED) ? 0 : wxTE_DONTWRAP;
	m_pInText      = new wxTextCtrl( panel, ID_INTEXT, "",
        wxDefaultPosition, wxSize(195,200), wxTE_MULTILINE|wxTE_READONLY|notwrapped);

	m_pConvert     = new wxButton  (panel, ID_CONVERT, wxT("&CONVERT"),
				wxDefaultPosition,wxSize(80,40));
	MakeBold(m_pConvert);
	m_pConvert->SetToolTip("Do conversion (Alt C)");

	m_pOutFormat   = new wxChoice(panel,ID_OUTFORMAT,wxDefaultPosition,wxDefaultSize,
		0, NULL);
	m_pOutInfo     = new wxButton  (panel, ID_OUTINFO, wxT("Format Info"),
				wxDefaultPosition);
	m_pOutFilename = new wxTextCtrl(panel, ID_OUTFILENAME,wxEmptyString,
				wxDefaultPosition,wxSize(150,20));
	m_pOutFiles    = new wxButton  (panel, ID_OUTGETFILES, wxT("..."),
				wxDefaultPosition,wxSize(35,20));
	m_pNoOutFile   = new wxCheckBox(panel,ID_NOOUTFILE,
				wxT("Output below only (no output file)"));
	notwrapped = viewMenu->IsChecked(ID_OUTWRAPPED) ? 0 : wxTE_DONTWRAP;

	//Output windows: splitter with messages(clog) and output text(cout) 
	m_pSplitter = new wxSplitterWindow(panel,wxID_ANY,wxDefaultPosition,wxSize(1955,200));
	m_pMessages = new wxTextCtrl(m_pSplitter,ID_MESSAGES,wxEmptyString,
		wxDefaultPosition,wxDefaultSize,wxTE_MULTILINE|wxTE_READONLY|wxNO_BORDER);
	m_pMessages->SetToolTip("Message window. Drag the divider down to make bigger");
	m_pfixedFont = new wxFont(
		8, //int pointSize
		wxFONTFAMILY_MODERN, //wxFontFamily family
		wxFONTSTYLE_NORMAL, //style
		wxFONTWEIGHT_NORMAL); //wxFontWeight weight
	m_pMessages->SetFont(*m_pfixedFont);
	
	m_pOutText = new wxTextCtrl(m_pSplitter, ID_OUTTEXT, "",
        wxDefaultPosition, wxDefaultSize, wxTE_MULTILINE|wxTE_READONLY|notwrapped|wxNO_BORDER);
	
	m_pSplitter->SplitHorizontally(m_pMessages, m_pOutText, 20);
	m_pSplitter->SetMinimumPaneSize(20);
	int messize;
	config.Read("MessageWindowSize",&messize,20);
	m_pSplitter->SetSashPosition(messize);	

//******************************************************
//************** Layout with Sizers ********************

	topSizer		 = new wxBoxSizer(wxHORIZONTAL);
	InSizer      = new wxBoxSizer(wxVERTICAL);
	OutSizer     = new wxBoxSizer(wxVERTICAL);
	OptionsSizer = new wxBoxSizer(wxVERTICAL);
	wxBoxSizer *InFilesSizer = new wxBoxSizer(wxHORIZONTAL);
	wxBoxSizer *OutFilesSizer = new wxBoxSizer(wxHORIZONTAL);
	wxBoxSizer *InFormatSizer = new wxBoxSizer(wxHORIZONTAL);
	wxBoxSizer *OutFormatSizer = new wxBoxSizer(wxHORIZONTAL);
		

	InFormatSizer->Add(m_pInFormat,1,wxEXPAND);
	InFormatSizer->Add(m_pInInfo,0,wxLEFT,5);
	
	InFilesSizer->Add(m_pInFilename,1,wxEXPAND);
	InFilesSizer->Add(m_pInFiles,0,wxLEFT,5);

	OutFilesSizer->Add(m_pOutFilename,1,wxEXPAND);
	OutFilesSizer->Add(m_pOutFiles,0,wxLEFT,5);

	OutFormatSizer->Add(m_pOutFormat,1,wxEXPAND);
	OutFormatSizer->Add(m_pOutInfo,0,wxLEFT,5);
	wxStaticText* pStatic = new wxStaticText(panel,wxID_STATIC,wxT("        ---- INPUT FORMAT ----"));
	MakeBold(pStatic);
	InSizer->Add(pStatic);
	InSizer->Add(InFormatSizer,0, wxLEFT|wxTOP,5);
	InSizer->Add(m_pForceInFormat,0,wxLEFT|wxBOTTOM,5);
	InSizer->Add(m_pInPath,0,wxEXPAND|wxLEFT|wxTOP,5);
	InSizer->Add(InFilesSizer,0,wxEXPAND|wxALL,5);
	InSizer->Add(m_pInputHere,0,wxLEFT|wxBOTTOM,5);
	InSizer->Add(m_pInText,
      1,        // make vertically stretchable
      wxEXPAND| // make horizontally stretchable
      wxALL,    // and make border all around
      5 );     // set border width to 10

	m_pGenOptsPanel = new DynOptionswx(panel, OptionsSizer);
	m_pAPIOptsPanel = new DynOptionswx(panel, OptionsSizer);
	m_pConvOptsPanel = new DynOptionswx(panel, OptionsSizer);
	m_pInOptsPanel = new DynOptionswx(panel, OptionsSizer);
	m_pOutOptsPanel = new DynOptionswx(panel, OptionsSizer);

	pStatic = new wxStaticText(panel,wxID_STATIC,wxT("        ---- OUTPUT FORMAT ----"));
	MakeBold(pStatic);
	OutSizer->Add(pStatic);
	OutSizer->Add(OutFormatSizer,0, wxLEFT|wxTOP,5);

	//gap where checkbox used to be
	int chkwidth,chkheight;
	m_pForceInFormat->GetSize(&chkwidth, &chkheight);
	OutSizer->AddSpacer(chkheight+3);

	OutSizer->Add(new wxStaticText(panel,wxID_STATIC,wxT("Output file")),0,wxLEFT|wxTOP,5);
	OutSizer->Add(OutFilesSizer,0,wxEXPAND|wxALL,5);
	OutSizer->Add(m_pNoOutFile,0, wxLEFT|wxBOTTOM,5);
	OutSizer->Add(m_pSplitter, 1, wxEXPAND | wxALL, 5 );

	OptionsSizer->Add(m_pConvert, 0, wxALL|wxALIGN_CENTER_HORIZONTAL,10);
	OptionsSizer->Add(new wxStaticLine(panel),0,wxBOTTOM|wxEXPAND,5);

	topSizer->Add(InSizer,1,wxEXPAND);
	topSizer->Add(OptionsSizer,0.3,wxEXPAND);
	topSizer->Add(OutSizer,1,wxEXPAND);
	
	panel->SetSizer( topSizer );     // use the sizer for layout
	topSizer->Fit( panel );          // fit the dialog to the contents
	topSizer->SetSizeHints( panel ); // set hints to honor min size

	GetAvailableFormats();
	wxString inExt = config.Read("InExt","smi");
	wxString outExt = config.Read("OutExt","smi");
	SetChoice(m_pInFormat, inExt);
	SetChoice(m_pOutFormat,outExt);

	m_InFileBasePath = config.Read("InputPath", wxGetCwd());
	m_pInPath->SetLabel(ShortenedPath(m_InFileBasePath, *m_pInPath));
	::wxSetWorkingDirectory(m_InFileBasePath);

	config.Read("InputHere",&chk,false);
	m_pInputHere->SetValue(chk);
	ChangeInputHere(chk);

	//Display the options
	wxCommandEvent dum;
	OnChangeFormat(dum);
}

//******************************************************
//********** Event handlers ****************************

void OBGUIFrame::OnClose(wxCloseEvent& event)
{
  //Save the window size, in and out formats, and various options, for use next time.
	wxConfig config("OpenBabelGUI");

	int width, height, left, top;
	GetPosition(&left,&top);
	GetSize(&width, &height);
	config.Write("Left",left);
	config.Write("Top",top);
	config.Write("Width",width);
	config.Write("Height",height);
	config.Write("MessageWindowSize",m_pSplitter->GetSashPosition());

	config.Write("InputPath", m_InFileBasePath);
	
	wxString ext = m_pInFormat->GetStringSelection();
	int pos = ext.find_first_of(" \t-");
	config.Write("InExt", ext.substr(0,pos));
	
	ext = m_pOutFormat->GetStringSelection();
	pos = ext.find_first_of(" \t-");
	config.Write("OutExt", ext.substr(0,pos));

	config.Write("ShowConvOptions",viewMenu->IsChecked(ID_SHOWCONVOPTIONS));
	config.Write("ShowAPIOptions",viewMenu->IsChecked(ID_SHOWAPIOPTIONS));
	config.Write("ShowObjOptions1",viewMenu->IsChecked(ID_SHOWOBJOPTIONS1));
	config.Write("ShowObjOptions2",viewMenu->IsChecked(ID_SHOWOBJOPTIONS2));
	config.Write("ShowInOptions",viewMenu->IsChecked(ID_SHOWINOPTIONS));
	config.Write("ShowOutOptions",viewMenu->IsChecked(ID_SHOWOUTOPTIONS));
	config.Write("InWrapped",viewMenu->IsChecked(ID_INWRAPPED));
	config.Write("OutWrapped",viewMenu->IsChecked(ID_OUTWRAPPED));
	config.Write("InputHere", m_pInputHere->IsChecked());
	
	m_ActiveFormats.WriteConfig(config);
	config.Write("UseRestrictedFormats", viewMenu->IsChecked(ID_RESTRICTFORMATS));

	delete m_pGenOptsPanel;
	delete m_pAPIOptsPanel;
	delete m_pConvOptsPanel;
	delete m_pInOptsPanel;
	delete m_pOutOptsPanel;
	delete m_pfixedFont;
	this->Destroy();
}

void OBGUIFrame::OnQuit(wxCommandEvent& WXUNUSED(event))
{
	// true is to force the frame to close
   Close(true);
}

void OBGUIFrame::OnSaveInputText(wxCommandEvent& WXUNUSED(event))
{
	wxString filter = GetFilter(m_pInFormat) +_T("All files(*.*)|*.*||");
	wxFileDialog dialog(this, _T("Save input text"), m_InFileBasePath,
		_T("FromOB"), filter, wxSAVE|wxOVERWRITE_PROMPT);
	if(dialog.ShowModal() == wxID_OK)
	{
		wxString filename = dialog.GetPath();
		if(!filename.empty())
			m_pInText->SaveFile(filename);
	}
}

void OBGUIFrame::OnHelp(wxCommandEvent& WXUNUSED(event))
{
//  // 1) Get File Type
	wxMimeTypesManager mimeType;
    wxFileType *c_type;
	 c_type = mimeType.GetFileTypeFromExtension(_T("html"));
    if(!c_type) return ; //Couldn't find association
   
    // 2) Get Open Message
    wxString command;
   
		command = c_type->GetOpenCommand(((OBGUIApp*)wxTheApp)->HelpFile);
    if(!command) return; //No default program
   
    // 3) Execute message
    wxExecute(command);
   
    delete c_type;
//	::wxLaunchDefaultBrowser(_T("OpenBabelGUIHelp.html"));
	//wxExecute(_T("OpenBabelGUI.html"));
}
void OBGUIFrame::OnAbout(wxCommandEvent& WXUNUSED(event))
{
    wxString msg(_T("OpenBabelGUI (C) 2006 by Chris Morley\n\nThis program is part of\
the OpenBabel project,\nwhich is released under the GNU General Public License.\n\
See: http://openbabel.sourceforge.net/wiki/Main_Page\n\n\
For more detailed information see:\n\
http://openbabel.sourceforge.net/wiki/Windows_GUI\n\n \
OpenBabel version "));
	 msg << BABEL_VERSION;
	msg << " development snapshot March 2005";
    wxMessageBox(msg, _T("About OpenBabelGUI"), wxOK | wxICON_INFORMATION | wxCENTER, this);
}
///////////////////////////////////////////

void OBGUIFrame::OnSelectFormats(wxCommandEvent& event)
{
	if(m_ActiveFormats.SelectFormats())
	{
		//When formats have been selected always set the Use Restricted Set check
		viewMenu->Check(ID_RESTRICTFORMATS,true);
		GetAvailableFormats();
	}
}
void OBGUIFrame::OnRestrictFormats(wxCommandEvent& event)
{
	GetAvailableFormats();
}
void OBGUIFrame::OnConvert(wxCommandEvent& WXUNUSED(event))
{
	wxBusyCursor cw;
	m_pOutText->Clear();
	m_pMessages->Clear();

	//Default input is from input text box;
	std::stringstream ss(m_pInText->GetValue().c_str());
//	char ch = ss.peek();
	//Default output is cout, which is redirected to Output text box
	OBConversion Conv(&ss, &std::cout);

	int iSel = m_pInFormat->GetSelection();
	if((iSel)<0) return;
	int oSel = m_pOutFormat->GetSelection();
	if((oSel)<0) return;
	
	OBFormat* pInFormat = NULL;
	OBFormat* pOutFormat = (OBFormat*)m_pOutFormat->GetClientData(oSel);
	if(m_pForceInFormat->IsChecked() || m_pInputHere->IsChecked())
		pInFormat = (OBFormat*)m_pInFormat->GetClientData(iSel);
	Conv.SetInAndOutFormats( pInFormat,pOutFormat);

	DoOptions(Conv);

	//Setup output file
	std::string stdOutputFileName;
	wxString OutFileName = m_pOutFilename->GetValue();
	//If you are trying to output with no filename what you really wanted
	//was to output to the OutConsole
	if(m_pNoOutFile->IsChecked() || OutFileName.IsEmpty())
		m_pNoOutFile->SetValue(true);
	else
	{
		//If output filename has no path, use input path
		wxFileName fn(OutFileName);
		if(!fn.IsAbsolute())
		{
			fn.MakeAbsolute(m_InFileBasePath);
			OutFileName = fn.GetFullPath();
		}
		stdOutputFileName = OutFileName;
	}
	OBFormat* pOutFileFormat = Conv.FormatFromExt(OutFileName);
	if(!m_pNoOutFile->IsChecked() && pOutFileFormat && (pOutFileFormat!=pOutFormat))
		if(wxMessageBox("The output file name extension does not correspond \
with the output format.\nDo you wish to continue the conversion?",
			"Is the output filename correct?", wxOK | wxCANCEL)!=wxOK)
			return;

	//Setup input file
	std::vector<std::string> FileList, OutputFileList;
	if(!m_pInputHere->IsChecked())
		m_pInFilename->Expand(FileList);

  //redirect cerr & clog & cout
		wxStreamToTextRedirector cerrCapture(m_pMessages, &std::cerr);
		wxStreamToTextRedirector clogCapture(m_pMessages, &std::clog);
		wxStreamToTextRedirector coutCapture(m_pOutText);

	m_pOutText->Freeze();//Otherwise seems to be redrawn after each char from cout
	
	int count = Conv.FullConvert(FileList, stdOutputFileName, OutputFileList);
	
	m_pOutText->Thaw();

	//Get the last word on the first line of the description which should
	//be "molecules", "reactions", etc and remove the s if only one object converted
	std::string objectname(pOutFormat->TargetClassDescription());
	int pos = objectname.find('\n');
	if(count==1) --pos;
	objectname.erase(pos);
	pos = objectname.rfind(' ');
	if(pos==std::string::npos)
		pos=0;
	std::clog << count << objectname.substr(pos) << " converted";
  if(OutputFileList.size()>1)
  {
		std::clog << '\n' << OutputFileList.size() 
			<< " files output. The first is " << OutputFileList[0];
  }
	
	if(count>0 && !m_pNoOutFile->IsChecked())
	{
		//Read back file and add to output console
		m_pOutText->Clear();
		m_pOutText->LoadFile(OutputFileList[0].c_str());
	}
}

///////////////////////////////////////////
void OBGUIFrame::OnInFormatInfo(wxCommandEvent& WXUNUSED(event)) 
{
	int nSel=m_pInFormat->GetSelection();
	if(nSel<0) return;
	OBFormat* pFormat = (OBFormat*)m_pInFormat->GetClientData(nSel);
	wxString mes(pFormat->Description());
	wxString url(pFormat->SpecificationURL());
	if(!url.IsEmpty())
		mes += "\nURL for specification: " + url;
	wxMessageBox(mes, "Format info");		
}

///////////////////////////////////////////
void OBGUIFrame::OnOutFormatInfo(wxCommandEvent& WXUNUSED(event)) 
{
	int nSel=m_pOutFormat->GetSelection();
	if(nSel<0) return;
	OBFormat* pFormat = (OBFormat*)m_pOutFormat->GetClientData(nSel);
	wxString mes(pFormat->Description());
	wxString url(pFormat->SpecificationURL());
	if(!url.IsEmpty())
		mes += "\nURL for specification: " + url;
	wxMessageBox(mes, "Format info");		
}

void OBGUIFrame::OnGetInputFile(wxCommandEvent& WXUNUSED(event))
{
	wxFileDialog dialog(this,_T("Choose Input File"),m_InFileBasePath,_T(""),
			GetFilter(m_pInFormat) + InputFilterString,
			wxMULTIPLE | wxFILE_MUST_EXIST | wxHIDE_READONLY );
	if(dialog.ShowModal() == wxID_OK)
	{
//		m_pInFilename->Clear();
		wxArrayString filepatharray;
		dialog.GetPaths(filepatharray);
		DisplayInputFiles(filepatharray);
	}
}
void OBGUIFrame::DisplayInputFiles(wxArrayString filepatharray)
{
	int i, endsel=0, startsel=0;
	if(!wxGetKeyState(WXK_CONTROL))
	{
		m_pInFilename->Clear();
		wxFileName filenamewx(filepatharray[0]);
		m_InFileBasePath = filenamewx.GetVolume() + filenamewx.GetVolumeSeparator()
			+ filenamewx.GetPath(wxPATH_GET_SEPARATOR );
		m_pInPath->SetLabel(ShortenedPath(m_InFileBasePath, *m_pInPath));
		::wxSetWorkingDirectory(m_InFileBasePath);
	}
	else
	{
		if(!m_pInFilename->GetValue().empty())
			m_pInFilename->AppendText(';');
		startsel  = m_pInFilename->GetLastPosition();
	}

	for(i=0;i<filepatharray.GetCount();++i)
	{
		if(i==1)
			endsel = m_pInFilename->GetLastPosition();
		if(i!=0)
			m_pInFilename->AppendText(';');
		wxFileName fnamewx(filepatharray[filepatharray.GetCount()-1-i]);
		fnamewx.MakeRelativeTo(m_InFileBasePath);
		m_pInFilename->AppendText(fnamewx.GetFullPath());

		if(wxGetKeyState(WXK_CONTROL) && i==0)
			endsel = m_pInFilename->GetLastPosition();
	}
	m_pInFilename->SetSelection(startsel,endsel);
	
	m_pInText->Clear();
	m_pInText->LoadFile(filepatharray[filepatharray.GetCount()-1]);
	m_pInFilename->SetFocus();
}

void OBGUIFrame::OnGetOutputFile(wxCommandEvent& WXUNUSED(event))
{
	wxFileDialog dialog(this,_T("Choose Output File"),_T(""),_T(""),
			GetFilter(m_pOutFormat) + OutputFilterString,	
			wxSAVE | wxOVERWRITE_PROMPT);
	if(dialog.ShowModal() == wxID_OK)
	{
		wxString filepath = dialog.GetPath();
		m_pOutFilename->Clear();
		m_pOutFilename->AppendText(filepath);
	}
}

void OBGUIFrame::OnOutFileNameUpdate(wxUpdateUIEvent& event)
{
	event.Enable(!m_pNoOutFile->IsChecked());	
}
void OBGUIFrame::OnInFileNameUpdate(wxUpdateUIEvent& event)
{
	event.Enable(!m_pInputHere->IsChecked());	
}

void OBGUIFrame::OnChangeInputHere(wxCommandEvent& event)
{
	ChangeInputHere(event.IsChecked());
}

void OBGUIFrame::ChangeInputHere(bool chk)
{
	wxColour bg = chk ? wxColour(250,255,210) : wxNullColour; 
	m_pInText->SetEditable(chk);
	m_pInText->SetBackgroundColour(bg);
	m_pInText->Refresh();
	m_pInText->SetFocus();
	m_pInFilename->Enable(chk);
}

void OBGUIFrame::OnChangeFormat(wxCommandEvent& WXUNUSED(event))
{
	//Display the options
	m_pAPIOptsPanel->Clear();
	m_pConvOptsPanel->Clear();
	m_pGenOptsPanel->Clear();
	m_pInOptsPanel->Clear();
	m_pOutOptsPanel->Clear();

	if(viewMenu->IsChecked(ID_SHOWAPIOPTIONS))
	{
		OBFormat* pAPI= OBConversion::FindFormat("obapi");
		if(pAPI)
			m_pAPIOptsPanel->Construct(pAPI->Description());
	}

	if(viewMenu->IsChecked(ID_SHOWCONVOPTIONS))
		m_pConvOptsPanel->Construct(OBConversion::Description());

	OBFormat* pInFormat = (OBFormat*)m_pInFormat->GetClientData(m_pInFormat->GetSelection());
	OBFormat* pOutFormat = (OBFormat*)m_pOutFormat->GetClientData(m_pOutFormat->GetSelection());
	if(!pInFormat || !pOutFormat)
		return;
	if(viewMenu->IsChecked(ID_SHOWOBJOPTIONS1)) //Single char
		m_pGenOptsPanel->Construct(pOutFormat->TargetClassDescription(),NULL,1);
	if(viewMenu->IsChecked(ID_SHOWOBJOPTIONS2)) //Multi char
		m_pGenOptsPanel->Construct(pOutFormat->TargetClassDescription(),NULL,2);

	if(viewMenu->IsChecked(ID_SHOWINOPTIONS))
	{
		if(pInFormat && !m_pInOptsPanel->Construct(pInFormat->Description(),"input"))
			m_pInOptsPanel->Construct(pInFormat->Description(),"read");//try again
	}
	if(viewMenu->IsChecked(ID_SHOWOUTOPTIONS))
	{
		if(pOutFormat && !m_pOutOptsPanel->Construct(pOutFormat->Description(),"output"))
			m_pOutOptsPanel->Construct(pOutFormat->Description());//try again without the "output"
	}
	topSizer->Layout();
}

void OBGUIFrame::SetInitialFocus()
{
	m_pInFilename->SetFocus();
}

//*******************************************************
//********** Local functions ****************************

/**	 Set the selection of the control to the string which starts with
     the file extension (case independent). The extension .gz is ignored.
		 The string parameter can be the extension itself. 
		 Returns false if no match found.
 **/
bool OBGUIFrame::SetChoice(wxChoice* pChoice, const wxString& FileName)
{
	wxString ext(FileName);
	int pos = FileName.rfind('.');
	if(pos!=-1)
		wxString ext = FileName.substr(pos+1);
	if(FileName.substr(pos)==".gz")
		pos = FileName.rfind('.',pos-1);

	for(int iSel=0;iSel<pChoice->GetCount();iSel++)
	{
		wxString txt(pChoice->GetString(iSel).substr(0,ext.length()));
		if(txt.MakeUpper()==ext.MakeUpper())
		{
			pChoice->SetSelection(iSel);
			return true;
		}
	}
	return false;
}
///////////////////////////////////////////////////
wxString OBGUIFrame::GetFilter(wxChoice* pChoice)
{
	//Uses text from the window (input or output combo) to construct filter string
	wxString  txt = pChoice->GetStringSelection();
	int pos1 = txt.find(" ");
	int pos2 = txt.find_first_not_of( " :-\t", pos1);
	return txt.substr(pos2) + " (*." + txt.substr(0, pos1) + ")|*." + txt.substr(0, pos1) + "|";
}
/////////////////////////////////////////////////////
wxString OBGUIFrame::ShortenedPath(const wxString& path, const wxWindow& wnd, int wndwidth)
{
	/* Tries to leave first and last folder names, so that with at least
	four separators for the path to be shortened like:

	c:\My Documents\projects\wxWindows\test\readme.txt
	to
	c:\My Documents\...\wxWindows\test\readme.txt
	or
	c:\My Documents\...\test\readme.txt
	
	If the width is still too large or there are 3 or fewer separators,
	the path is returned unshortened.:
	*/
	int txtwidth, txtheight, wndheight;
	wnd.GetTextExtent(path, &txtwidth, &txtheight);
	if(wndwidth<0)
		wnd.GetClientSize(&wndwidth,&wndheight);
	if(txtwidth <= wndwidth) return path;

	// Set pos1 at the second separator or return unchanged string if not possible
	size_t pos1 = path.find_first_of("/\\");
	if(pos1==-1 || pos1+1 == path.length()) return path; 
	pos1= path.find_first_of("/\\", pos1+1);
	if(pos1==-1 || pos1+1 == path.length()) return path; 

	wxString tpath(path);
	size_t pos2 = pos1;
	while(txtwidth > wndwidth)
	{
		pos2= tpath.find_first_of("/\\", pos2+1);
		if(pos2==-1 || pos2+1 == tpath.length()) return path; 
		//pos2 now has next separator
		//Ensure that there is at least one further directory
		if(tpath.find_first_of("/\\", pos2+1)==-1) return path;

		tpath.replace(pos1+1, pos2-pos1-1,"...");
		pos2=pos1+4;
		wnd.GetTextExtent(tpath, &txtwidth, &txtheight);
	}
	return tpath;
}

void OBGUIFrame::DoOptions(OpenBabel::OBConversion& Conv)
{
	// Is a API directive, e.g.---errorlevel
	//Send to the pseudoformat "obapi" (without any leading -)
	OBFormat* pAPI= OpenBabel::OBConversion::FindFormat("obapi");
	if(pAPI)
	{
		OBConversion apiConv;
		if(m_pAPIOptsPanel->SetOptions(apiConv, OBConversion::GENOPTIONS))
		{
			apiConv.SetOutFormat(pAPI);
			apiConv.Write(NULL);
		}
	}
	m_pGenOptsPanel->SetOptions(Conv, OBConversion::GENOPTIONS);
	m_pConvOptsPanel->SetOptions(Conv, OBConversion::GENOPTIONS);
	m_pInOptsPanel->SetOptions(Conv, OBConversion::INOPTIONS);
	m_pOutOptsPanel->SetOptions(Conv, OBConversion::OUTOPTIONS);
}

void OBGUIFrame::GetAvailableFormats()
{
	//Get data on available formats and add to comboboxes and to filter string
	OBConversion dummy; //needed for OBConversion to load format classes
	const char* str=NULL;
	OBFormat* pFormat;
	int nInSel=0,nOutSel=0;
	m_pInFormat->Clear();
	m_pOutFormat->Clear();
	InputFilterString="All Chemical Formats|*.";
	OutputFilterString = InputFilterString;
	m_ActiveFormats.Clear();
	Formatpos pos;
	while(OBConversion::GetNextFormat(pos,str,pFormat))
	{
		if(!str || !pFormat) break; //no formats available
		if((pFormat->Flags() & NOTWRITABLE) && (pFormat->Flags() & NOTREADABLE))
			continue;

		wxString txt(str);
		int pos = txt.find('[');
		if(pos!=wxString::npos)
			txt.erase(pos);
		//Check whether is in restricted set of formats, if this is in use
		if(!m_ActiveFormats.Add(txt) && viewMenu->IsChecked(ID_RESTRICTFORMATS))
			continue;
		int n;
		if(!(pFormat->Flags() & NOTREADABLE))
		{
			n = m_pInFormat->Append(txt,pFormat);
			InputFilterString+=txt.Left(txt.Find(" "));
			InputFilterString+=";*.";
		}
		if(!(pFormat->Flags() & NOTWRITABLE))
		{
			n = m_pOutFormat->Append(txt,pFormat);
			OutputFilterString+=txt.Left(txt.Find(" "));
			OutputFilterString+=";*.";
		}		
	}
	if(m_pInFormat->GetCount()==0)
	{
		m_pInFormat->Append("No input format in the selected set");
		m_pInFormat->SetClientData(0,NULL);
	}
	if(m_pOutFormat->GetCount()==0)
	{
		m_pOutFormat->Append("No output format in the selected set");
		m_pOutFormat->SetClientData(0,NULL);
	}
	m_pInFormat->SetSelection(nInSel);
	m_pOutFormat->SetSelection(nOutSel);

	InputFilterString = InputFilterString.Left(InputFilterString.Length()-3); //remove unneeded ;*.
	OutputFilterString = OutputFilterString.Left(OutputFilterString.Length()-3); //remove unneeded ;*.
	InputFilterString+="|AllFiles(*.*)|*.*||";
	OutputFilterString+="|AllFiles(*.*)|*.*||";
}	

void OBGUIFrame::DisplayInFile(wxString filename)
{
	m_pInText->Clear();
	wxFileName fn(filename);
	fn.MakeAbsolute(m_InFileBasePath);
	m_pInText->LoadFile(fn.GetFullPath());
}

void OBGUIFrame::MakeBold(wxWindow* pWnd)
{
	wxFont font = pWnd->GetFont();
	font.SetWeight(wxFONTWEIGHT_BOLD );
	pWnd->SetFont(font);
}

void OBGUIFrame::OnDropFiles(wxDropFilesEvent& event)
{
//	m_pInFilename->Clear();	
	DisplayInputFiles(wxArrayString(event.GetNumberOfFiles(),event.GetFiles()));
}

//**********************************************
BEGIN_EVENT_TABLE(CFilenames,wxTextCtrl)
	EVT_LEFT_DCLICK(CFilenames::OnDblClick)
	EVT_CHAR(CFilenames::OnKeyPress)
END_EVENT_TABLE()

void CFilenames::OnDblClick(wxMouseEvent& event)
{
	//extract double-clicked filename
	OBGUIFrame* frame = static_cast<OBGUIFrame*>(GetParent()->GetParent());
	wxASSERT(frame);
	frame->DisplayInFile(SelectFilename());
}
/// Highlights and returns the filename at the current cursor position
wxString CFilenames::SelectFilename()
{
	wxString fname(GetValue());
	int endsel;
	long n = GetInsertionPoint();
	int pos = fname.find(';',n);
	endsel=pos;
	if(pos!=-1)
		fname.erase(pos);
	else
		endsel=GetLastPosition();
	pos = fname.rfind(';',n);
	if(pos!=-1)
		fname.erase(0,pos+1);
	SetSelection(pos+1,endsel);

	fname.Trim();
	fname.Trim(false);
	return fname;
}

void OBGUIFrame::OnMouseWheel(wxMouseEvent& event)
{
	int delta = event.GetWheelRotation() / event.GetWheelDelta();
	m_pInFilename->ToNextFile(-delta); //is direction correct?
}
void CFilenames::OnKeyPress(wxKeyEvent& event)
{
	int delta=1;
	OBGUIFrame* frame = static_cast<OBGUIFrame*>(GetParent()->GetParent());
	switch (event.GetKeyCode())
	{
	case 313: 
	case WXK_PAGEDOWN: //why is code not correct?
		delta=-1;
	case 312: 
	case WXK_PAGEUP:
		ToNextFile(delta);
		break;
	case WXK_TAB:
		if(event.ShiftDown())
			ToNextFile(-1);
		else
			ToNextFile(+1);
		break;

	case WXK_RETURN:
		if(event.ShiftDown())
			SetValue(nameWithWildcard);
		else
		{
			nameWithWildcard = GetValue();
			if(nameWithWildcard.find_first_of("*?")==-1)
			{
				wxCommandEvent dum;
				frame->OnConvert(dum);
			}
			else
			{
				//expand wildcards in each filename
				std::vector<std::string> filelist;
				int count = Expand(filelist);
				wxString mes;
				mes << count << _T(" files found");
				frame->DisplayMessage(mes);

				Clear();
				int i, endsel=0;
				for(i=0;i<filelist.size();++i)
				{
					if(i!=0)
					{
						if(i==1)
							endsel = GetLastPosition();
						AppendText(';');
					}
					wxString name(filelist[i].c_str());
					wxFileName fnamewx(name);
					fnamewx.MakeRelativeTo(frame->GetInFileBasePath());
					AppendText(fnamewx.GetFullPath());
				}
				SetSelection(0,endsel);
			}
		}
		break;
	default:
		event.Skip();
	}
}

int  CFilenames::Expand(std::vector<std::string>& filelist)
{
	//Adds full path names of all input files to filelist, expanding wildcards they are present. 
	OBGUIFrame* frame = static_cast<OBGUIFrame*>(GetParent()->GetParent());
	int namestart=0, count=0;
	wxString txt(GetValue());
	do //for each filename
	{
		int nameend = txt.find(';',namestart);
		wxString name = txt.substr(namestart, nameend-namestart);
		name.Trim().Trim(false);
		if(name.IsEmpty())
			break;
		wxFileName fn(name);
		namestart=nameend+1; // 0 at end
		fn.MakeAbsolute(frame->GetInFileBasePath());
		count += DLHandler::findFiles(filelist, std::string(fn.GetFullPath()));
	}while(namestart);

	return count;
}
bool CFilenames::ToNextFile(int delta)
{
	wxString fname(GetValue());
	if (fname.IsEmpty())
		return false;
	long n = GetInsertionPoint();
	int pos;
	if(delta>0)
	{
		pos = fname.find(';',n);
		if(pos==-1)
			pos=GetLastPosition();
		else
			++pos;
	}
	else
		pos = fname.rfind(';',n);
	if(pos!=-1)
	{
		SetInsertionPoint(pos);
		OBGUIFrame* frame = static_cast<OBGUIFrame*>(GetParent()->GetParent());
		wxString nxtname = SelectFilename();
		if(nxtname.IsEmpty())
			return false;
		frame->DisplayInFile(nxtname);
		return true;
	}
	else
		return false;
}
//***********************************************



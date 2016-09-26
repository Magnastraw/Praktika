///////////////////////////////////////////////////////////////////////////////
// main.cpp
// ========
// main driver
//
//  AUTHOR: Song Ho Ahn (song.ahn@gmail.com)
// CREATED: 2006-06-28
// UPDATED: 2016-05-11
///////////////////////////////////////////////////////////////////////////////

#define WIN32_LEAN_AND_MEAN             // exclude rarely-used stuff from Windows headers

#include <windows.h>
#include "Window.h"
#include "DialogWindow.h"
#include "ControllerGL.h"
#include "ControllerForm.h"
#include "ControllerMain.h"
#include "ModelGL.h"
#include "ViewGL.h"
#include "ViewForm.h"
#include "resource.h"


// function declarations
int mainMessageLoop(HACCEL hAccelTable=0);




///////////////////////////////////////////////////////////////////////////////
// main function of a windows application
///////////////////////////////////////////////////////////////////////////////
int WINAPI WinMain(HINSTANCE hInst, HINSTANCE hPrevInst, LPSTR cmdArgs, int cmdShow)
{

	Win::ControllerMain mainCtrl;
    // create "controller" component by specifying what are "model" and "view"

	Win::Window mainWin(hInst, L"Elements", 0, &mainCtrl);
	mainWin.setWidth(600);
	mainWin.setHeight(560);
	mainWin.setWindowStyleEx(WS_EX_WINDOWEDGE);
	mainWin.create();
	
	// instantiate model and view components, so "controller" component can reference them
	ModelGL model;
	Win::ViewGL view;   // under "Win" namespace because it is Windows specific view component.
	Win::ControllerGL glCtrl(&model, &view);
	
    // create window with given controller
    Win::Window glWin(hInst, L"glwin", mainWin.getHandle(), &glCtrl);
	glWin.setClassStyle(CS_OWNDC);
    glWin.setWindowStyle(WS_CHILD | WS_VISIBLE | WS_CLIPSIBLINGS | WS_CLIPCHILDREN);
    glWin.setWidth(600);
    glWin.setHeight(400);
    glWin.create();

    // create sub dialog to hold some controls (buttons, editbox, etc.)
    Win::ViewForm formView;                             // view component of dialog
    Win::ControllerForm formCtrl(&model, &formView);    // controller component of dialog
    Win::DialogWindow formDialog(hInst, IDD_FORMVIEW, mainWin.getHandle(), &formCtrl);
    formDialog.create();


	mainCtrl.setGLHandle(glWin.getHandle());
	mainCtrl.setFormHandle(formDialog.getHandle());

	::SetWindowPos(formDialog.getHandle(), 0, 0, 600, 300, 200, SWP_NOZORDER);


	// compute height of all sub-windows
	int height = 0;
	RECT rect;
	::GetWindowRect(glWin.getHandle(), &rect);      // get size of glWin
	height += rect.bottom - rect.top;
	::GetWindowRect(formDialog.getHandle(), &rect);   // get size of glDialog
	height += rect.bottom - rect.top;

	// resize main window, so all sub windows are fit into the client area of main window
	DWORD style = (DWORD)::GetWindowLongPtr(mainWin.getHandle(), GWL_STYLE);       // get current window style
	DWORD styleEx = (DWORD)::GetWindowLongPtr(mainWin.getHandle(), GWL_EXSTYLE);   // get current extended window style
	rect.left = 0;
	rect.right = 600;
	rect.top = 0;
	rect.bottom = height;
	::AdjustWindowRectEx(&rect, style, TRUE, styleEx);
	::SetWindowPos(mainWin.getHandle(), 0, 100, 100, rect.right - rect.left, rect.bottom - rect.top, SWP_NOZORDER);
	glWin.show();
	formDialog.show();
	mainWin.show();

    // main message loop //////////////////////////////////////////////////////
    int exitCode;
	HACCEL hAccelTable = 0;
    exitCode = mainMessageLoop(hAccelTable);

    return exitCode;
}



///////////////////////////////////////////////////////////////////////////////
// main message loop
///////////////////////////////////////////////////////////////////////////////
int mainMessageLoop(HACCEL hAccelTable)
{
	HWND activeHandle;
    MSG msg;

    while(::GetMessage(&msg, 0, 0, 0) > 0)  // loop until WM_QUIT(0) received
    {
		activeHandle = GetActiveWindow();
		if (::GetWindowLongPtr(activeHandle, GWL_EXSTYLE) & WS_EX_CONTROLPARENT) // WS_EX_CONTROLPARENT is automatically added by CreateDialogBox()
		{
			if (::IsDialogMessage(activeHandle, &msg))
				continue;   // message handled, back to while-loop
		}

        // now, handle window messages
        if(!::TranslateAccelerator(msg.hwnd, hAccelTable, &msg))
        {
            ::TranslateMessage(&msg);
            ::DispatchMessage(&msg);
        }
    }
    return (int)msg.wParam;                 // return nExitCode of PostQuitMessage()
}

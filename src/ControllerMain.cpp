///////////////////////////////////////////////////////////////////////////////
// ControllerMain.cpp
// ==================
// Derived Controller class for main window
//
//  AUTHOR: Song Ho Ahn (song.ahn@gmail.com)
// CREATED: 2006-07-09
// UPDATED: 2014-01-01
///////////////////////////////////////////////////////////////////////////////

#include <windows.h>
#include <commctrl.h>                   // common controls
#include <sstream>
#include "ControllerMain.h"
#include "resource.h"
using namespace Win;


// handle events(messages) on all child windows that belong to the parent window.
// For example, close all child windows when the parent got WM_CLOSE message.
// lParam can be used to specify a event or message.
bool CALLBACK enumerateChildren(HWND childHandle, LPARAM lParam);



ControllerMain::ControllerMain() : glHandle(0), formHandle(0)
{
}



int ControllerMain::command(int id, int cmd, LPARAM msg)
{

	return 0;
}



int ControllerMain::close()
{

	// close all child windows first
	::EnumChildWindows(handle, (WNDENUMPROC)enumerateChildren, (LPARAM)WM_CLOSE);

	::DestroyWindow(handle);    // close itself
	return 0;
}



int ControllerMain::destroy()
{
	::PostQuitMessage(0);       // exit the message loop

	return 0;
}



int ControllerMain::create()
{
	return 0;
}

int ControllerMain::size(int width, int height, WPARAM wParam)
{
	RECT rect;

	// get height of glDialog
	::GetWindowRect(formHandle, &rect);
	int formHeight = rect.bottom - rect.top;

	int glHeight = height - formHeight;
	::SetWindowPos(glHandle, 0, 0, 0, width, glHeight, SWP_NOZORDER);
	::SetWindowPos(formHandle, 0, 0, glHeight, width, formHeight, SWP_NOZORDER);
	::InvalidateRect(formHandle, 0, TRUE);      // force to repaint



	return 0;
}


///////////////////////////////////////////////////////////////////////////////
// enumerate all child windows
///////////////////////////////////////////////////////////////////////////////
bool CALLBACK enumerateChildren(HWND handle, LPARAM lParam)
{
	if (lParam == WM_CLOSE)
	{
		::SendMessage(handle, WM_CLOSE, 0, 0);      // close child windows
	}

	return true;
}

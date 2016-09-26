///////////////////////////////////////////////////////////////////////////////
// ControllerForm.h
// ================
// Derived Controller class for FormView window
//
//  AUTHOR: Song Ho Ahn (song.ahn@gamil.com)
// CREATED: 2016-05-10
// UPDATED: 2016-05-11
///////////////////////////////////////////////////////////////////////////////

#ifndef WIN_CONTROLLER_FORM_H
#define WIN_CONTROLLER_FORM_H

#include "Controller.h"
#include "ViewForm.h"
#include "ModelGL.h"


namespace Win
{
    class ControllerForm : public Controller
    {
    public:
        ControllerForm(ModelGL* model, ViewForm* view); // ctor with params
        ~ControllerForm() {};                           // dtor

        int create();                               // for WM_CREATE
        int command(int id, int cmd, LPARAM msg);   // for WM_COMMAND

        void setParentWindowHandle(HWND handle)     { parentHandle = handle; }

    private:
        ModelGL* model;                             //
        ViewForm* view;                             //
        HWND parentHandle;                          // handle of parent window

		void openFileDialog();
		void saveFileDialog();
    };
}

#endif

///////////////////////////////////////////////////////////////////////////////
// ViewForm.h
// ==========
// View component of FormView window
//
//  AUTHORL Song Ho Ahn (song.ahn@gmail.com)
// CREATED: 2016-05-10
// UPDATED: 2016-05-10
///////////////////////////////////////////////////////////////////////////////

#ifndef VIEW_FORM_H
#define VIEW_FORM_H

#include <windows.h>
#include "Controls.h"

namespace Win
{
    class ViewForm
    {
    public:
        ViewForm();                             // ctor
        ~ViewForm();                            // dtor

        void initControls(HWND handle);         // init all controls
        float getValue(EditBox);


		EditBox editE;
		EditBox editMu;
		EditBox editNx;
		EditBox editNy;
		EditBox editNz;
		EditBox editA;
		EditBox editB;
		EditBox editC;
		EditBox editP;
		EditBox editdX;
		EditBox editdY;
		EditBox editdZ;
		EditBox editCNx;
		EditBox editCNy;
		EditBox editCNz;
		CheckBox checkX;
		CheckBox checkY;
		CheckBox checkZ;
		CheckBox checkP;
    protected:

    private:
        Button  buttonBuild;
		Button  buttonNextIter;
		Button  buttonSave;
		Button  buttonLoad;
		//EditBox editE;
		//EditBox editMu;
		//EditBox editNx;
		//EditBox editNy;
		//EditBox editNz;
    };
}

#endif

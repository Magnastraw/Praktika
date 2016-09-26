///////////////////////////////////////////////////////////////////////////////
// ViewForm.cpp
// ==========
// View component of FormView window
//
//  AUTHORL Song Ho Ahn (song.ahn@gmail.com)
// CREATED: 2016-05-10
// UPDATED: 2016-05-11
///////////////////////////////////////////////////////////////////////////////

#include "ViewForm.h"
#include "resource.h"
using namespace Win;

///////////////////////////////////////////////////////////////////////////////
// default ctor
///////////////////////////////////////////////////////////////////////////////
ViewForm::ViewForm()
{
}


///////////////////////////////////////////////////////////////////////////////
// default dtor
///////////////////////////////////////////////////////////////////////////////
ViewForm::~ViewForm()
{
}



///////////////////////////////////////////////////////////////////////////////
// initialize all controls
///////////////////////////////////////////////////////////////////////////////
void ViewForm::initControls(HWND handle)
{
    // set all controls
    buttonBuild.set(handle, IDC_BUILD);
	buttonNextIter.set(handle, IDC_NEXTITERATION);
	buttonSave.set(handle, IDC_BUTTON_SAVE);
	buttonLoad.set(handle, IDC_BUTTON_LOAD);
    editE.set(handle, IDC_EDIT_E);
	editMu.set(handle, IDC_EDIT_MU);
	editCNx.set(handle, IDC_EDIT_CNX);
	editCNy.set(handle, IDC_EDIT_CNY);
	editCNz.set(handle, IDC_EDIT_CNZ);
	editNx.set(handle, IDC_EDIT_NX);
	editNy.set(handle, IDC_EDIT_NY);
	editNz.set(handle, IDC_EDIT_NZ);
	editA.set(handle, IDC_EDIT_A);
	editB.set(handle, IDC_EDIT_B);
	editC.set(handle, IDC_EDIT_C);
	editP.set(handle, IDC_EDIT_P);
	editdX.set(handle, IDC_EDIT_DX);
	editdY.set(handle, IDC_EDIT_DY);
	editdZ.set(handle, IDC_EDIT_DZ);
	checkX.set(handle, IDC_CHECK_X);
	checkY.set(handle, IDC_CHECK_Y);
	checkZ.set(handle, IDC_CHECK_Z);
	checkP.set(handle, IDC_CHECK_P);

}


///////////////////////////////////////////////////////////////////////////////
// return value of the editbox as float
///////////////////////////////////////////////////////////////////////////////
float ViewForm::getValue(EditBox edit)
{
    const int MAX_LENGTH = 1024;
    wchar_t wcharBuff[MAX_LENGTH];
    char charBuff[MAX_LENGTH];

    edit.getText(wcharBuff, MAX_LENGTH);
    wcstombs(charBuff, wcharBuff, MAX_LENGTH);  // copy wchar_t as char
    charBuff[MAX_LENGTH-1] = '\0';              // in case when source exceeded max length

    return (float)atof(charBuff);
}


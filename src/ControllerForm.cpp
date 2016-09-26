///////////////////////////////////////////////////////////////////////////////
// ControllerForm.cpp
// ==================
// Derived Controller class for FormView window
//
//  AUTHOR: Song Ho Ahn (song.ahn@gamil.com)
// CREATED: 2016-05-10
// UPDATED: 2016-05-11
///////////////////////////////////////////////////////////////////////////////

#include "ControllerForm.h"
#include "resource.h"
using namespace Win;

int CNx, CNy, CNz, Nx, Ny, Nz;
double A_, B_, C_;
double Pp, DX, DY, DZ;
bool move_X, move_Y, move_Z;

const int WCHAR_MAX_COUNT = 16;                                 // max number of string buffers
const int WCHAR_MAX_LENGTH = 1024;                              // max string length per buffer
static wchar_t wchar_wideStr[WCHAR_MAX_COUNT][WCHAR_MAX_LENGTH];// circular buffer for wchar_t*
static char wchar_str[WCHAR_MAX_COUNT][WCHAR_MAX_LENGTH];       // circular buffer for char*
static int wchar_indexWchar = 0;                                // current index of circular buffer
static int wchar_indexChar = 0;                                 // current index of circular buffer
///////////////////////////////////////////////////////////////////////////////
// default contructor
///////////////////////////////////////////////////////////////////////////////
ControllerForm::ControllerForm(ModelGL* model, ViewForm* view) : model(model), view(view)
{
    parentHandle = 0;
}



///////////////////////////////////////////////////////////////////////////////
// handle WM_CREATE
///////////////////////////////////////////////////////////////////////////////
int ControllerForm::create()
{
    // initialize all controls
    view->initControls(handle);
    return 0;
}



///////////////////////////////////////////////////////////////////////////////
// handle WM_COMMAND
///////////////////////////////////////////////////////////////////////////////
int ControllerForm::command(int id, int command, LPARAM msg)
{
    switch(id)
    {
    case IDC_BUILD:
        if(command == BN_CLICKED)
        {
			model->createObject(CNx ,CNy ,CNz ,A_ ,B_ ,C_, Nx, Ny, Nz);
			model->initParametrs();
        }
        break;
	case IDC_NEXTITERATION:
		if (command == BN_CLICKED)
		{
			model->setP(Pp);
			model->updateTransform();
		}
		break;
	case IDC_BUTTON_SAVE:
		if (command == BN_CLICKED)
		{
			saveFileDialog();
		}
		break;
	case IDC_BUTTON_LOAD:
		if (command == BN_CLICKED)
		{
			openFileDialog();
			model->runFlag = true;
		}
		break;
    case IDC_EDIT_E:
        if(command == EN_CHANGE)
        {
            float E = view->getValue(view->editE);
			model->setE(E);
        }
        break;
	case IDC_EDIT_MU:
		if (command == EN_CHANGE)
		{
			float mu = view->getValue(view->editMu);
			model->setMu(mu);
		}
		break;
	case IDC_EDIT_CNX:
		if (command == EN_CHANGE)
		{
			CNx = (int)view->getValue(view->editCNx);
		}
		break;
	case IDC_EDIT_CNY:
		if (command == EN_CHANGE)
		{
			CNy = (int)view->getValue(view->editCNy);
		}
		break;
	case IDC_EDIT_CNZ:
		if (command == EN_CHANGE)
		{
			CNz = (int)view->getValue(view->editCNz);
		}
		break;
	case IDC_EDIT_NX:
		if (command == EN_CHANGE)
		{
			Nx = (int)view->getValue(view->editNx);
		}
		break;
	case IDC_EDIT_NY:
		if (command == EN_CHANGE)
		{
			Ny = (int)view->getValue(view->editNy);
		}
		break;
	case IDC_EDIT_NZ:
		if (command == EN_CHANGE)
		{
			Nz = (int)view->getValue(view->editNz);
		}
		break;
	case IDC_EDIT_A:
		if (command == EN_CHANGE)
		{
			A_ = view->getValue(view->editA);
		}
		break;
	case IDC_EDIT_B:
		if (command == EN_CHANGE)
		{
			B_ = view->getValue(view->editB);
		}
		break;
	case IDC_EDIT_C:
		if (command == EN_CHANGE)
		{
			C_ = view->getValue(view->editC);
		}
		break;
	case IDC_EDIT_P:
		if (command == EN_CHANGE)
		{
			Pp = view->getValue(view->editP);
			model->setP(Pp);
		}
		break;
	case IDC_EDIT_DX:
		if (command == EN_CHANGE)
		{
			DX = view->getValue(view->editdX);
			model->setdX(DX);
		}
		break;
	case IDC_EDIT_DY:
		if (command == EN_CHANGE)
		{
			DY = view->getValue(view->editdY);
			model->setdY(DY);
		}
		break;
	case IDC_EDIT_DZ:
		if (command == EN_CHANGE)
		{
			DZ = view->getValue(view->editdZ);
			model->setdZ(DZ);
		}
		break;
	case IDC_CHECK_X:
		if (view->checkX.isChecked()) {
			model->setMoveX(true);
		}
		else if (!view->checkX.isChecked()) {
			model->setMoveX(false);
		}
		break;
	case IDC_CHECK_Y:
		if (view->checkY.isChecked()) {
			model->setMoveY(true);
		}
		else if (!view->checkY.isChecked()) {
			model->setMoveY(false);
		}
		break;
	case IDC_CHECK_Z:
		if (view->checkZ.isChecked()) {
			model->setMoveZ(true);
		}
		else if (!view->checkZ.isChecked()) {
			model->setMoveZ(false);
		}
		break;
	case IDC_CHECK_P:
		if (view->checkP.isChecked()) {
			model->setP_bool(true);
		}
		else if (!view->checkP.isChecked()) {
			model->setP_bool(false);
		}
		break;
    }
    return 0;
}

///////////////////////////////////////////////////////////////////////////////
// convert wchar_t* string to char* string
///////////////////////////////////////////////////////////////////////////////
const char* toChar(const wchar_t* src)
{
	wchar_indexChar = (++wchar_indexChar) % WCHAR_MAX_COUNT;    // circulate index

	wcstombs(wchar_str[wchar_indexChar], src, WCHAR_MAX_LENGTH);// copy string as char
	wchar_str[wchar_indexChar][WCHAR_MAX_LENGTH - 1] = '\0';      // in case when source exceeded max length

	return wchar_str[wchar_indexChar];                          // return string as char
}

///////////////////////////////////////////////////////////////////////////////
// open a dialog to open obj file
///////////////////////////////////////////////////////////////////////////////
void ControllerForm::openFileDialog()
{
    const int MAX_CHARS = 512;

    OPENFILENAME ofn;
    memset(&ofn, 0, sizeof(ofn));

    wchar_t fileName[MAX_CHARS] = L"";  // file path to open

    ofn.lStructSize = sizeof(ofn);
    ofn.hwndOwner = 0;
    ofn.lpstrFilter = L"TXT files (*.txt)\0*.txt\0All Files (*.*)\0*.*\0\0";
    ofn.lpstrFile = fileName;
	ofn.lpstrFile[0] = '\0';
    ofn.nMaxFile = MAX_CHARS;
    ofn.Flags = OFN_EXPLORER | OFN_FILEMUSTEXIST | OFN_HIDEREADONLY;
    ofn.lpstrDefExt = L"txt";

    if(::GetOpenFileName(&ofn))
    {
        // change cursor
        HCURSOR oldCursor = ::GetCursor();
        HCURSOR cursor = ::LoadCursor(0, IDC_WAIT);
        ::SetCursor(cursor);

        model->initParametrsFromFile(toChar(fileName));

        // restore cursor
        ::SetCursor(oldCursor);
    }
}

///////////////////////////////////////////////////////////////////////////////
// open a dialog to open obj file
///////////////////////////////////////////////////////////////////////////////
void ControllerForm::saveFileDialog()
{
	const int MAX_CHARS = 512;

	OPENFILENAME ofn;
	memset(&ofn, 0, sizeof(ofn));

	wchar_t fileName[MAX_CHARS] = L"";  // file path to open

	ofn.lStructSize = sizeof(ofn);
	ofn.hwndOwner = 0;
	ofn.lpstrFilter = L"TXT files (*.txt)\0*.txt\0All Files (*.*)\0*.*\0\0";
	ofn.lpstrFile = fileName;
	ofn.lpstrFile[0] = '\0';
	ofn.nMaxFile = MAX_CHARS;
	ofn.Flags = OFN_EXPLORER | OFN_FILEMUSTEXIST | OFN_HIDEREADONLY;
	ofn.lpstrDefExt = L"txt";

	if (::GetSaveFileName(&ofn))
	{
		// change cursor
		HCURSOR oldCursor = ::GetCursor();
		HCURSOR cursor = ::LoadCursor(0, IDC_WAIT);
		::SetCursor(cursor);

		model->save(toChar(fileName));

		// restore cursor
		::SetCursor(oldCursor);
	}
}

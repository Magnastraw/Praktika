///////////////////////////////////////////////////////////////////////////////
// ControllerGL.cpp
// ================
// Derived Controller class for OpenGL window
// It is the controller of OpenGL rendering window. It initializes DC and RC,
// when WM_CREATE called, then, start new thread for OpenGL rendering loop.
//
// When this class is constructed, it gets the pointers to model and view
// components.
//
//  AUTHOR: Song Ho Ahn (song.ahn@gamil.com)
// CREATED: 2006-07-09
// UPDATED: 2016-05-11
///////////////////////////////////////////////////////////////////////////////

#include <process.h>                                // for _beginthreadex()
#include "ControllerGL.h"
using namespace Win;



///////////////////////////////////////////////////////////////////////////////
// default contructor
///////////////////////////////////////////////////////////////////////////////
ControllerGL::ControllerGL(ModelGL* model, ViewGL* view) : modelGL(model), viewGL(view),
                                                           threadHandle(0), threadId(0),
                                                           loopFlag(false), resizeFlag(false),
                                                           clientWidth(0), clientHeight(0)
{

}



///////////////////////////////////////////////////////////////////////////////
// handle WM_CLOSE
///////////////////////////////////////////////////////////////////////////////
int ControllerGL::close()
{
    loopFlag = false;
    ::WaitForSingleObject(threadHandle, INFINITE);  // wait for rendering thread is terminated

    // close OpenGL Rendering context
    viewGL->closeContext(handle);

    ::DestroyWindow(handle);
    return 0;
}



///////////////////////////////////////////////////////////////////////////////
// handle WM_DESTROY
///////////////////////////////////////////////////////////////////////////////
int ControllerGL::destroy()
{
    ::PostQuitMessage(0);       // exit the message loop
    return 0;
}



///////////////////////////////////////////////////////////////////////////////
// handle WM_CREATE
///////////////////////////////////////////////////////////////////////////////
int ControllerGL::create()
{
    // create a OpenGL rendering context
    if(!viewGL->createContext(handle, 32, 24, 8))
    {
        //Win::log(L"[ERROR] Failed to create OpenGL rendering context from ControllerGL::create().");
        return -1;
    }

    // create a thread for OpenGL rendering
    // The params of _beginthreadex() are security, stackSize, functionPtr, argPtr, initFlag, threadId.
    threadHandle = (HANDLE)_beginthreadex(0, 0, (unsigned (__stdcall *)(void *))threadFunction, this, 0, &threadId);
    if(threadHandle)
    {
        loopFlag = true;
        //Win::log(L"Created a rendering thread for OpenGL.");
    }
    else
    {
        ;//Win::log(L"[ERROR] Failed to create rendering thread from ControllerGL::create().");
    }

    return 0;
}



///////////////////////////////////////////////////////////////////////////////
// handle WM_PAINT
///////////////////////////////////////////////////////////////////////////////
int ControllerGL::paint()
{
    return 0;
}

///////////////////////////////////////////////////////////////////////////////
// handle WM_COMMAND
///////////////////////////////////////////////////////////////////////////////
int ControllerGL::command(int id, int cmd, LPARAM msg)
{
    return 0;
}



///////////////////////////////////////////////////////////////////////////////
// route to worker thread
// The function prototype is:
// unsigned int (__stdcall *)(void *)
///////////////////////////////////////////////////////////////////////////////
void ControllerGL::threadFunction(void* param)
{
	
    ((ControllerGL*)param)->runThread();
}



///////////////////////////////////////////////////////////////////////////////
// rendering thread
// initialize OpenGL states and start rendering loop
///////////////////////////////////////////////////////////////////////////////
void ControllerGL::runThread()
{
    // set the current RC in this thread
    ::wglMakeCurrent(viewGL->getDC(), viewGL->getRC());

    // initialize OpenGL states
    modelGL->init();

	//while (modelGL->runFlag == false)
	//{
	//	Sleep(10);
	//}

	//modelGL->initParametrs();
    // cofigure projection matrix
    RECT rect;
    ::GetClientRect(handle, &rect);
    modelGL->resizeWindow(rect.right, rect.bottom);

    // rendering loop
    while(loopFlag)
    {
        ::Sleep(1);                    // yield to other processes or threads

        modelGL->draw();
        viewGL->swapBuffers();
    }

    // terminate rendering thread
    ::wglMakeCurrent(0, 0);             // unset RC
    ::CloseHandle(threadHandle);
}



///////////////////////////////////////////////////////////////////////////////
// handle Left mouse down
///////////////////////////////////////////////////////////////////////////////
int ControllerGL::lButtonDown(WPARAM state, int x, int y)
{
    // update mouse position
    modelGL->setMousePosition(x, y);

    if(state == MK_LBUTTON)
    {
        modelGL->setMouseLeft(true);
    }

	::SetFocus(handle);

    return 0;
}



///////////////////////////////////////////////////////////////////////////////
// handle Left mouse up
///////////////////////////////////////////////////////////////////////////////
int ControllerGL::lButtonUp(WPARAM state, int x, int y)
{
    // update mouse position
    modelGL->setMousePosition(x, y);

    modelGL->setMouseLeft(false);

    return 0;
}



///////////////////////////////////////////////////////////////////////////////
// handle reft mouse down
///////////////////////////////////////////////////////////////////////////////
int ControllerGL::rButtonDown(WPARAM state, int x, int y)
{
    // update mouse position
    modelGL->setMousePosition(x, y);

    if(state == MK_RBUTTON)
    {
        modelGL->setMouseRight(true);
    }

	::SetFocus(handle);
    return 0;
}



///////////////////////////////////////////////////////////////////////////////
// handle reft mouse up
///////////////////////////////////////////////////////////////////////////////
int ControllerGL::rButtonUp(WPARAM state, int x, int y)
{
    // update mouse position
    modelGL->setMousePosition(x, y);

    modelGL->setMouseRight(false);

    return 0;
}



///////////////////////////////////////////////////////////////////////////////
// handle WM_MOUSEMOVE
///////////////////////////////////////////////////////////////////////////////
int ControllerGL::mouseMove(WPARAM state, int x, int y)
{
    if(state == MK_LBUTTON)
    {
        modelGL->rotateCamera(x, y);
    }
    if(state == MK_RBUTTON)
    {
       /* modelGL->cameraForward();*/
    }

    return 0;
}



///////////////////////////////////////////////////////////////////////////////
// handle WM_KEYDOWN
///////////////////////////////////////////////////////////////////////////////
int ControllerGL::keyDown(int key, LPARAM lParam)
{
    if(key == VK_ESCAPE)
    {
        ::PostMessage(handle, WM_CLOSE, 0, 0);
    }
	if (key == VK_SPACE)
	{
		modelGL->addSloi();
	}
	if (key == VK_UP)
	{
		modelGL->goUp();
	}
	if (key == VK_DOWN)
	{
		modelGL->goDown();
	}
	if (key == VK_LEFT)
	{
		modelGL->goLeft();
	}
	if (key == VK_RIGHT)
	{
		modelGL->goRight();
	}
	if (key == VK_NUMPAD1)
	{
		modelGL->goToUpSloi();
	}
	if (key == VK_NUMPAD3)
	{
		modelGL->goToDownSloi();
	}
	if (key == VK_NUMPAD0)
	{
		modelGL->goToStartPosition();
	}
	if (key == VK_RETURN)
	{
		modelGL->setEnter(true);
	}
	if (key == VK_BACK)
	{
		modelGL->deleteSloi();
	}
	if (key == 0x57)
	{
		modelGL->cameraForward();
	}
	if (key == 0x53)
	{
		modelGL->cameraBackward();
	}
	if (key == 0x41)
	{
		modelGL->cameraStrafeLeft();
	}
	if (key == 0x44)
	{
		modelGL->cameraStrafeRight();
	}
    return 0;
}

int ControllerGL::keyUp(int key, LPARAM lParam)
{
	if (key == VK_RETURN)
	{
		modelGL->setEnter(false);
		modelGL->Rect();
	}

	return 0;
}


///////////////////////////////////////////////////////////////////////////////
// handle WM_SIZE notification
// Note that the input param, width and height is for client area only.
// It excludes non-client area.
///////////////////////////////////////////////////////////////////////////////
int ControllerGL::size(int w, int h, WPARAM wParam)
{
	modelGL->resizeWindow(w, h);
	return 0;
}
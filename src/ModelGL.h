///////////////////////////////////////////////////////////////////////////////
// ModelGL.h
// =========
// Model component of OpenGL
// 
//  AUTHOR: Song Ho Ahn (song.ahn@gmail.com)
// CREATED: 2006-07-10
// UPDATED: 2006-07-10
///////////////////////////////////////////////////////////////////////////////

#include <string>

#ifndef MODEL_GL_H
#define MODEL_GL_H

class ModelGL
{
public:
	ModelGL();                                      // ctor
	~ModelGL();                                     // dtor

	void init();                                    // initialize OpenGL states
	void initParametrs();
	void initParametrsFromFile(const std::string fileName);
	void setCamera(float posX, float posY, float posZ, float targetX, float targetY, float targetZ);
	void setViewport(int width, int height);
	void draw();
	void resizeWindow(int width, int height);

	void setMouseLeft(bool flag) { mouseLeftDown = flag; };
	void setMouseRight(bool flag) { mouseRightDown = flag; };
	void setMousePosition(int x, int y) { mouseX = x; mouseY = y; };
	void setEnter(bool flag) { enterDown = flag; };
	void setE(float);
	void setMu(float);
	void setdX(double);
	void setdY(double);
	void setdZ(double);
	void setMoveX(bool);
	void setMoveY(bool);
	void setMoveZ(bool);
	void setP(double);
	void setP_bool(bool);
	void Rect();

	void rotateCamera(int x, int y);
	void zoomCamera(int dist);
	void cameraForward();
	void cameraBackward();
	void cameraStrafeRight();
	void cameraStrafeLeft();
	void addSloi();
	void goUp();
	void goDown();
	void goLeft();
	void goRight();
	void goToUpSloi();
	void goToDownSloi();
	void goToStartPosition();
	void updateTransform();
	void deleteSloi();
	void createObject(int,int,int, double, double, double, int, int, int);
	void save(const std::string fileName);

	wchar_t ep[8];
	wchar_t np[8];
	wchar_t color1[8];
	wchar_t color2[8];
	wchar_t color3[8];
	bool runFlag=false;
	bool buttonClicked = false;
protected:

private:
	// member functions
	void initLights();                              // add a white light ti scene

													// members
	bool mouseLeftDown;
	bool mouseRightDown;
	bool enterDown;
	int mouseX;
	int mouseY;
	float cameraAngleX;
	float cameraAngleY;
	float cameraDistance;
	float x_pos=0.0f;
	float z_pos=0.0f;
	int windowWidth;
	int windowHeight;
	bool windowResized;

};
#endif





#if _WIN32
#include <windows.h>    // include windows.h to avoid thousands of compile errors even though this class is not depending on Windows
#endif
#include <math.h>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <string>
#include "ModelGL.h"
#include <GL/gl.h>
#include <GL/glu.h>
#include <stdio.h>
#include <vector>

#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp>

using namespace std;

bool flag_enter;

double X0, Y0, Z0;
double A, B, C;

int NX, NY, NZ;
double hx, hy, hz;

int cNX, cNY, cNZ;
double a, b, c;

double dX, dY, dZ;
bool moveble_X, moveble_Y, moveble_Z;

GLdouble xxx;
GLdouble yyy;
GLdouble zzz;

float E, Mu;

int index_tetra;
int sloi;

double P=0; 

bool flag_enter_P = false;

// Initial position : on +Z
glm::vec3 position = glm::vec3(3, 5, 30);
glm::vec3 direction = glm::vec3(0, 0, 0);
glm::vec3 right_c = glm::vec3(0, 0, 0);

struct Cube {
	double coordinate[3];
	double dimension[3];
	GLfloat colour[3] = { 0.0f,0.0f,1.0f };
	double E;
	double mu;
	bool currentCub;
};

struct Node {
	double coordinate[3];
	double displacement[3];
	bool boundary[3];
	GLfloat colour[3] = { 0.0f,0.0f,1.0f };

	bool currentNode = false;
	bool flag;
	bool currentNodeRect = false;
};

struct Element {
	int coordinate[4];
	double centr[3];
	double colour[3] = { 0.0f,0.0f,1.0f };
	double volume;
	double E;
	double mu;
};

struct Triangle
{
	int coordinate[3];
	double normal[3];
	double normal0[3];
};

vector<Node> N_rect;

Triangle *tr;
Cube *cub;				//кубики для графики
Element *e;			    //элементы - тетраэдры
Node *n;				//узлы

int iteration = 0;

void print_M(double **MATRIX, int Ni, int Nj, string name) {
	ofstream F;
	F.open(name);

	for (int i = 0; i < Ni; i++) {
		F << i << ":	";
		for (int j = 0; j < Nj; j++) {
			F << setw(10) << MATRIX[i][j] << "	";
		}
		F << "	:" << i << endl;
	}

	F.close();
}

double func_b(Node *n, int m[7], int ii) {
	return (n[m[ii + 2]].coordinate[1] * n[m[ii + 3]].coordinate[2] +
		n[m[ii + 1]].coordinate[1] * n[m[ii + 2]].coordinate[2] +
		n[m[ii + 1]].coordinate[2] * n[m[ii + 3]].coordinate[1] -
		n[m[ii + 1]].coordinate[2] * n[m[ii + 2]].coordinate[1] -
		n[m[ii + 2]].coordinate[2] * n[m[ii + 3]].coordinate[1] -
		n[m[ii + 1]].coordinate[1] * n[m[ii + 3]].coordinate[2]) * pow(-1.0, ii + 1);
}

double func_c(Node *n, int m[7], int ii) {
	return (n[m[ii + 2]].coordinate[0] * n[m[ii + 3]].coordinate[2] +
		n[m[ii + 1]].coordinate[0] * n[m[ii + 2]].coordinate[2] +
		n[m[ii + 1]].coordinate[2] * n[m[ii + 3]].coordinate[0] -
		n[m[ii + 1]].coordinate[2] * n[m[ii + 2]].coordinate[0] -
		n[m[ii + 2]].coordinate[2] * n[m[ii + 3]].coordinate[0] -
		n[m[ii + 1]].coordinate[0] * n[m[ii + 3]].coordinate[2]) * pow(-1.0, ii);
}

double func_d(Node *n, int m[7], int ii) {
	return (n[m[ii + 2]].coordinate[0] * n[m[ii + 3]].coordinate[1] +
		n[m[ii + 1]].coordinate[0] * n[m[ii + 2]].coordinate[1] +
		n[m[ii + 1]].coordinate[1] * n[m[ii + 3]].coordinate[0] -
		n[m[ii + 1]].coordinate[1] * n[m[ii + 2]].coordinate[0] -
		n[m[ii + 2]].coordinate[1] * n[m[ii + 3]].coordinate[0] -
		n[m[ii + 1]].coordinate[0] * n[m[ii + 3]].coordinate[1]) * pow(-1.0, ii + 1);
}

void matrix_D(Element e, double D[6][6]) {
	double k0 = e.E * (1 - e.mu) / ((1 + e.mu) * (1 - 2 * e.mu));
	double k1 = e.mu / (1 - e.mu);
	double k2 = (1 - 2 * e.mu) / (2 * (1 - e.mu));
	double d[6][6] = {
		{ k0,		k0 * k1,	k0 * k1, 0,			0,			0 },
		{ k0 * k1,	k0,			k0 * k1, 0,			0,			0 },
		{ k0 * k1,	k0 * k1,	k0,		 0,			0,			0 },
		{ 0,        0,          0,		 k0 * k2,	0,			0 },
		{ 0,        0,          0,		 0,			k0 * k2,	0 },
		{ 0,        0,          0,		 0,			0,			k0 * k2 },
	};

	for (int i = 0; i < 6; i++) {
		for (int j = 0; j < 6; j++) {
			D[i][j] = d[i][j];
		}
	};
}

void matrix_B(double B[6][12], Node *n, int m[7]) {
	double b[4];
	double c[4];
	double d[4];

	for (int ii = 0; ii < 4; ii++) {
		b[ii] = func_b(n, m, ii);
		c[ii] = func_c(n, m, ii);
		d[ii] = func_d(n, m, ii);
	}

	double bb[6][12] = {
		{ b[0],		0,			0,			b[1],		0,			0,			b[2],		0,			0,			b[3],		0,			0 },
		{ 0,		c[0],		0,			0,			c[1],		0,			0,			c[2],		0,			0,			c[3],		0 },
		{ 0,		0,			d[0],		0,			0,			d[1],		0,			0,			d[2],		0,			0,			d[3] },
		{ c[0],		b[0],		0,			c[1],		b[1],		0,			c[2],		b[2],		0,			c[3],		b[3],		0 },
		{ d[0],		0,			b[0],		d[1],		0,			b[1],		d[2],		0,			b[2],		d[3],		0,			b[3] },
		{ 0,		d[0],		c[0],		0,			d[1],		c[1],		0,			d[2],		c[2],		0,			d[3],		c[3] },
	};

	for (int i = 0; i < 6; i++) {
		for (int j = 0; j < 12; j++) {
			B[i][j] = bb[i][j];
		}
	};
}

void Boundary(Node *n, double **K, double *F, int N) {
	for (int i = 0; i < N; i++) {    //учет граничных условий
		for (int k = 0; k < 3; k++) {
			if (n[i].boundary[k]) {
				double count = K[3 * i + k][3 * i + k];
				for (int j = 0; j < 3 * N; j++) {
					K[3 * i + k][j] = 0;
					F[j] -= n[i].displacement[k] * K[j][3 * i + k];
					K[j][3 * i + k] = 0;
				}
				K[3 * i + k][3 * i + k] = count;
				F[3 * i + k] = n[i].displacement[k] * count;
			}
		}
	}
}

void GaussianMethod(double **K, double *F, Node *n, int N) {
	double *U = new double[N];
	for (int s = 0; s < N; s++) {
		U[s] = F[s];
	}

	//приведение матрицы к треугольному виду
	double m;
	for (int k = 0; k < N; k++) {
		for (int j = k + 1; j < N; j++) {
			m = K[j][k] / K[k][k];
			for (int i = 0; i < N; i++) {
				K[j][i] -= m * K[k][i];
			}
			U[j] -= m * U[k];
		}
	}

	for (int i = 0; i < N; i++) {
		F[i] = U[i];
	}

	for (int i = N - 1; i >= 0; i--) {
		for (int j = i + 1; j < N; j++) {
			U[i] -= i != j ? K[i][j] * U[j] : 0;
		}
		U[i] = U[i] / K[i][i];
	}

	for (int i = 0; i < N / 3; i++)
	{
		for (int j = 0; j < 3; j++)
		{
			n[i].displacement[j] = U[3 * i + j];
		}
	}
}

void generate_Cube(Cube *cub, int cNX, int cNY, int cNZ, double X0, double Y0, double Z0, double a, double b, double c) {
	for (int i = 0; i < cNX; i++) {
		for (int j = 0; j < cNY; j++) {
			for (int k = 0; k < cNZ; k++) {
				int index = i * cNY * cNZ + j * cNZ + k;

				cub[index].coordinate[0] = X0 + i * a;
				cub[index].coordinate[1] = Y0 + j * b;
				cub[index].coordinate[2] = Z0 + k * c;

				cub[index].dimension[0] = a;
				cub[index].dimension[1] = b;
				cub[index].dimension[2] = c;

				cub[index].colour[0] = 0;
				cub[index].colour[1] = 0;
				cub[index].colour[2] = 1.0f;

				cub[index].E = 10;
				cub[index].mu = 0.49;

				cub[index].currentCub = false;
			}
		}
	}
}

void ijk(Triangle *tr, int NX, int NY, int NZ)
{
	int ind = 0;
	int sd = 0;
	int sd1 = 0;

	sd = NX * (NY + 1) * (NZ + 1);
	for (int iy = 0; iy < NY; iy++)
	{
		for (int iz = 0; iz < NZ; iz++)
		{
			tr[4 * ind].coordinate[0] = iy * (NZ + 1) + iz + (NZ + 1);
			tr[4 * ind].coordinate[1] = iy * (NZ + 1) + iz;
			tr[4 * ind].coordinate[2] = iy * (NZ + 1) + iz + (NZ + 1) + 1;

			tr[4 * ind + 1].coordinate[1] = iy * (NZ + 1) + iz + (NZ + 1) + 1;
			tr[4 * ind + 1].coordinate[0] = iy * (NZ + 1) + iz + 1;
			tr[4 * ind + 1].coordinate[2] = iy * (NZ + 1) + iz;

			tr[4 * ind + 2].coordinate[0] = tr[4 * ind].coordinate[1] + sd;
			tr[4 * ind + 2].coordinate[1] = tr[4 * ind].coordinate[0] + sd;
			tr[4 * ind + 2].coordinate[2] = tr[4 * ind].coordinate[2] + sd;

			tr[4 * ind + 3].coordinate[0] = tr[4 * ind + 1].coordinate[0] + sd;
			tr[4 * ind + 3].coordinate[1] = tr[4 * ind + 1].coordinate[2] + sd;
			tr[4 * ind + 3].coordinate[2] = tr[4 * ind + 1].coordinate[1] + sd;

			ind++;
		}
	}
	//+++

	sd1 = 4 * NY * NZ;
	sd = NY * (NZ + 1);
	ind = 0;
	for (int ix = 0; ix < NX; ix++)
	{
		for (int iz = 0; iz < NZ; iz++)
		{
			tr[sd1 + 4 * ind].coordinate[1] = ix * (NY + 1) * (NZ + 1) + iz;
			tr[sd1 + 4 * ind].coordinate[0] = ix * (NY + 1) * (NZ + 1) + iz + (NY + 1) * (NZ + 1);
			tr[sd1 + 4 * ind].coordinate[2] = ix * (NY + 1) * (NZ + 1) + iz + (NY + 1) * (NZ + 1) + 1;

			tr[sd1 + 4 * ind + 1].coordinate[0] = ix * (NY + 1) * (NZ + 1) + iz;
			tr[sd1 + 4 * ind + 1].coordinate[2] = ix * (NY + 1) * (NZ + 1) + iz + (NY + 1) * (NZ + 1) + 1;
			tr[sd1 + 4 * ind + 1].coordinate[1] = ix * (NY + 1) * (NZ + 1) + iz + 1;

			tr[sd1 + 4 * ind + 2].coordinate[0] = tr[sd1 + 4 * ind].coordinate[1] + sd;
			tr[sd1 + 4 * ind + 2].coordinate[1] = tr[sd1 + 4 * ind].coordinate[0] + sd;
			tr[sd1 + 4 * ind + 2].coordinate[2] = tr[sd1 + 4 * ind].coordinate[2] + sd;

			tr[sd1 + 4 * ind + 3].coordinate[0] = tr[sd1 + 4 * ind + 1].coordinate[0] + sd;
			tr[sd1 + 4 * ind + 3].coordinate[1] = tr[sd1 + 4 * ind + 1].coordinate[2] + sd;
			tr[sd1 + 4 * ind + 3].coordinate[2] = tr[sd1 + 4 * ind + 1].coordinate[1] + sd;

			ind++;
		}
	}

	sd1 += 4 * NX * NZ;
	sd = NZ;
	ind = 0;
	for (int ix = 0; ix < NX; ix++)
	{
		for (int iy = 0; iy < NY; iy++)
		{
			tr[sd1 + 4 * ind].coordinate[1] = ix * (NY + 1) * (NZ + 1) + iy * (NZ + 1);
			tr[sd1 + 4 * ind].coordinate[0] = ix * (NY + 1) * (NZ + 1) + iy * (NZ + 1) + (NY + 1) * (NZ + 1);
			tr[sd1 + 4 * ind].coordinate[2] = ix * (NY + 1) * (NZ + 1) + iy * (NZ + 1) + (NZ + 1);

			tr[sd1 + 4 * ind + 1].coordinate[0] = ix * (NY + 1) * (NZ + 1) + iy * (NZ + 1) + (NY + 1) * (NZ + 1);
			tr[sd1 + 4 * ind + 1].coordinate[2] = ix * (NY + 1) * (NZ + 1) + iy * (NZ + 1) + (NY + 1) * (NZ + 1) + (NZ + 1);
			tr[sd1 + 4 * ind + 1].coordinate[1] = ix * (NY + 1) * (NZ + 1) + iy * (NZ + 1) + (NZ + 1);

			tr[sd1 + 4 * ind + 2].coordinate[0] = tr[sd1 + 4 * ind].coordinate[1] + sd;
			tr[sd1 + 4 * ind + 2].coordinate[1] = tr[sd1 + 4 * ind].coordinate[0] + sd;
			tr[sd1 + 4 * ind + 2].coordinate[2] = tr[sd1 + 4 * ind].coordinate[2] + sd;

			tr[sd1 + 4 * ind + 3].coordinate[0] = tr[sd1 + 4 * ind + 1].coordinate[0] + sd;
			tr[sd1 + 4 * ind + 3].coordinate[1] = tr[sd1 + 4 * ind + 1].coordinate[2] + sd;
			tr[sd1 + 4 * ind + 3].coordinate[2] = tr[sd1 + 4 * ind + 1].coordinate[1] + sd;

			ind++;
		}
	}
}

void ijkp(Element *e, int NX, int NY, int NZ) {
	int i = 0;
	int index = 0;

	for (int ix = 0; ix < NX; ix++) {
		for (int iy = 0; iy < NY; iy++) {
			for (int iz = 0; iz < NZ; iz++) {
				e[index + 0].coordinate[0] = i + (NZ + 1) * (NY + 1) + 1;
				e[index + 0].coordinate[1] = i + NZ + 1 + 1;
				e[index + 0].coordinate[2] = i + (NZ + 1) * (NY + 2) + 1;
				e[index + 0].coordinate[3] = i + (NZ + 1) * (NY + 1);

				e[index + 1].coordinate[0] = i + NZ + 1;
				e[index + 1].coordinate[1] = i + (NZ + 1) * (NY + 2);
				e[index + 1].coordinate[2] = i + (NZ + 1) * (NY + 2) + 1;
				e[index + 1].coordinate[3] = i + (NZ + 1) * (NY + 1);

				e[index + 2].coordinate[0] = i + NZ + 1;
				e[index + 2].coordinate[1] = i + (NZ + 1) * (NY + 2) + 1;
				e[index + 2].coordinate[2] = i + NZ + 1 + 1;
				e[index + 2].coordinate[3] = i + (NZ + 1) * (NY + 1);

				e[index + 3].coordinate[0] = i + 1;
				e[index + 3].coordinate[1] = i + NZ + 1 + 1;
				e[index + 3].coordinate[2] = i + (NZ + 1) * (NY + 1) + 1;
				e[index + 3].coordinate[3] = i;

				e[index + 4].coordinate[0] = i + (NZ + 1) * (NY + 1);
				e[index + 4].coordinate[1] = i + (NZ + 1) * (NY + 1) + 1;
				e[index + 4].coordinate[2] = i + NZ + 1 + 1;
				e[index + 4].coordinate[3] = i;

				e[index + 5].coordinate[0] = i + (NZ + 1) * (NY + 1);
				e[index + 5].coordinate[1] = i + NZ + 1 + 1;
				e[index + 5].coordinate[2] = i + NZ + 1;
				e[index + 5].coordinate[3] = i;

				i++;
				index += 6;
			}
			i++;
		}
		i += NZ + 1;
	}

	index_tetra = index;
}

void xyz(Node *n, int NX, int NY, int NZ, double X0, double Y0, double Z0, double hx, double hy, double hz) {
	for (int i = 0; i <= NX; i++) {
		for (int j = 0; j <= NY; j++) {
			for (int k = 0; k <= NZ; k++) {
				int index = i * (NY + 1) * (NZ + 1) + j * (NZ + 1) + k;
				n[index].coordinate[0] = X0 + i * hx;
				n[index].coordinate[1] = Y0 + j * hy;
				n[index].coordinate[2] = Z0 + k * hz;

				if ((i == 0) || (i == NX) || (j == 0) || (j == NY) || (k == 0) || (k == NZ)) {
					for (int t = 0; t < 3; t++) {
						n[index].displacement[t] = 0;
						n[index].boundary[t] = true; ///////////////////////////////////////////!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
					}
				}
				else {
					for (int t = 0; t < 3; t++) {
						n[index].boundary[t] = false;
					}
				}
			}
		}
	}
}

void xyz1(Node *n, int NX, int NY, int NZ) {
	for (int i = 0; i <= NX; i++) {
		for (int j = 0; j <= NY; j++) {
			for (int k = 0; k <= NZ; k++) {
				int index = i * (NY + 1) * (NZ + 1) + j * (NZ + 1) + k;

				if ((j == NY / 2) && (k == NZ)) {
					n[index].displacement[0] = 0;
					n[index].displacement[1] = 0;
					n[index].displacement[2] = -0.2;
					n[index].boundary[0] = false;
					n[index].boundary[1] = false;
					n[index].boundary[2] = true;
				}
				else {
					n[index].boundary[0] = false;
					n[index].boundary[1] = false;
					n[index].boundary[2] = false;
				}

				if ((j == 0) || (j == NY)) {
					n[index].displacement[0] = 0;
					n[index].displacement[1] = 0;
					n[index].displacement[2] = 0;
					n[index].boundary[0] = true;
					n[index].boundary[1] = true;
					n[index].boundary[2] = true;
				}
			}
		}
	}
}

void size(Element *e, Node *n, int NX, int NY, int NZ) {
	for (int i = 0; i < 6 * NX * NY * NZ; i++) {
		e[i].volume = ((n[e[i].coordinate[1]].coordinate[0] - n[e[i].coordinate[0]].coordinate[0]) * (n[e[i].coordinate[2]].coordinate[1] - n[e[i].coordinate[0]].coordinate[1]) * (n[e[i].coordinate[3]].coordinate[2] - n[e[i].coordinate[0]].coordinate[2]) +
			(n[e[i].coordinate[1]].coordinate[1] - n[e[i].coordinate[0]].coordinate[1]) * (n[e[i].coordinate[2]].coordinate[2] - n[e[i].coordinate[0]].coordinate[2]) * (n[e[i].coordinate[3]].coordinate[0] - n[e[i].coordinate[0]].coordinate[0]) +
			(n[e[i].coordinate[1]].coordinate[2] - n[e[i].coordinate[0]].coordinate[2]) * (n[e[i].coordinate[2]].coordinate[0] - n[e[i].coordinate[0]].coordinate[0]) * (n[e[i].coordinate[3]].coordinate[1] - n[e[i].coordinate[0]].coordinate[1]) -
			(n[e[i].coordinate[1]].coordinate[2] - n[e[i].coordinate[0]].coordinate[2]) * (n[e[i].coordinate[2]].coordinate[1] - n[e[i].coordinate[0]].coordinate[1]) * (n[e[i].coordinate[3]].coordinate[0] - n[e[i].coordinate[0]].coordinate[0]) -
			(n[e[i].coordinate[1]].coordinate[0] - n[e[i].coordinate[0]].coordinate[0]) * (n[e[i].coordinate[2]].coordinate[2] - n[e[i].coordinate[0]].coordinate[2]) * (n[e[i].coordinate[3]].coordinate[1] - n[e[i].coordinate[0]].coordinate[1]) -
			(n[e[i].coordinate[1]].coordinate[1] - n[e[i].coordinate[0]].coordinate[1]) * (n[e[i].coordinate[2]].coordinate[0] - n[e[i].coordinate[0]].coordinate[0]) * (n[e[i].coordinate[3]].coordinate[2] - n[e[i].coordinate[0]].coordinate[2])) / 6;
	}
}

void centr(Element *e, Node *n, int NX, int NY, int NZ) {
	for (int i = 0; i < 6 * NX * NY * NZ; i++) {
		e[i].centr[0] = (n[e[i].coordinate[0]].coordinate[0] + n[e[i].coordinate[1]].coordinate[0] + n[e[i].coordinate[2]].coordinate[0] + n[e[i].coordinate[3]].coordinate[0]) / 4;
		e[i].centr[1] = (n[e[i].coordinate[0]].coordinate[1] + n[e[i].coordinate[1]].coordinate[1] + n[e[i].coordinate[2]].coordinate[1] + n[e[i].coordinate[3]].coordinate[1]) / 4;
		e[i].centr[2] = (n[e[i].coordinate[0]].coordinate[2] + n[e[i].coordinate[1]].coordinate[2] + n[e[i].coordinate[2]].coordinate[2] + n[e[i].coordinate[3]].coordinate[2]) / 4;
	}
}

void characteristic(Element *e, Node *n, Cube *cub, int NX, int NY, int NZ, int cNX, int cNY, int cNZ) {
	for (int i = 0; i < 6 * NX * NY * NZ; i++) {
		for (int j = 0; j < cNX * cNY * cNZ; j++) {
			if ((e[i].centr[0] > cub[j].coordinate[0]) && (e[i].centr[0] <= (cub[j].coordinate[0] + cub[j].dimension[0])) &&
				(e[i].centr[1] > cub[j].coordinate[1]) && (e[i].centr[1] <= (cub[j].coordinate[1] + cub[j].dimension[1])) &&
				(e[i].centr[2] > cub[j].coordinate[2]) && (e[i].centr[2] <= (cub[j].coordinate[2] + cub[j].dimension[2]))) {
				e[i].mu = cub[j].mu;
				e[i].E = cub[j].E;
				for (int k = 0; k < 3; k++) {
					e[i].colour[k] = cub[j].colour[k];
				}
			}
		}
	}
}

void Generate_K(Element *e, Node *n, double **K, double *F, int NX, int NY, int NZ) {
	for (int i = 0; i < 6 * NX * NY * NZ; i++) {
		int m[7] = { e[i].coordinate[0], e[i].coordinate[1], e[i].coordinate[2], e[i].coordinate[3], e[i].coordinate[0], e[i].coordinate[1], e[i].coordinate[2] };

		double D[6][6];
		matrix_D(e[i], D);

		double B[6][12];
		matrix_B(B, n, m);

		double Bt[12][6];
		for (int ii = 0; ii < 12; ii++) {
			for (int jj = 0; jj < 6; jj++) {
				Bt[ii][jj] = B[jj][ii];
			}
		}

		double c1[12][6];
		for (int ii = 0; ii < 12; ii++) {
			for (int jj = 0; jj < 6; jj++) {
				c1[ii][jj] = 0;
				for (int kk = 0; kk < 6; kk++) {
					c1[ii][jj] += Bt[ii][kk] * D[kk][jj];
				}
			}
		}
		double c2[12][12];
		for (int ii = 0; ii < 12; ii++) {
			for (int jj = 0; jj < 12; jj++) {
				c2[ii][jj] = 0;
				for (int kk = 0; kk < 6; kk++) {
					c2[ii][jj] += c1[ii][kk] * B[kk][jj];
				}
			}
		}

		for (int p = 0; p < 4; p++) {
			for (int q = 0; q < 4; q++) {
				for (int ii = 0; ii < 3; ii++) {
					for (int jj = 0; jj < 3; jj++) {
						K[3 * m[p] + ii][3 * m[q] + jj] += c2[3 * p + ii][3 * q + jj] / 36 / e[i].volume;
					}
				}
			}
		}
	}
}

void Generate_F(Triangle *tr,Node *n, double *F) {
	for (int i = 0; i < 4 * (NX * NY + NX * NZ + NY * NZ); i++) {
			int count = 0;
			for (int j = 0; j <(NX + 1) * (NY + 1) * (NZ + 1); j++) {
				if (n[j].flag) {
					if (tr[i].coordinate[0] == j) count++;
					if (tr[i].coordinate[1] == j) count++;
					if (tr[i].coordinate[2] == j) count++;
				}
			}

			if (count == 3)
			{
				//формирование правой части
				for (int k = 0; k < 3; k++)
				{
					for (int j = 0; j < 3; j++)
					{
						F[tr[i].coordinate[k]*3 + j] += P * tr[i].normal[j] / 6;
					}
			}
		}
	}
}

void Normal(Triangle *tr, Node *n, int NX, int NY, int NZ)
{
	for (int i = 0; i < 4 * (NX * NY + NX * NZ + NY * NZ); i++)
	{
		tr[i].normal[0] = (n[tr[i].coordinate[1]].coordinate[1] - n[tr[i].coordinate[0]].coordinate[1]) * (n[tr[i].coordinate[2]].coordinate[2] - n[tr[i].coordinate[0]].coordinate[2]) - (n[tr[i].coordinate[1]].coordinate[2] - n[tr[i].coordinate[0]].coordinate[2]) * (n[tr[i].coordinate[2]].coordinate[1] - n[tr[i].coordinate[0]].coordinate[1]);
		tr[i].normal[1] = (n[tr[i].coordinate[1]].coordinate[0] - n[tr[i].coordinate[0]].coordinate[0]) * (n[tr[i].coordinate[2]].coordinate[2] - n[tr[i].coordinate[0]].coordinate[2]) - (n[tr[i].coordinate[1]].coordinate[2] - n[tr[i].coordinate[0]].coordinate[2]) * (n[tr[i].coordinate[2]].coordinate[0] - n[tr[i].coordinate[0]].coordinate[0]);
		tr[i].normal[2] = (n[tr[i].coordinate[1]].coordinate[0] - n[tr[i].coordinate[0]].coordinate[0]) * (n[tr[i].coordinate[2]].coordinate[1] - n[tr[i].coordinate[0]].coordinate[1]) - (n[tr[i].coordinate[1]].coordinate[1] - n[tr[i].coordinate[0]].coordinate[1]) * (n[tr[i].coordinate[2]].coordinate[0] - n[tr[i].coordinate[0]].coordinate[0]);
		
		tr[i].normal0[0] = (n[tr[i].coordinate[0]].coordinate[0] + n[tr[i].coordinate[1]].coordinate[0]+ n[tr[i].coordinate[2]].coordinate[0])/3;
		tr[i].normal0[1] = (n[tr[i].coordinate[0]].coordinate[1] + n[tr[i].coordinate[1]].coordinate[1] + n[tr[i].coordinate[2]].coordinate[1]) / 3;
		tr[i].normal0[2] = (n[tr[i].coordinate[0]].coordinate[2] + n[tr[i].coordinate[1]].coordinate[2] + n[tr[i].coordinate[2]].coordinate[2]) / 3;
	}

}

ModelGL::ModelGL() : mouseLeftDown(false), mouseRightDown(false), enterDown(false),
cameraAngleX(-3), cameraAngleY(-0.5), cameraDistance(20)
{
}

ModelGL::~ModelGL()
{
}

void ModelGL::init() {
	glShadeModel(GL_SMOOTH);                        // shading mathod: GL_SMOOTH or GL_FLAT
	glPixelStorei(GL_UNPACK_ALIGNMENT, 4);          // 4-byte pixel alignment

													//enable /disable features
	glHint(GL_PERSPECTIVE_CORRECTION_HINT, GL_NICEST);
	//glHint(GL_LINE_SMOOTH_HINT, GL_NICEST);
	//glHint(GL_POLYGON_SMOOTH_HINT, GL_NICEST);
	glEnable(GL_DEPTH_TEST);
	//glEnable(GL_LIGHTING);
	//glEnable(GL_TEXTURE_2D);
	//glEnable(GL_CULL_FACE);
	//glEnable(GL_BLEND);

	//track material ambient and diffuse from surface color, call it before glEnable(GL_COLOR_MATERIAL)
	glColorMaterial(GL_FRONT_AND_BACK, GL_AMBIENT_AND_DIFFUSE);
	glEnable(GL_COLOR_MATERIAL);

	glClearColor(0, 0, 0, 0);                       // background color
	glClearStencil(0);                              // clear stencil buffer
	glClearDepth(1.0f);                             // 0 is near, 1 is far
	glDepthFunc(GL_LEQUAL);

	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();
}

void ModelGL::initParametrs() {
	delete[] cub;
	delete[] e;
	delete[] n;
	delete[] tr;
	flag_enter = false;

	//нижний левый угол
	X0 = 0;
	Y0 = 0;
	Z0 = 0;

	//шаги по осям
	hx = A / (double)NX;
	hy = B / NY;
	hz = C / NZ;
	hx;
	//размеры одного кубика
	a = A / cNX;
	b = B / cNY;
	c = C / cNZ;

	xxx = a / 2;
	yyy = b / 2;
	zzz = c / 2;

	sloi = 1;

	cub = new Cube[cNX * cNY * cNZ];									//кубики для графики
	e = new Element[6 * NX * NY * NZ];								//элементы - тетраэдры
	n = new Node[(NX + 1) * (NY + 1) * (NZ + 1)];						//узлы
	tr = new Triangle[4 * (NX * NY + NX * NZ + NY * NZ)];

	generate_Cube(cub, cNX, cNY, cNZ, X0, Y0, Z0, a, b, c);					//формирование сетки для графики	
	ijkp(e, NX, NY, NZ);													//формирование таблицы смежности
	ijk(tr, NX, NY, NZ);
	xyz(n, NX, NY, NZ, X0, Y0, Z0, hx, hy, hz);								//формирование сетки и свойств узлов
	Normal(tr, n, NX, NY, NZ);
	centr(e, n, NX, NY, NZ);
	characteristic(e, n, cub, NX, NY, NZ, cNX, cNY, cNZ);						//заполнение свойств тетраэдров, согласно расположению
}

void ModelGL::initParametrsFromFile(const std::string name) {
	delete[] cub;
	delete[] e;
	delete[] n;
	delete[] tr;

	string buff;

	flag_enter = false;


	ifstream fin;
	fin.open(name); // (ВВЕЛИ НЕ КОРРЕКТНОЕ ИМЯ ФАЙЛА)

	fin >> NX >> NY >> NZ >> cNX >> cNY >> cNZ;
	fin >> A >> B >> C;
	fin >> index_tetra;

	//нижний левый угол
	X0 = 0;
	Y0 = 0;
	Z0 = 0;

	//шаги по осям
	hx = A / NX;
	hy = B / NY;
	hz = C / NZ;

	//размеры одного кубика
	a = A / cNX;
	b = B / cNY;
	c = C / cNZ;

	xxx = a / 2;
	yyy = b / 2;
	zzz = c / 2;

	sloi = 1;

	cub = new Cube[cNX * cNY * cNZ];									//кубики для графики
	e = new Element[6 * NX * NY * NZ];								//элементы - тетраэдры
	n = new Node[(NX + 1) * (NY + 1) * (NZ + 1)];						//узлы
	tr = new Triangle[4 * (NX * NY + NX * NZ + NY * NZ)];


	for (int i = 0; i < cNX * cNY * cNZ; i++) {
		fin >> buff >> cub[i].coordinate[0] >> cub[i].coordinate[1] >> cub[i].coordinate[2];
		fin >> buff >> cub[i].dimension[0] >> cub[i].dimension[1] >> cub[i].dimension[2];
		fin >> buff >> cub[i].colour[0] >> cub[i].colour[1] >> cub[i].colour[2];
		fin >> buff >> cub[i].E >> cub[i].mu;
		fin >> buff >> cub[i].currentCub;
	}

	for (int i = 0; i < (NX + 1) * (NY + 1) * (NZ + 1); i++) {
		fin >> buff >> n[i].coordinate[0] >> n[i].coordinate[1] >> n[i].coordinate[2];
		fin >> buff >> n[i].displacement[0] >> n[i].displacement[1] >> n[i].displacement[2];
		fin >> buff >> n[i].boundary[0] >> n[i].boundary[1] >> n[i].boundary[2];
		fin >> buff >> n[i].colour[0] >> n[i].colour[1] >> n[i].colour[2];
		fin >> buff >> n[i].currentNode;
		fin >> buff >> n[i].flag;
		fin >> buff >> n[i].currentNodeRect;
	}

	for (int i = 0; i < 6 * NX * NY * NZ; i++) {
		fin >> buff >> e[i].coordinate[0] >> e[i].coordinate[1] >> e[i].coordinate[2] >> e[i].coordinate[3];
		fin >> buff >> e[i].centr[0] >> e[i].centr[1] >> e[i].centr[2];
		fin >> buff >> e[i].colour[0] >> e[i].colour[1] >> e[i].colour[2];
		fin >> buff >> e[i].volume;
		fin >> buff >> e[i].E >> e[i].mu;
	}

	for (int i = 0; i < 4 * (NX * NY + NX * NZ + NY * NZ); i++) {
		fin >> buff >> tr[i].coordinate[0] >> tr[i].coordinate[1] >> tr[i].coordinate[2];
		fin >> buff >> tr[i].normal[0] >> tr[i].normal[1] >> tr[i].normal[2];
		fin >> buff >> tr[i].normal0[0] >> tr[i].normal0[1] >> tr[i].normal0[2];
	}

	fin.close(); // закрываем файл

	//characteristic(e, n, cub, NX, NY, NZ, cNX, cNY, cNZ);						//заполнение свойств тетраэдров, согласно расположению
}

void ModelGL::setCamera(float posX, float posY, float posZ, float targetX, float targetY, float targetZ) {

	direction = glm::vec3(
		cos(cameraAngleY) * sin(cameraAngleX),
		sin(cameraAngleY),
		cos(cameraAngleY) * cos(cameraAngleX)
	);

	// Right vector
	right_c = glm::vec3(
		sin(cameraAngleX- 3.14f / 2.0f),
		0,
		cos(cameraAngleX - 3.14f / 2.0f)
	);

	// Up vector
	glm::vec3 up = glm::cross(right_c, direction);


	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();
	gluLookAt(position.x,position.y,position.z, position.x+direction.x, position.y + direction.y, position.z + direction.z, up.x, up.y, up.z); // eye(x,y,z), focal(x,y,z), up(x,y,z)
}

void ModelGL::setViewport(int w, int h) {

	windowWidth = w;
	windowHeight = h;

	// set viewport to be the entire window
	glViewport(0, 0, (GLsizei)w, (GLsizei)h);

	// set perspective viewing frustum
	float aspectRatio = (float)w / h;
	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
	gluPerspective(50.0f, (float)(w) / h, 0.1f, 100.0f); // FOV, AspectRatio, NearClip, FarClip
	glMatrixMode(GL_MODELVIEW);	// switch to modelview matrix in order to set scene
}

void ModelGL::resizeWindow(int w, int h) {
	// assign the width/height of viewport
	windowWidth = w;
	windowHeight = h;
	windowResized = true;
}

GLvoid Line(GLfloat X0, GLfloat Y0, GLfloat Z0, GLfloat x1, GLfloat y1, GLfloat z1) {
	glLineWidth((GLfloat)2.0);
	glBegin(GL_LINE_STRIP);
	glColor4f(0.0f, 1.0f, 0.0f, 0.0f);
	glVertex3f(X0, Y0, Z0);
	glVertex3f(x1, y1, z1);
	glEnd();
}

GLvoid Cub() {
	for (int ii = 0; ii < cNY; ii++)
		for (int jj = 0; jj < cNZ; jj++)
			for (int kk = 0; kk < sloi; kk++) {
				int i = kk * (cNY) * (cNZ)+ii * (cNZ)+jj;

				double x[8] = { cub[i].coordinate[0], cub[i].coordinate[0] + cub[i].dimension[0], cub[i].coordinate[0] + cub[i].dimension[0], cub[i].coordinate[0],
					cub[i].coordinate[0], cub[i].coordinate[0] + cub[i].dimension[0], cub[i].coordinate[0] + cub[i].dimension[0], cub[i].coordinate[0] };

				double y[8] = { cub[i].coordinate[1], cub[i].coordinate[1], cub[i].coordinate[1] + cub[i].dimension[1], cub[i].coordinate[1] + cub[i].dimension[1],
					cub[i].coordinate[1], cub[i].coordinate[1], cub[i].coordinate[1] + cub[i].dimension[1], cub[i].coordinate[1] + cub[i].dimension[1] };

				double z[8] = { cub[i].coordinate[2], cub[i].coordinate[2], cub[i].coordinate[2], cub[i].coordinate[2],
					cub[i].coordinate[2] + cub[i].dimension[2], cub[i].coordinate[2] + cub[i].dimension[2], cub[i].coordinate[2] + cub[i].dimension[2], cub[i].coordinate[2] + cub[i].dimension[2] };

				glBegin(GL_QUADS);			// Рисуем куб		

				glEnable(GL_BLEND);
				glDisable(GL_DEPTH_TEST);

				glColor3f(cub[i].colour[0], cub[i].colour[1], cub[i].colour[2]);
				glVertex3f(x[0], y[0], z[0]);		// Право верх квадрата (Верх)
				glVertex3f(x[1], y[1], z[1]);		// Лево верх
				glVertex3f(x[2], y[2], z[2]);		// Лево низ
				glVertex3f(x[3], y[3], z[3]);		// Право низ*

				glColor3f(cub[i].colour[0], cub[i].colour[1], cub[i].colour[2]);
				glVertex3f(x[4], y[4], z[4]);		// Верх право квадрата (Низ)
				glVertex3f(x[5], y[5], z[5]);		// Верх лево
				glVertex3f(x[6], y[6], z[6]);		// Низ лево
				glVertex3f(x[7], y[7], z[7]);		// Низ право

				glColor3f(cub[i].colour[0], cub[i].colour[1], cub[i].colour[2]);
				glVertex3f(x[0], y[0], z[0]);		// Верх право квадрата (Перед)
				glVertex3f(x[3], y[3], z[3]);		// Верх лево
				glVertex3f(x[7], y[7], z[7]);		// Низ лево
				glVertex3f(x[4], y[4], z[4]);		// Низ право

				glColor3f(cub[i].colour[0], cub[i].colour[1], cub[i].colour[2]);
				glVertex3f(x[1], y[1], z[1]);		// Верх право квадрата (Зад)
				glVertex3f(x[2], y[2], z[2]);		// Верх лево
				glVertex3f(x[6], y[6], z[6]);		// Низ лево
				glVertex3f(x[5], y[5], z[5]);		// Низ право

				glColor3f(cub[i].colour[0], cub[i].colour[1], cub[i].colour[2]);
				glVertex3f(x[0], y[0], z[0]);		// Верх право квадрата (Лево)
				glVertex3f(x[1], y[1], z[1]);		// Верх лево
				glVertex3f(x[5], y[5], z[5]);		// Низ лево
				glVertex3f(x[4], y[4], z[4]);		// Низ право

				glColor3f(cub[i].colour[0], cub[i].colour[1], cub[i].colour[2]);
				glVertex3f(x[3], y[3], z[3]);		// Верх право квадрата (Право)
				glVertex3f(x[2], y[2], z[2]);		// Верх лево
				glVertex3f(x[6], y[6], z[6]);		// Низ лево
				glVertex3f(x[7], y[7], z[7]);		// Низ право.

				glEnd();

				Line(x[0], y[0], z[0], x[1], y[1], z[1]);
				Line(x[1], y[1], z[1], x[2], y[2], z[2]);
				Line(x[2], y[2], z[2], x[3], y[3], z[3]);
				Line(x[3], y[3], z[3], x[0], y[0], z[0]);
				Line(x[4], y[4], z[4], x[5], y[5], z[5]);
				Line(x[5], y[5], z[5], x[6], y[6], z[6]);
				Line(x[6], y[6], z[6], x[7], y[7], z[7]);
				Line(x[7], y[7], z[7], x[4], y[4], z[4]);
				Line(x[0], y[0], z[0], x[4], y[4], z[4]);
				Line(x[1], y[1], z[1], x[5], y[5], z[5]);
				Line(x[2], y[2], z[2], x[6], y[6], z[6]);
				Line(x[3], y[3], z[3], x[7], y[7], z[7]);
			}
}

GLvoid Cub_tetra() {
	for (int i = 0; i < index_tetra; i++) {
		glBegin(GL_POLYGON);
		glEnable(GL_BLEND);
		glDisable(GL_DEPTH_TEST);
		glColor3f(e[i].colour[0], e[i].colour[1], e[i].colour[2]);
		glVertex3f(n[e[i].coordinate[0]].coordinate[0], n[e[i].coordinate[0]].coordinate[1], n[e[i].coordinate[0]].coordinate[2]);
		glVertex3f(n[e[i].coordinate[1]].coordinate[0], n[e[i].coordinate[1]].coordinate[1], n[e[i].coordinate[1]].coordinate[2]);
		glVertex3f(n[e[i].coordinate[2]].coordinate[0], n[e[i].coordinate[2]].coordinate[1], n[e[i].coordinate[2]].coordinate[2]);
		glEnd();

		glBegin(GL_POLYGON);
		glEnable(GL_BLEND);
		glDisable(GL_DEPTH_TEST);
		glColor3f(e[i].colour[0], e[i].colour[1], e[i].colour[2]);
		glVertex3f(n[e[i].coordinate[1]].coordinate[0], n[e[i].coordinate[1]].coordinate[1], n[e[i].coordinate[1]].coordinate[2]);
		glVertex3f(n[e[i].coordinate[2]].coordinate[0], n[e[i].coordinate[2]].coordinate[1], n[e[i].coordinate[2]].coordinate[2]);
		glVertex3f(n[e[i].coordinate[3]].coordinate[0], n[e[i].coordinate[3]].coordinate[1], n[e[i].coordinate[3]].coordinate[2]);
		glEnd();

		glBegin(GL_POLYGON);
		glEnable(GL_BLEND);
		glDisable(GL_DEPTH_TEST);
		glColor3f(e[i].colour[0], e[i].colour[1], e[i].colour[2]);
		glVertex3f(n[e[i].coordinate[2]].coordinate[0], n[e[i].coordinate[2]].coordinate[1], n[e[i].coordinate[2]].coordinate[2]);
		glVertex3f(n[e[i].coordinate[3]].coordinate[0], n[e[i].coordinate[3]].coordinate[1], n[e[i].coordinate[3]].coordinate[2]);
		glVertex3f(n[e[i].coordinate[0]].coordinate[0], n[e[i].coordinate[0]].coordinate[1], n[e[i].coordinate[0]].coordinate[2]);
		glEnd();

		glBegin(GL_POLYGON);
		glEnable(GL_BLEND);
		glDisable(GL_DEPTH_TEST);
		glColor3f(e[i].colour[0], e[i].colour[1], e[i].colour[2]);
		glVertex3f(n[e[i].coordinate[3]].coordinate[0], n[e[i].coordinate[3]].coordinate[1], n[e[i].coordinate[3]].coordinate[2]);
		glVertex3f(n[e[i].coordinate[0]].coordinate[0], n[e[i].coordinate[0]].coordinate[1], n[e[i].coordinate[0]].coordinate[2]);
		glVertex3f(n[e[i].coordinate[1]].coordinate[0], n[e[i].coordinate[1]].coordinate[1], n[e[i].coordinate[1]].coordinate[2]);
		glEnd();


		Line(n[e[i].coordinate[0]].coordinate[0], n[e[i].coordinate[0]].coordinate[1], n[e[i].coordinate[0]].coordinate[2], n[e[i].coordinate[1]].coordinate[0], n[e[i].coordinate[1]].coordinate[1], n[e[i].coordinate[1]].coordinate[2]);
		Line(n[e[i].coordinate[1]].coordinate[0], n[e[i].coordinate[1]].coordinate[1], n[e[i].coordinate[1]].coordinate[2], n[e[i].coordinate[2]].coordinate[0], n[e[i].coordinate[2]].coordinate[1], n[e[i].coordinate[2]].coordinate[2]);
		Line(n[e[i].coordinate[0]].coordinate[0], n[e[i].coordinate[0]].coordinate[1], n[e[i].coordinate[0]].coordinate[2], n[e[i].coordinate[2]].coordinate[0], n[e[i].coordinate[2]].coordinate[1], n[e[i].coordinate[2]].coordinate[2]);
		Line(n[e[i].coordinate[1]].coordinate[0], n[e[i].coordinate[1]].coordinate[1], n[e[i].coordinate[1]].coordinate[2], n[e[i].coordinate[3]].coordinate[0], n[e[i].coordinate[3]].coordinate[1], n[e[i].coordinate[3]].coordinate[2]);
		Line(n[e[i].coordinate[3]].coordinate[0], n[e[i].coordinate[3]].coordinate[1], n[e[i].coordinate[3]].coordinate[2], n[e[i].coordinate[2]].coordinate[0], n[e[i].coordinate[2]].coordinate[1], n[e[i].coordinate[2]].coordinate[2]);
		Line(n[e[i].coordinate[0]].coordinate[0], n[e[i].coordinate[0]].coordinate[1], n[e[i].coordinate[0]].coordinate[2], n[e[i].coordinate[3]].coordinate[0], n[e[i].coordinate[3]].coordinate[1], n[e[i].coordinate[3]].coordinate[2]);
	}
}

GLvoid Triangle_draw()
{
	for (int i = 0; i < 4 * (NX * NY + NX * NZ + NY * NZ); i++) {

		glBegin(GL_POLYGON);
		glEnable(GL_BLEND);
		glDisable(GL_DEPTH_TEST);
		glColor3f(1.0f, 0.0f, 0.0f);
		glVertex3f(n[tr[i].coordinate[0]].coordinate[0], n[tr[i].coordinate[0]].coordinate[1], n[tr[i].coordinate[0]].coordinate[2]);
		glVertex3f(n[tr[i].coordinate[1]].coordinate[0], n[tr[i].coordinate[1]].coordinate[1], n[tr[i].coordinate[1]].coordinate[2]);
		glVertex3f(n[tr[i].coordinate[2]].coordinate[0], n[tr[i].coordinate[2]].coordinate[1], n[tr[i].coordinate[2]].coordinate[2]);
		glEnd();

		Line(n[tr[i].coordinate[0]].coordinate[0], n[tr[i].coordinate[0]].coordinate[1], n[tr[i].coordinate[0]].coordinate[2], n[tr[i].coordinate[1]].coordinate[0], n[tr[i].coordinate[1]].coordinate[1], n[tr[i].coordinate[1]].coordinate[2]);
		Line(n[tr[i].coordinate[1]].coordinate[0], n[tr[i].coordinate[1]].coordinate[1], n[tr[i].coordinate[1]].coordinate[2], n[tr[i].coordinate[2]].coordinate[0], n[tr[i].coordinate[2]].coordinate[1], n[tr[i].coordinate[2]].coordinate[2]);
		Line(n[tr[i].coordinate[0]].coordinate[0], n[tr[i].coordinate[0]].coordinate[1], n[tr[i].coordinate[0]].coordinate[2], n[tr[i].coordinate[2]].coordinate[0], n[tr[i].coordinate[2]].coordinate[1], n[tr[i].coordinate[2]].coordinate[2]);
	}
}

GLvoid Points(Node node_) {
	glPointSize(5);
	glBegin(GL_POINTS);
	glEnable(GL_BLEND);
	glDisable(GL_DEPTH_TEST);
	glColor3f(node_.colour[0], node_.colour[1], node_.colour[2]);
	glVertex3f(node_.coordinate[0], node_.coordinate[1], node_.coordinate[2]);
	glEnd();
}

void ModelGL::Rect()
{
		for (int i = 0; i < NY + 1; i++)
			for (int j = 0; j < NZ + 1; j++)
				for (int k = 0; k < NX + 1; k++) {
					if (n[k * (NY + 1) * (NZ + 1) + i * (NZ + 1) + j].currentNodeRect)
					{
						N_rect.push_back(n[k * (NY + 1) * (NZ + 1) + i * (NZ + 1) + j]);
					}
				}
	
	if (N_rect.size() == 3) {
	Node node_0 = N_rect.at(1);
	Node node_1 = N_rect.at(2);

	if (node_0.coordinate[0] > node_1.coordinate[0])
	{
		double node_temp = node_0.coordinate[0];
		node_0.coordinate[0] = node_0.coordinate[0] - abs(node_1.coordinate[0] - node_0.coordinate[0]);
		node_1.coordinate[0] = node_temp;
	}
	if (node_0.coordinate[1] > node_1.coordinate[1])
	{
		double node_temp = node_0.coordinate[1];
		node_0.coordinate[1] = node_0.coordinate[1] - abs(node_1.coordinate[1] - node_0.coordinate[1]);
		node_1.coordinate[1] = node_temp;
	}
	if (node_0.coordinate[2] > node_1.coordinate[2])
	{
		double node_temp = node_0.coordinate[2];
		node_0.coordinate[2] = node_0.coordinate[2] - abs(node_1.coordinate[2] - node_0.coordinate[2]);
		node_1.coordinate[2] = node_temp;
	}

	for (int i = 0; i < NY + 1; i++)
		for (int j = 0; j < NZ + 1; j++)
			for (int k = 0; k < NX + 1; k++) {
					n[k * (NY + 1) * (NZ + 1) + i * (NZ + 1) + j].currentNodeRect = false;
					if (((n[k * (NY + 1) * (NZ + 1) + i * (NZ + 1) + j].coordinate[0] >= node_0.coordinate[0]) && (n[k * (NY + 1) * (NZ + 1) + i * (NZ + 1) + j].coordinate[0] <= node_1.coordinate[0])) &&
						((n[k * (NY + 1) * (NZ + 1) + i * (NZ + 1) + j].coordinate[1] >= node_0.coordinate[1]) && (n[k * (NY + 1) * (NZ + 1) + i * (NZ + 1) + j].coordinate[1] <= node_1.coordinate[1])) &&
						((n[k * (NY + 1) * (NZ + 1) + i * (NZ + 1) + j].coordinate[2] >= node_0.coordinate[2]) && (n[k * (NY + 1) * (NZ + 1) + i * (NZ + 1) + j].coordinate[2] <= node_1.coordinate[2]))) {		
						if (flag_enter_P)
						{
							n[k * (NY + 1) * (NZ + 1) + i * (NZ + 1) + j].flag = true;
							n[k * (NY + 1) * (NZ + 1) + i * (NZ + 1) + j].currentNode = true;
							n[k * (NY + 1) * (NZ + 1) + i * (NZ + 1) + j].colour[0] = 1.0f;
							n[k * (NY + 1) * (NZ + 1) + i * (NZ + 1) + j].colour[1] = 0.0f;
							n[k * (NY + 1) * (NZ + 1) + i * (NZ + 1) + j].colour[2] = 0.0f;
						} else{
							n[k * (NY + 1) * (NZ + 1) + i * (NZ + 1) + j].currentNode = true;
							n[k * (NY + 1) * (NZ + 1) + i * (NZ + 1) + j].colour[0] = 1.0f;
							n[k * (NY + 1) * (NZ + 1) + i * (NZ + 1) + j].colour[1] = 1.0f;
							n[k * (NY + 1) * (NZ + 1) + i * (NZ + 1) + j].colour[2] = 0.8f;
						}
						n[k * (NY + 1) * (NZ + 1) + i * (NZ + 1) + j].displacement[0] = dX;
						n[k * (NY + 1) * (NZ + 1) + i * (NZ + 1) + j].displacement[1] = dY;
						n[k * (NY + 1) * (NZ + 1) + i * (NZ + 1) + j].displacement[2] = dZ;
						n[k * (NY + 1) * (NZ + 1) + i * (NZ + 1) + j].boundary[0] = moveble_X;
						n[k * (NY + 1) * (NZ + 1) + i * (NZ + 1) + j].boundary[1] = moveble_Y;
						n[k * (NY + 1) * (NZ + 1) + i * (NZ + 1) + j].boundary[2] = moveble_Z;
					}
				}
			N_rect.erase(N_rect.begin(),N_rect.end());
			}
}

void ModelGL::draw() {

	if (windowResized) {
		setViewport(windowWidth, windowHeight);
		windowResized = false;
	}

	setCamera(0, 0, 0, 0, 0, 0);


	// clear buffer
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT | GL_STENCIL_BUFFER_BIT);

	// save the initial ModelView matrix before modifying ModelView matrix
	glPushMatrix();

	// tramsform camera
	glTranslatef(0, 0, cameraDistance);
	//glRotatef(cameraAngleX, 1, 0, 0);   // pitch
	//glRotatef(cameraAngleY, 0, 1, 0);   // heading

	//X ось
	glLineWidth((GLfloat)2.0);
	glBegin(GL_LINE_STRIP);
	glColor4f(1.0f, 0.0f, 0.0f, 0.0f);
	glVertex3f(0.0, 0.0, 0.0);
	glVertex3f(100.0, 0.0, 0.0);
	glEnd();

	//Y ось
	glLineWidth((GLfloat)2.0);
	glBegin(GL_LINE_STRIP);
	glColor4f(1.0f, 1.0f, 0.0f, 0.0f);
	glVertex3f(0.0, 0.0, 0.0);
	glVertex3f(0.0, 100.0, 0.0);
	glEnd();

	//Z ось
	glLineWidth((GLfloat)2.0);
	glBegin(GL_LINE_STRIP);
	glColor4f(0.0f, 0.0f, 1.0f, 0.0f);
	glVertex3f(0.0, 0.0, 0.0);
	glVertex3f(0.0, 0.0, 100.0);
	glEnd();

	//for (int i = 0; i < 4 * (NX * NY + NX * NZ + NY * NZ); i++)
	//{
	//	Line(tr[i].normal0[0], tr[i].normal0[1], tr[i].normal0[2], tr[i].normal0[0] + tr[i].normal[0], tr[i].normal0[1]+ tr[i].normal[1], tr[i].normal0[2]+ tr[i].normal[2]);
	//}

	if (!flag_enter) {
		/*Triangle_draw();*/
		iteration = 0;
		for (int i = 0; i < cNY; i++)
			for (int j = 0; j < cNZ; j++)
				for (int k = 0; k < sloi; k++) {
					int ind = k * (cNY * cNZ) + i * cNZ + j;

					if ((xxx >= cub[ind].coordinate[0]) &&
						(xxx < cub[ind].coordinate[0] + cub[ind].dimension[0]) &&

						(yyy >= cub[ind].coordinate[1]) &&
						(yyy < cub[ind].coordinate[1] + cub[ind].dimension[1]) &&

						(zzz >= cub[ind].coordinate[2]) &&
						(zzz < cub[ind].coordinate[2] + cub[ind].dimension[2])) {
						if (!cub[ind].currentCub) {
							cub[ind].colour[0] = 1.0f;
							cub[ind].colour[1] = 1.0f;
							cub[ind].colour[2] = 1.0f;
						}
						if (enterDown) {
							cub[ind].currentCub = true;
							cub[ind].E = E;
							cub[ind].mu = Mu;
							cub[ind].colour[0] = 1.0f;
							cub[ind].colour[1] = 1.0f - Mu;
							cub[ind].colour[2] = 1.0f - 0.25*Mu;
						}
					}
					else if (!cub[ind].currentCub) {
						cub[ind].colour[0] = 0.0f;
						cub[ind].colour[1] = 0.0f;
						cub[ind].colour[2] = 1.0f;
					}
				}
		Cub();
	}
	else {
		Cub_tetra();
		/*Triangle_draw();*/
		for (int i = 0; i < NY + 1; i++)
			for (int j = 0; j < NZ + 1; j++)
				for (int k = 0; k < NX + 1; k++) {
					if (((xxx >= n[k * (NY + 1) * (NZ + 1)].coordinate[0]-0.1) && (xxx <= n[k * (NY + 1) * (NZ + 1)].coordinate[0] + 0.1)) &&
						((yyy >= n[k * (NY + 1) * (NZ + 1) + i * (NZ + 1) + j].coordinate[1]-0.1) && (yyy <= n[k * (NY + 1) * (NZ + 1) + i * (NZ + 1) + j].coordinate[1] + 0.1)) &&
						((zzz >= n[k * (NY + 1) * (NZ + 1) + i * (NZ + 1) + j].coordinate[2]-0.1) && (zzz <= n[k * (NY + 1) * (NZ + 1) + i * (NZ + 1) + j].coordinate[2] + 0.1))) {
						
							n[k * (NY + 1) * (NZ + 1) + i * (NZ + 1) + j].colour[0] = 1.0f;
							n[k * (NY + 1) * (NZ + 1) + i * (NZ + 1) + j].colour[1] = 1.0f;
							n[k * (NY + 1) * (NZ + 1) + i * (NZ + 1) + j].colour[2] = 0.0f;
						
						if (enterDown) {
							n[k * (NY + 1) * (NZ + 1) + i * (NZ + 1) + j].currentNode = true;
							n[k * (NY + 1) * (NZ + 1) + i * (NZ + 1) + j].currentNodeRect = true;
							n[k * (NY + 1) * (NZ + 1) + i * (NZ + 1) + j].colour[0] = 1.0f;
							n[k * (NY + 1) * (NZ + 1) + i * (NZ + 1) + j].colour[1] = 0.0f;
							n[k * (NY + 1) * (NZ + 1) + i * (NZ + 1) + j].colour[2] = 0.0f;
						}
					}
					else if (n[k * (NY + 1) * (NZ + 1) + i * (NZ + 1) + j].currentNode)
					{
						if (n[k * (NY + 1) * (NZ + 1) + i * (NZ + 1) + j].flag == true)
						{
							n[k * (NY + 1) * (NZ + 1) + i * (NZ + 1) + j].colour[0] = 0.5f;
							n[k * (NY + 1) * (NZ + 1) + i * (NZ + 1) + j].colour[1] = 0.35f;
							n[k * (NY + 1) * (NZ + 1) + i * (NZ + 1) + j].colour[2] = 0.05f;
						}
						else {
							n[k * (NY + 1) * (NZ + 1) + i * (NZ + 1) + j].colour[0] = 1.0f;
							n[k * (NY + 1) * (NZ + 1) + i * (NZ + 1) + j].colour[1] = 0.0f;
							n[k * (NY + 1) * (NZ + 1) + i * (NZ + 1) + j].colour[2] = 0.0f;
						}
					}
					else if (!n[k * (NY + 1) * (NZ + 1) + i * (NZ + 1) + j].currentNode)
					{
						n[k * (NY + 1) * (NZ + 1) + i * (NZ + 1) + j].colour[0] = 0.0f;
						n[k * (NY + 1) * (NZ + 1) + i * (NZ + 1) + j].colour[1] = 0.0f;
						n[k * (NY + 1) * (NZ + 1) + i * (NZ + 1) + j].colour[2] = 1.0f;
						n[k * (NY + 1) * (NZ + 1) + i * (NZ + 1) + j].flag = false;
					}
					
					Points(n[k * (NY + 1) * (NZ + 1) + i * (NZ + 1) + j]);
				}
	}
	glPopMatrix();
}

void ModelGL::rotateCamera(int x, int y) {
	if (mouseLeftDown) {
		cameraAngleX -= (x - mouseX)*0.005;
		cameraAngleY -= (y - mouseY)*0.005;
		mouseX = x;
		mouseY = y;
	}
}

void ModelGL::zoomCamera(int delta) {
	if (mouseRightDown) {
		cameraDistance += (delta - mouseY) * 0.05f;
		mouseY = delta;
	}
}

void ModelGL::addSloi() {
	if (sloi < cNX) {
		sloi++;
		ModelGL::goToStartPosition();
	}
}

void ModelGL::deleteSloi() {
	if (sloi > 1) {
		sloi--;
		ModelGL::goToStartPosition();
	}
}

void ModelGL::goUp() {
	if (!flag_enter) {
		if (yyy < (GLdouble)(cNY - 1) * b) {
			yyy += b;
		}
	}
	else {
		if (yyy < (GLdouble)(NY)* hy) {
			yyy += hy;
		}
	}
}

void ModelGL::goDown() {
	if (!flag_enter) {
		if (yyy >= b) {
			yyy -= b;
		}
	}
	else
	{
		if (yyy >= hy) {
			yyy -= hy;
		}
	}
}

void ModelGL::goRight() {
	if (!flag_enter) {
		if (zzz >= c) {
			zzz -= c;
		}
	}
	else {
		if (zzz >= hz) {
			zzz -= hz;
		}
	}
}

void ModelGL::goLeft() {
	if (!flag_enter) {
		if (zzz < (GLdouble)(cNZ - 1) * c) {
			zzz += c;
		}
	}
	else {
		if (zzz < (GLdouble)(NZ)* hz) {
			zzz += hz;
		}
	}
}

void ModelGL::goToUpSloi() {
	if (!flag_enter) {
		if ((xxx < (GLdouble)(sloi-1)* a) && (sloi > 1)) {
			xxx += a;
		}
	}
	else {
		if ((xxx < (GLdouble)(NX)* hx) && (NX > 1)) {
			xxx += hx;
		}
	}
}

void ModelGL::goToDownSloi() {
	if (!flag_enter) {
		if (xxx >= a) {
			xxx -= a;
		}
	}
	else {
		if (xxx >= hx) {
			xxx -= hx;
		}
	}
}

void ModelGL::goToStartPosition() {
	xxx = a * (sloi - 1) + a / 2;
	yyy = b / 2;
	zzz = c / 2;
}

void ModelGL::cameraForward() {
	position += direction * 0.1f;
}

void ModelGL::cameraBackward() {
	position -= direction * 0.1f;
}

void ModelGL::cameraStrafeRight() {
	position += right_c * 0.1f;
}

void ModelGL::cameraStrafeLeft() {
	position -= right_c * 0.1f;
}

void ModelGL::updateTransform()
{
	characteristic(e, n, cub, NX, NY, NZ, cNX, cNY, cNZ);					//заполнение свойств тетраэдров, согласно расположению
	xxx = 0;
	yyy = 0;
	zzz = 0;
	size(e, n, NX, NY, NZ);													//вычисление объемов, каждую итерацию(каждый клик кнопки), тела меняются	
	Normal(tr, n, NX, NY, NZ);
	flag_enter = true;
	if (iteration >= 1) {

		double *F = new double[3 * (NX + 1) * (NY + 1) * (NZ + 1)];				//вектор правой части, заполняется давлением
		double **K = new double *[3 * (NX + 1) * (NY + 1) * (NZ + 1)];			//глобальная матрица жесткости

		for (int i = 0; i < 3 * (NX + 1) * (NY + 1) * (NZ + 1); i++) {
			K[i] = new double[3 * (NX + 1) * (NY + 1) * (NZ + 1)];
		}
		for (int i = 0; i < 3 * (NX + 1) * (NY + 1) * (NZ + 1); i++) {
			for (int j = 0; j < 3 * (NX + 1) * (NY + 1) * (NZ + 1); j++) {
				K[i][j] = 0;
			}
			F[i] = 0;
		}

		Generate_K(e, n, K, F, NX, NY, NZ);										//формирование глобальной матрицы жесткости и вектора правой части
		Generate_F(tr, n, F);
		Boundary(n, K, F, (NX + 1) * (NY + 1) * (NZ + 1));						//учет граничных условий	
		GaussianMethod(K, F, n, 3 * (NX + 1) * (NY + 1) * (NZ + 1));			//поиск перемещений

		for (int i = 0; i < (NX + 1) * (NY + 1) * (NZ + 1); i++) {				//перемещение узлов для следующей итерации
			for (int j = 0; j < 3; j++) {
				n[i].coordinate[j] += n[i].displacement[j];
			}
		}

		delete[] F;
		for (int count = 0; count < 3 * (NX + 1) * (NY + 1) * (NZ + 1); count++)
			delete[] K[count];
	}

	iteration++;

}

void ModelGL::createObject(int cNx, int cNy, int cNz, double A_, double B_, double C_, int Nx_, int Ny_, int Nz_) {
	
	cNX = cNx;
	cNY = cNy;
	cNZ = cNz;

	A = A_;
	B = B_;
	C = C_;

	NX = Nx_;
	NY = Ny_;
	NZ = Nz_;

	runFlag = true;
}

void ModelGL::setE(float E_) {
	E = E_;
}

void ModelGL::setMu(float Mu_) {
	Mu = Mu_;
}

void ModelGL::setdX(double dX_) {
	dX = dX_;
}

void ModelGL::setdY(double dY_) {
	dY = dY_;
}

void ModelGL::setdZ(double dZ_) {
	dZ = dZ_;
}

void ModelGL::setMoveX(bool move_X_) {
	moveble_X = move_X_;
}

void ModelGL::setMoveY(bool move_Y_) {
	moveble_Y = move_Y_;
}

void ModelGL::setMoveZ(bool move_Z_) {
	moveble_Z = move_Z_;
}

void ModelGL::setP(double P_) {
	P = P_;
}

void ModelGL::setP_bool(bool P_bool) {
	flag_enter_P = P_bool;
}

void ModelGL::save(const std::string name)
{
	ofstream F;
	F.open(name);

	F << NX << " " << NY << " " << NZ << " " << cNX << " " << cNY << " " << cNZ << endl;
	F << A << " " << B << " "<< C << endl;
	F << index_tetra << endl;
	F << endl << endl;
	for (int i = 0; i < cNX * cNY * cNZ; i++) {
		F << "Coordinate"<< " " << cub[i].coordinate[0] << " " << cub[i].coordinate[1] << " " << cub[i].coordinate[2] << endl;
		F << "Dimension" << " " <<cub[i].dimension[0] << " " << cub[i].dimension[1] << " " << cub[i].dimension[2] <<endl;
		F << "Colour" << " " << cub[i].colour[0] << " " << cub[i].colour[1] << " " << cub[i].colour[2] << endl;
		F << "E,MU" << " " << cub[i].E << " " << cub[i].mu << endl;
		F << "CurrCubFlag" << " " << cub[i].currentCub << endl;
		F << endl;
	}

	F << endl << endl;

	for (int i = 0; i < (NX + 1) * (NY + 1) * (NZ + 1); i++) {
		F << "Coordinate" << " " << n[i].coordinate[0] << " " << n[i].coordinate[1] << " " << n[i].coordinate[2]  << endl;
		F << "Displacement" << " " << n[i].displacement[0] << " " << n[i].displacement[1] << " " << n[i].displacement[2] << endl;
		F << "Boundary" << " " << n[i].boundary[0] << " "<< n[i].boundary[1] << " " << n[i].boundary[2] << endl;
		F << "Colour" << " " << n[i].colour[0] << " " << n[i].colour[1] << " " << n[i].colour[2] <<endl;
		F << "CurrNodeFlag" << " " << n[i].currentNode << endl;
		F << "PFlag" << " " << n[i].flag << endl;
		F << "NodeRectFlag" << " " << n[i].currentNodeRect << endl;
		F << endl;
	}

	for (int i = 0; i < 6 * NX * NY * NZ; i++) {
		F << "Coordinate" << " " << e[i].coordinate[0] << " " << e[i].coordinate[1] << " " << e[i].coordinate[2] << " " << e[i].coordinate[3] << endl;
		F << "Centr" << " " << e[i].centr[0] << " " << e[i].centr[1] << " " << e[i].centr[2] << endl;
		F << "Colour" << " " << e[i].colour[0] << " " << e[i].colour[1] << " " << e[i].colour[2] << endl;
		F << "Volume" << " " << e[i].volume << endl;
		F << "E,MU" << " " << e[i].E << " " << e[i].mu << endl;
		F << endl;
	}

	F << endl << endl;

	for (int i = 0; i < 4 * (NX * NY + NX * NZ + NY * NZ); i++) {
		F << "Coordinate" << " " << tr[i].coordinate[0] << " " << tr[i].coordinate[1] << " " << tr[i].coordinate[2] << endl;
		F << "Normal" << " " << tr[i].normal[0] << " " << tr[i].normal[1] << " " << tr[i].normal[2] <<endl;
		F << "Normal0" << " " << tr[i].normal0[0] << " " << tr[i].normal0[1] << " " << tr[i].normal0[2] << endl;
		F << endl;
	}

	F.close();
}
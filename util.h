#include <stdio.h>
#define _USE_MATH_DEFINES	// for constant M_PI
#include <math.h>
#include <stdlib.h>
#include "cv.h"
#include "cxcore.h"
#include "highgui.h"
#include "allocate.h"

#define INF			1000000000
#define SCALE		3

typedef struct point
{
	int		x;
	int		y;
} Point;

typedef struct dpoint
{
	double		x;
	double		y;
} DPoint;

//#define delta(a)	(double)((a)-(int)(a))
//#define real_coord(f,a,b)	(double)(f[(int)(a)+1][(int)(b)+1])*(delta(a))*(delta(b)) -(double)(f[(int)(a)][(int)(b)+1])*(delta(a)-1)*(delta(b)) -(double)(f[(int)(a)+1][(int)(b)])*(delta(b)-1)*(delta(a)) +(double)(f[(int)(a)][(int)(b)])*(delta(b)-1)*(delta(a)-1)
#define delta(a)	((a)-floor(a))
#define real_coord(f,a,b)	((double)(f[(int)(a)+1][(int)(b)+1])*(delta(a))*(delta(b)) -(double)(f[(int)(a)][(int)(b)+1])*(delta(a)-1)*(delta(b)) -(double)(f[(int)(a)+1][(int)(b)])*(delta(b)-1)*(delta(a)) +(double)(f[(int)(a)][(int)(b)])*(delta(b)-1)*(delta(a)-1))

void LoadImageFromMemory(IplImage *image, unsigned char **y);
double dist(DPoint a, DPoint b);

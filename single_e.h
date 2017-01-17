#include <stdio.h>
#define _USE_MATH_DEFINES	// for constant M_PI
#include <math.h>
#include <stdlib.h>
#include "cv.h"
#include "cxcore.h"
#include "highgui.h"
//#include "cvcam.h"
#include "cvaux.h"
#include "allocate.h"
#include "random.h"
#include "tiff.h"
#include "util.h"

//#define TEST_SINGLE_OBJECT

// Display level
#define SHOW_WINDOW	1
#define LEVEL0		1
#define LEVEL1		1
#define LEVEL2		1
#define LEVEL3		1
// for single object test
#ifdef TEST_SINGLE_OBJECT
	#define LEVEL4		0
#else
	#define LEVEL4		0
#endif
#define LEVEL5		0

#define TEXT_NONE			-1
#define TEXT_SINGLE_E		0
#define TEXT_MULTIPLE_E		1
#define TEXT_TOTAL_E		2
#define TEXT_BOTH_E			3
#define TEXT_E0_E1_E2		4
#define TEXT_NUMBER			5
#define TEXT_NUM_E0_E1_E2	6


#define MAX_POINT_NUM		500
#ifdef TEST_POLY
	#define POINT_NUM			15
#else
	#define POINT_NUM			200 // 200	//15
#endif
#define DELTA_T				2*M_PI/POINT_NUM
#define RMIN				4
#define RMAX				10//25
#define CURVE_MAX			100
#define GD_ITER_MAX			1000//150//50
#define MAX_SINGLE_ENERGY	-0
#define INTERNAL_AREA_UNDEF	0
#define INTERNAL_AREA_OUT	1
#define INTERNAL_AREA_IN	2

typedef struct curve_point
{
	Point	r;
	char	type;
} Curve_point;

typedef enum
{
	STATE_NON_EXIST,
	STATE_NEW_BORN,
	STATE_EXIST,
} MP_State;

typedef struct curve
{
	int			num;
	MP_State	state;
	DPoint		center;
	DPoint		r[POINT_NUM];
	DPoint		min;
	DPoint		max;
	double		single_E;	// single energy
	double		multiple_E;	// multiple energy
	char		interior[CURVE_MAX*2+1][CURVE_MAX*2+1];
	double		e0;
	double		e1;
	double		e2;
} Curve;

void Internal_Area(Curve *omega);
void draw_curve_box(Curve *omega, IplImage *image);
void draw_fill_curve(Curve *omega, IplImage *image);
void draw_curve(Curve *omega, IplImage *image, int energy_type);
void Qbirth(Curve *omega, DPoint center, double delta_t, int rmin, int rmax, int cols, int rows);
int Single_Object(unsigned char **y, Curve *omega, double GD_step, double lambda_g, double lambda_G,
				  double GDepsilon,  double delta_t, double *mean, double *vari, int cols, int rows,
				  IplImage *image, const char* win_name, CvVideoWriter *video);
void test_single_object(unsigned char **yimg, double GD_step,  
						double lambda_g, double lambda_G, double GDepsilon,
						double *mu, double *vari, int cols, int rows, 
						IplImage *image, const char* win_name, CvVideoWriter *video);
void test_single_object2(unsigned char **yimg, double GD_step,  
						double lambda_g, double lambda_G, double GDepsilon,
						double *mean, double *vari, int cols, int rows, 
						IplImage *image, const char* win_name, CvVideoWriter *video);

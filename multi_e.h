#include <stdio.h>
#define _USE_MATH_DEFINES	// for constant M_PI
#include <math.h>
#include <stdlib.h>
#include "cv.h"
#include "cxcore.h"
#include "highgui.h"
#include "allocate.h"
#include "random.h"
#include "tiff.h"
#include "single_e.h"

#define HARDCORE_REPULSION 3

#define MAX_MKPNT_NUM		220		// Max marked point number
#define ITER_STEP		10

void Internal_Area(Curve *omega);

void as_mpp(unsigned char **yimg, double **lm, double *mean, double *vari,
			double GD_step, double lambda_g, double lambda_G, double GDepsilon,
			double delta, double b_zero, double alpha, double beta,
			double F, int cols, int rows, IplImage *image, const char* win_name,
			CvVideoWriter *video);

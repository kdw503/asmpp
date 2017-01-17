#include "multi_e.h"

// Fill
void asmpp_Qdelete(Curve *omega, int k, int *curve_num)
{
	int j;
	for (j = k; j < *curve_num ; j++)
			omega[j] = omega[j+1];
	(*curve_num)--;
}

void draw_all_curves(Curve *omega, int curve_num, IplImage *image, int energy_type)
{
	int k;
	CvScalar red = CV_RGB(255,0,0);

	for(k = 0; k<curve_num; k++)
		draw_curve(&(omega[k]), image, energy_type, red);
}

int compare_Single_E (const void * a, const void * b)
{
	if(((Curve*)a)->single_E > ((Curve*)b)->single_E) return -1;
	else if(((Curve*)a)->single_E < ((Curve*)b)->single_E) return 1;
	else return 0;
}

double asmpp_intersection_area(Curve *omega1, Curve *omega2)
{
	int num1, num2, in_num;
	int i,j,m,n;
	int size1x,size1y,size2x,size2y;

	num1 = 0;
	in_num = 0;
	size1x = (int)(omega1->max.x - omega1->min.x + 3);
	size1y = (int)(omega1->max.y - omega1->min.y + 3);
	size2x = (int)(omega2->max.x - omega2->min.x + 3);
	size2y = (int)(omega2->max.y - omega2->min.y + 3);
	for (i=0;i < size1y ;i++){
		for (j=0;j < size1x ;j++){
			if(omega1->interior[i][j] == INTERNAL_AREA_IN)
			{
				num1++;
				m = i + (int)(omega1->min.y + omega1->center.y - omega2->min.y - omega2->center.y);
				n = j + (int)(omega1->min.x + omega1->center.x - omega2->min.x - omega2->center.x);
				if((m>=0)&&(m<size2y)&&(n>=0)&&(n<size2x))
				{
					if(omega2->interior[m][n] == INTERNAL_AREA_IN)
						in_num++;
				}
			}
		}
	}
	num2 = 0;
	for (i=0;i < size2y ;i++){
		for (j=0;j < size2x ;j++){
			if(omega2->interior[i][j] == INTERNAL_AREA_IN)
			{
				num2++;
			}
		}
	}
	if(num2<num1)
		num1 = num2;
	return (double)in_num/(double)num1;
}

double asmpp_C_prior (Curve *omega, int k, int curve_num, double epsilon)
{	
	int i;
	double sum,dtmp1,dtmp2,dist1,dist2;

	dtmp1 = sqrt(omega[k].max.x * omega[k].max.x + omega[k].max.y * omega[k].max.y);
	dtmp2 = sqrt(omega[k].min.x * omega[k].min.x + omega[k].min.y * omega[k].min.y);
	if(dtmp2<dtmp1) dist1 = dtmp1;
	else dist1 = dtmp2;

	sum = 0.;
	for(i=0;i<curve_num;i++){
		if(i!=k){
			dtmp1 = sqrt(omega[i].max.x * omega[i].max.x + omega[i].max.y * omega[i].max.y);
			dtmp2 = sqrt(omega[i].min.x * omega[i].min.x + omega[i].min.y * omega[i].min.y);
			if(dtmp2<dtmp1) dist2 = dtmp1;
			else dist2 = dtmp2;
			dtmp1 = dist(omega[k].center,omega[i].center);
			if (dtmp1<epsilon)
				return INF;
			if(dtmp1<dist1+dist2){
				sum += asmpp_intersection_area(&(omega[k]), &(omega[i]));
			}
		}
	}
	return sum;
}

/******************************************************************************/
// as_mpp : 
// 
/******************************************************************************/
void as_mpp(unsigned char **yimg, unsigned char **lm, double *mean, double *vari,
			double GD_step, double lambda_g, double lambda_G, double GDepsilon,
			double delta, double b_zero, double alpha, double beta,
			double F, int cols, int rows, IplImage *image, const char* win_name, CvVideoWriter *video)
{
	double test_pointx[15] = {6,-39,-41,-33,-23,-23,-57,-49,-28,1,0,49,50,13,60};
	double test_pointy[15] = {-59,-38,-15,-23,-23,10,14,51,36,48,77,82,22,-13,-29};
	double test_centerx = 63;
	double test_centery = 97;
	unsigned char **MP_exist;
	DPoint r;
	int i, j, ii, jj, m, n, k, iter;
	double min, max, dtmp;
	double birth_rate, d_beta, d_rate;
	Curve omega[MAX_MKPNT_NUM];
	int curve_num, succeed;
	double epsilon=HARDCORE_REPULSION;
	int mpp_iter_stop;
	double delta_t = DELTA_T;
	CvScalar red = CV_RGB(255,0,0);
	IplImage *image2 = 0;
	char win_name2[256] = "Birth";
	int channel = 3;
	int width=cols*SCALE, height=rows*SCALE;
	CvVideoWriter *video0001= cvCreateVideoWriter("video0001.avi",
		0, 30, cvSize(image->width,image->height),1);


	image2 = cvCreateImage(cvSize(width, height), IPL_DEPTH_8U, channel);
	LoadImageFromMemory(image2, yimg);
	if(SHOW_WINDOW*LEVEL3){
		// create a window
		cvNamedWindow(win_name2, CV_WINDOW_AUTOSIZE);
		cvResizeWindow( win_name2, width, height);
		cvMoveWindow(win_name2, 1000, 500);
		cvShowImage(win_name2, image2);
		// wait for a key
		cvWaitKey(1);
	}

	MP_exist = (unsigned char **)get_img(cols, rows, sizeof(unsigned char));
	for (i = 0; i < rows; i++)
		for (j = 0; j < cols; j++){
			MP_exist[i][j] = 0;
		}
	for(i = 0; i < MAX_MKPNT_NUM; i++){
		omega[i].state = STATE_NON_EXIST;
	}
	birth_rate = delta*b_zero;
	curve_num = 0;
	iter = 0;
	for (m = 0;m<ITER_STEP;m+=2){
	for (n = 0;n<ITER_STEP;n+=2){
	//----------------------- Birth -------------------------//ITER_STEP
		for (ii = 3;ii<rows-3;ii+=ITER_STEP){
			i = ii+m;
			if(i>rows-3) i = rows-3;
			for (jj = 3;jj<cols-3;jj+=ITER_STEP){
			j = jj+m;
			if(j>cols-3) j = cols-3;
				if (MP_exist[i][j]==0){
					dtmp = ((double)rand()/RAND_MAX);
					if (dtmp < birth_rate*(double)(lm[i][j])){ // birth
						if (curve_num>MAX_MKPNT_NUM-2) 
							break;
						r.x = (double)j;
						r.y = (double)i;
						Qbirth(&(omega[curve_num]), r, delta_t, RMIN, RMAX, cols, rows);
						succeed = Single_Object(yimg, &(omega[curve_num]), GD_step, lambda_g,lambda_G, GDepsilon,
										delta_t, mean, vari, cols, rows, image2, win_name2, video0001);
						if ((succeed)&&(omega[curve_num].single_E < MAX_SINGLE_ENERGY))
						{
							MP_exist[(int)(omega[curve_num].center.y)][(int)(omega[curve_num].center.x)] ++;
							omega[curve_num].state = STATE_NEW_BORN;
							if(SHOW_WINDOW*LEVEL3){
								draw_curve(&(omega[curve_num]), image2, TEXT_SINGLE_E, red);
								cvShowImage(win_name2, image2);
								cvWaitKey((int)(1));
							}
							curve_num ++;
						}
						// Calculate Interior Area
						Internal_Area(&(omega[curve_num]));
					}
				}
			}
		}
		//----------------------- death -------------------------//
		qsort (omega, curve_num, sizeof(Curve), compare_Single_E);

		for(k = 0; k<curve_num; k++){
			omega[k].multiple_E = asmpp_C_prior (omega, k, curve_num, epsilon);
			if(SHOW_WINDOW*LEVEL2){
				draw_curve(&(omega[k]), image, TEXT_BOTH_E, red);
			}
			if(omega[k].multiple_E>=INF) 
				d_rate = 1.;
			else{
				d_beta = exp(-beta*(-alpha*omega[k].single_E-omega[k].multiple_E));
				d_rate = (delta*d_beta)/(1.0+delta*d_beta);
//					printf("k=%d, dp=%1.2f\n",k,d_rate);
			}

			dtmp = ((double)rand()/RAND_MAX);
			if (dtmp < d_rate){
				MP_exist[(int)omega[k].center.y][(int)omega[k].center.x] --;
				asmpp_Qdelete(omega, k, &curve_num);
				k--;
			}
			else
				omega[k].num = k;
		}
		if(SHOW_WINDOW*LEVEL2){
			cvShowImage(win_name2, image2);
			cvWaitKey((int)(1));
			LoadImageFromMemory(image2, yimg);
			LoadImageFromMemory(image, yimg);
			draw_all_curves(omega, curve_num, image, TEXT_BOTH_E);
			cvShowImage(win_name, image);
			cvWaitKey((int)(1));
			cvWriteFrame( video, image );
			cvWriteFrame( video, image );
			cvWriteFrame( video, image );
			cvWriteFrame( video, image );
			cvWriteFrame( video, image );
			cvWriteFrame( video0001, image2 );
		}
		//---------------- Set up new parameters ----------------//
		delta = delta*F;
		beta = beta/F;
		printf("iter = %d\n",iter);
		iter ++;
	}
	}
	cvReleaseVideoWriter( &video0001 );

	if(SHOW_WINDOW*LEVEL1){
		LoadImageFromMemory(image, yimg);
		draw_all_curves(omega, curve_num, image, TEXT_BOTH_E);
		cvShowImage(win_name, image);
		cvWaitKey((int)(1));
	}
	free_img((void **)MP_exist);
}

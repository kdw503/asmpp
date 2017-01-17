#include <stdio.h>
#define _USE_MATH_DEFINES	// for constant M_PI
#include <math.h>
#include <stdlib.h>
#include "asmpp.h"


/******************************************************************************/
// main
// return : 
// TODO
// - change the mean an variance with the ground truth image or JDS emmpm
// - r : change to gaussian distribution
// - beta control when mpp are overlapped
/******************************************************************************/
int main( int argc , char** argv)
{
	struct TIFF_img input_img, input_gt_img, output_img;
	unsigned char **yimg, **laplacian;
	char infileName[1024], outfileName[1024];
	char segfileName[1024], outfilePrefix[1024];
	char win_name[256] = "input data";
	FILE *fp;
	IplImage *image = 0;
	Curve omega[MAX_MKPNT_NUM];
	int height, width, channel=3;
	int i, j, rows, cols;
	double birth_rate, lambda_g, lambda_G, GDepsilon;
	double GD_step, delta_beta;
//	double mu0 = 35.9, mu1 = 195.0, vari0 = 758.5, vari1 = 1970.5; // for synthetic_large.tiff
//	double mu0 = 20, mu1 = 230, vari0 = 100, vari1 = 100; // for synthetic_large.tiff
	double mu[CLASSES], mean, vari[CLASSES], variance;
	int cnt[CLASSES];
	// for EM/MPM
	double beta[MAX_CLASSES], gamma[MAX_CLASSES];
	int mpmiter = 10, emiter = 30, classes = 2;
	unsigned char **xt, **gt, xt255;  /* output : entropy image */
	double **blur, sigma = 0., dsum, sum[CLASSES], di, dj, misclassed;
	int blur_size = 5, enable_blur = 0;
	int run_emmpm = 0;
	double b_zero = 1.0;
	double delta = 0.9;
	double betampp = 10;
	double F = 0.98;
	double alpha = 0.2;
	CvVideoWriter *video0000;

	//*****************************************************************************
	//		Initialize parameters
	//*****************************************************************************
	if (argc != 9) {
		printf("usage:  %s infileprefix t_clust en_suscept\n", argv[0]);
		//exit(1);
		argv[1] = "img1_t_crop";   // synthetic_large				"synthetic_blur";	"slice_crop017";  
		argv[2] = "30";				 // lambda_g(2.9)				"20";				"40";								
		argv[3] = "400";			 // lambda_G(2.9)				"100";				"270";						
		argv[4] = "0.000000921";		 // GDepsilon = 1 (0.000921)	"0.0000921";		"0.0000921";				
		argv[5] = "1.007";			 // delta_beta = 1/0.993		"1.007";			"1.007";						
		argv[6] = "0.001";			 // birth_rate					"0.001";			"0.001";					
		argv[7] = "0.00010";		 // GD_step = 0.99 (0.0004)		"0.0005";			"0.00025";	(100point,center=(120,95),r=10)			
		argv[8] = "0.1";				 // alpha (0.25) 			
	}
	sprintf(infileName, "%s.tiff",argv[1]);
	sprintf(segfileName, "%s_seg.tiff",argv[1]);
	sprintf(outfilePrefix, "%s",argv[1]);
	lambda_g = atof(argv[2]);
	lambda_G = atof(argv[3]);
	GDepsilon = atof(argv[4]);
	delta_beta = atof(argv[5]);
	birth_rate = atof(argv[6]);
	GD_step = atof(argv[7]);
	alpha = atof(argv[8]);

	readseed();

	//*****************************************************************************
	//		read input image
	//*****************************************************************************
	if ((fp = fopen(infileName, "rb")) == NULL) {
		printf("Cannot open file %s\n", infileName);
		exit(1);
	}
	else {
		if (read_TIFF(fp, &input_img)) { 	/* read image */
			printf("Error reading file %s\n", infileName);
			exit(1);
		}
		cols = input_img.width;
		rows = input_img.height;
		height = rows*SCALE;	// for window
		width = cols*SCALE;
		printf ("cols = %d, rows = %d\n",cols, rows);
		yimg = (unsigned char **)get_img(cols, rows, sizeof(unsigned char));
		image = cvCreateImage(cvSize(width, height), IPL_DEPTH_8U, channel);
		laplacian = (unsigned char **)get_img(cols, rows, sizeof(unsigned char));
		gt = (unsigned char **)get_img(cols, rows, sizeof(char));
		xt = (unsigned char **)get_img(cols, rows, sizeof(char));
		blur = (double **)get_img(blur_size, blur_size, sizeof(double));
		get_TIFF(&output_img, rows, cols, 'g');
#ifdef TEST_SINGLE_OBJECT
		video0000 = cvCreateVideoWriter("video0000.avi", 0, 30, cvSize(image->width,image->height),1);
#else
		video0000 = cvCreateVideoWriter("video0000.avi", 0, 5, cvSize(image->width,image->height),1);
#endif
		/* close image file */
		fclose(fp);

		/* check the type of image data */
		if (input_img.TIFF_type != 'g') {
			printf("Error:  Image must be grayscale.\n");
			exit(1);
		}
		for (i = 0; i < rows; i++)
			for (j = 0; j < cols; j++){
				yimg[i][j] = input_img.mono[i][j];
			}
//		for (i = 1; i < rows-1; i++)
//			for (j = 1; j < cols-1; j++){
//				tmp = yimg[i+1][j]+yimg[i-1][j]+yimg[i][j+1]+yimg[i][j-1]-4*yimg[i][j];
//				if(tmp==0)
//					laplacian[i][j] = 0;
//				else
//					laplacian[i][j] = tmp+128;
//			}
	}

	LoadImageFromMemory(image, yimg);
	if(SHOW_WINDOW*LEVEL0){
		// create a window
		cvNamedWindow(win_name, CV_WINDOW_AUTOSIZE);
		cvResizeWindow( win_name, width, height);
		cvMoveWindow(win_name, 1000, 50);
		cvShowImage(win_name, image);
		printf("pause\n");
		// wait for a key
		cvWaitKey(0);
	}

	if ((fp = fopen(segfileName, "rb")) == NULL) {
		printf("Cannot open file %s\n", infileName);
		printf("Execute EM/MPM Segmentation..\n");
		run_emmpm = 1;
	}
	else {
		if (read_TIFF(fp, &input_img)) { 	/* read image */
			printf("Error reading file %s\n", infileName);
			printf("Execute EM/MPM Segmentation..\n");
			run_emmpm = 1;
		}
		else{
			for (i = 0; i < rows; i++)
				for (j = 0; j < cols; j++){
					// segmentation result file to use as lebesgue measure
					xt[i][j] = (int)(input_img.mono[i][j]*(classes - 1)/255) ; 
				}
			run_emmpm = 0;
			if(SHOW_WINDOW*LEVEL1){
				LoadImageFromMemory(image, input_img.mono);
				// create a window
				cvShowImage(win_name, image);
				printf("pause\n");
				// wait for a key
				cvWaitKey(0);
			}
		}
	}

	//*****************************************************************************
	//		EM/MPM for Lebesgue measuring and learn mu0 vari0 mu1 vari1
	//*****************************************************************************
	if (run_emmpm){
		// initialize for EM/MPM
		//max_entropy = log10(classes)/log(2);

		for (i=0; i<CLASSES; i++){
			beta[i] = 3.0;
			gamma[i] = 0;
		}
		for(i = 0; i < MAX_MKPNT_NUM; i++){
			omega[i].state = STATE_NON_EXIST;
		}

		// calculate blurring matrix 
		/*	Aexp(-x/sigma)											*/
		/*	A(amplitude), x(distance from center), sigma(width)		*/
		if(sigma == 0.)
			enable_blur = 0;
		else{
			enable_blur = 1;
			dsum = 0.;
			for (i = 0; i < blur_size; i++){
				for (j = 0; j < blur_size; j++){
					di = i-(blur_size-1)/2.;
					dj = j-(blur_size-1)/2.;
					blur[i][j] = exp(-sqrt(di*di+dj*dj)/sigma);
					dsum += blur[i][j];
				}
			}
			for (i = 0; i < blur_size; i++){
				for (j = 0; j < blur_size; j++){
					blur[i][j] = blur[i][j]/dsum;
					printf("%1.4f  ", blur[i][j]);
				}
				printf("\n");
			}
		}

		emmpm(yimg, xt, beta, gamma, emiter, mpmiter, rows, cols, classes, blur, blur_size, enable_blur);

		for (i=0; i<rows; i++)
			for (j=0; j<cols; j++)
				output_img.mono[i][j] = (int)xt[i][j] * 255 / (classes - 1);
		if ((fp = fopen(segfileName, "wb")) == NULL ) {
			printf("Cannot open file %s\n", outfileName);
			exit(1);
		}
		if (write_TIFF(fp, &output_img)) {
			printf("Error writing TIFF file %s\n", outfileName);
			exit(1);
		}
		fclose(fp);
		if(SHOW_WINDOW*LEVEL1){
			LoadImageFromMemory(image, output_img.mono);
			// create a window
			cvShowImage(win_name, image);
			printf("pause\n");
			// wait for a key
			cvWaitKey(0);
			LoadImageFromMemory(image, yimg);
		}
	}
	
	dsum = 0.;
	for (i=0; i<CLASSES; i++){
		sum[i] = 0.;
		cnt[i] = 0;
	}
	for (i=0; i<rows; i++){
		for (j=0; j<cols; j++){
			sum[xt[i][j]] += yimg[i][j];
			cnt[xt[i][j]] ++;
			dsum += yimg[i][j];
		}
	}
	mean = dsum/(double)(rows*cols);
	dsum = 0.;
	for (i=0; i<CLASSES; i++){
		mu[i] = sum[i]/(double)cnt[i];
		sum[i] = 0.;
	}
	for (i=0; i<rows; i++){
		for (j=0; j<cols; j++){
			sum[xt[i][j]] += (yimg[i][j]-mu[xt[i][j]])*(yimg[i][j]-mu[xt[i][j]]);
			dsum += (yimg[i][j]-mean)*(yimg[i][j]-mean);
		}
	}
	variance = dsum/(double)(rows*cols);
	for (i=0; i<CLASSES; i++){
		vari[i] = sum[i]/(double)cnt[i];
		printf ("mu[%d] = %1.3f, vari[%d] = %1.3f\n",i,mu[i],i,vari[i]);
	}
	
	//*****************************************************************************
	//		Arbitrary Shape MPP
	//*****************************************************************************
	vari[0] = variance;
	vari[1] = variance;
#ifdef TEST_SINGLE_OBJECT
	test_single_object2(yimg, GD_step, lambda_g,lambda_G, GDepsilon, mu, vari, cols, rows, image, win_name, video0000);
#else
	as_mpp(yimg, xt, mu, vari,GD_step, lambda_g, lambda_G, GDepsilon,
			delta, b_zero, alpha, betampp, F, cols, rows, image, win_name, video0000);
#endif

	// printf("curve_num = %d ",curve_num);

	//*****************************************************************************
	//		Evaluate the result
	//*****************************************************************************
	// PMP
	misclassed = 0;
	for (i = 0; i < rows; i++)
		for (j = 0; j < cols; j++){
			xt255 = (unsigned char)xt[i][j] * 255 / (classes - 1);
			if(xt255!=gt[i][j]) misclassed += 1.0;
		}
	misclassed = misclassed*100/(cols*rows); // percentage of misclassified pixels
	printf("misclassed = %1.4f\n", misclassed);
	for (i=0; i<rows; i++)
		for (j=0; j<cols; j++)
			output_img.mono[i][j] = (int)xt[i][j] * 255 / (classes - 1);
	sprintf(outfileName, "%s_asmpp_p%1.2f.tiff",outfilePrefix, misclassed);
	if ((fp = fopen(outfileName, "wb")) == NULL ) {
		printf("Cannot open file %s\n", outfileName);
		exit(1);
	}
	if (write_TIFF(fp, &output_img)) {
		printf("Error writing TIFF file %s\n", outfileName);
		exit(1);
	}
	fclose(fp);


	if(SHOW_WINDOW*LEVEL0){
		// create a window
//		cvNamedWindow("input data", CV_WINDOW_AUTOSIZE);
//		cvResizeWindow( "input data", width, height);
//		cvMoveWindow("input data", 200, 200);
		cvShowImage("input data", image);
		printf("pause1\n");
		// wait for a key
		cvWaitKey(0);
/*
		sprintf(filename, "%s_jij.jpg",argv[1]);
		if(!cvSaveImage(filename, image, 0)){
			printf("Could not save: %s\n",filename);
		}
*/
		cvReleaseImage(&image);
		cvDestroyWindow("input data");
	}
	cvReleaseVideoWriter( &video0000 );
	free_TIFF(&output_img);
	free_TIFF(&input_img);
	free_img((void **)yimg);
	free_img((void **)laplacian);
	free_img((void **)gt);
	free_img((void **)xt);
	free_img((void **)blur);

	writeseed();

	return 0;
}


#if 0
#include <stdio.h>
 #include "cv.h"
 #include "highgui.h"

 int main( int argc, char** argv )
 {
 IplImage *frame;
 int key;
 assert( argc == 2 );
 CvCapture *capture = cvCaptureFromAVI( argv[1] );
 if( !capture ) return 1; /*THE CODE EXITS AT THIS POINT*/
 int fps = ( int )cvGetCaptureProperty( capture, CV_CAP_PROP_FPS );
 cvNamedWindow( "video", 0 );
 while( key != 'q' ) {
 frame = cvQueryFrame( capture );
 if( !frame ){printf("\n no frame\n"); break;}
 cvShowImage( "video", frame );
 key = cvWaitKey( 1000 / fps );
 }
 cvReleaseCapture( &capture );
 cvDestroyWindow( "video" );

 return 0;
 }
#endif
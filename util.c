#include "util.h"

// Load Image to the Image memory
void LoadImageFromMemory(IplImage *image, unsigned char **y)
{
	int i,j,k;
	int step, channel, height, width;
	uchar *imgdata;
	double di,dj;

	step = image->widthStep;
	imgdata = (uchar *)image->imageData;
	channel = image->nChannels;
	height = image->height;
	width = image->width;
	for(i = 0; i < height-SCALE; i++){
		for(j = 0; j < width-SCALE; j++){
			for(k = 0; k < channel; k++){
				di = (double)i/SCALE;
				dj = (double)j/SCALE;
				imgdata[i*step + j*channel + k] = (uchar)(real_coord(y,di,dj));
			}
		}
	}
}

double dist(DPoint a, DPoint b)
{
	return sqrt((a.x - b.x)*(a.x - b.x)+(a.y - b.y)*(a.y - b.y));
}

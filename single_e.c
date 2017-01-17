#include "single_e.h"


void fill_external(Curve *omega,int i,int j, int sizey, int sizex)
{
	if(i+1<sizey){
		if (omega->interior[i+1][j]==INTERNAL_AREA_UNDEF){
			omega->interior[i+1][j] = INTERNAL_AREA_OUT;
			fill_external(omega,i+1,j,sizey,sizex);
		}
	}
	if(i-1>=0){
		if (omega->interior[i-1][j]==INTERNAL_AREA_UNDEF){
			omega->interior[i-1][j] = INTERNAL_AREA_OUT;
			fill_external(omega,i-1,j,sizey,sizex);
		}
	}
	if(j+1<sizex){
		if (omega->interior[i][j+1]==INTERNAL_AREA_UNDEF){
			omega->interior[i][j+1] = INTERNAL_AREA_OUT;
			fill_external(omega,i,j+1,sizey,sizex);
		}
	}
	if(j-1>=0){
		if (omega->interior[i][j-1]==INTERNAL_AREA_UNDEF){
			omega->interior[i][j-1] = INTERNAL_AREA_OUT;
			fill_external(omega,i,j-1,sizey,sizex);
		}
	}
}

// Calculate Internal Area
void Internal_Area_old(Curve *omega)
{	
	int i,j,k,sizex,sizey;
	double dtmp,x,y;

	// Initialize
	sizex = (int)(omega->max.x - omega->min.x + 3);
	sizey = (int)(omega->max.y - omega->min.y + 3);
	for (i=0;i < sizey ;i++){
		for (j=0;j < sizex ;j++){
			omega->interior[i][j] = INTERNAL_AREA_UNDEF;
		}
	}

	for (k= 0; k < POINT_NUM -1; k++){
		if(omega->r[k+1].x == omega->r[k].x){
			if(omega->r[k].y < omega->r[k+1].y){
				for (y = omega->r[k].y; y <= omega->r[k+1].y; y++){
					i = (int)(y-omega->min.y+1);
					j = (int)(omega->r[k].x-omega->min.x+1);
					omega->interior[i][j] = INTERNAL_AREA_IN;
//					data_pnt.y = y+omega->center.y;
//					data_pnt.y = height - data_pnt.y;
//					data_pnt.x = omega->r[k].x+omega->center.x;
//					cvCircle(img, data_pnt, 1, CV_RGB(255,0,0), 2, 8, 0);
				}
			}
			else{
				for (y = omega->r[k+1].y; y <= omega->r[k].y; y++){
					i = (int)(y-omega->min.y+1);
					j = (int)(omega->r[k].x-omega->min.x+1);
					omega->interior[i][j] = INTERNAL_AREA_IN;
				}
			}
		}
		else{
			dtmp = (omega->r[k+1].y-omega->r[k].y)/(omega->r[k+1].x-omega->r[k].x);
			if(fabs(dtmp)<1){
				if(omega->r[k].x < omega->r[k+1].x){
					for (x = omega->r[k].x; x <= omega->r[k+1].x; x++){
						i = (int)(dtmp*(x-omega->r[k].x)+omega->r[k].y-omega->min.y+1);
						j = (int)(x-omega->min.x+1);
						omega->interior[i][j] = INTERNAL_AREA_IN;
					}
				}
				else{
					for (x = omega->r[k+1].x; x <= omega->r[k].x; x++){
						i = (int)(dtmp*(x-omega->r[k].x)+omega->r[k].y-omega->min.y+1);
						j = (int)(x-omega->min.x+1);
						omega->interior[i][j] = INTERNAL_AREA_IN;
					}
				}
			}
			else{
				if(omega->r[k].y < omega->r[k+1].y){
					for (y = omega->r[k].y; y <= omega->r[k+1].y; y++){
						i = (int)(y-omega->min.y+1);
						j = (int)(1/dtmp*(y-omega->r[k].y)+omega->r[k].x-omega->min.x+1);
						omega->interior[i][j] = INTERNAL_AREA_IN;
					}
				}
				else{
					for (y = omega->r[k+1].y; y <= omega->r[k].y; y++){
						i = (int)(y-omega->min.y+1);
						j = (int)(1/dtmp*(y-omega->r[k].y)+omega->r[k].x-omega->min.x+1);
						omega->interior[i][j] = INTERNAL_AREA_IN;
					}
				}
			}
		}
	}
	if(omega->r[0].x == omega->r[k].x){
		if(omega->r[k].y < omega->r[0].y){
			for (y = omega->r[k].y; y <= omega->r[0].y; y++){
				i = (int)(y-omega->min.y+1);
				j = (int)(omega->r[k].x-omega->min.x+1);
				omega->interior[i][j] = INTERNAL_AREA_IN;
			}
		}
		else{
			for (y = omega->r[0].y; y <= omega->r[k].y; y++){
				i = (int)(y-omega->min.y+1);
				j = (int)(omega->r[k].x-omega->min.x+1);
				omega->interior[i][j] = INTERNAL_AREA_IN;
			}
		}
	}
	else{
		dtmp = (omega->r[0].y-omega->r[k].y)/(omega->r[0].x-omega->r[k].x);
		if(fabs(dtmp)<1){
			if(omega->r[k].x < omega->r[0].x){
				for (x = omega->r[k].x; x <= omega->r[0].x; x++){
					i = (int)(dtmp*(x-omega->r[k].x)+omega->r[k].y-omega->min.y+1);
					j = (int)(x-omega->min.x+1);
					omega->interior[i][j] = INTERNAL_AREA_IN;
				}
			}
			else{
				for (x = omega->r[0].x; x <= omega->r[k].x; x++){
					i = (int)(dtmp*(x-omega->r[k].x)+omega->r[k].y-omega->min.y+1);
					j = (int)(x-omega->min.x+1);
					omega->interior[i][j] = INTERNAL_AREA_IN;
				}
			}
		}
		else{
			if(omega->r[k].y < omega->r[0].y){
				for (y = omega->r[k].y; y <= omega->r[0].y; y++){
					i = (int)(y-omega->min.y+1);
					j = (int)(1/dtmp*(y-omega->r[k].y)+omega->r[k].x-omega->min.x+1);
					omega->interior[i][j] = INTERNAL_AREA_IN;
				}
			}
			else{
				for (y = omega->r[0].y; y <= omega->r[k].y; y++){
					i = (int)(y-omega->min.y+1);
					j = (int)(1/dtmp*(y-omega->r[k].y)+omega->r[k].x-omega->min.x+1);
					omega->interior[i][j] = INTERNAL_AREA_IN;
				}
			}
		}
	}

	omega->interior[0][0] = INTERNAL_AREA_OUT;
	fill_external(omega,0,0,sizey,sizex);

	for (i=0;i<sizey;i++){
		for (j=0;j<sizex;j++){
			if(omega->interior[i][j] == INTERNAL_AREA_UNDEF){
				omega->interior[i][j] = INTERNAL_AREA_IN;
//				data_pnt.y = (int)((i+omega->min.y-1+omega->center.y)*(double)SCALE);
//				data_pnt.x = (int)((j+omega->min.x-1+omega->center.x)*(double)SCALE);
//				cvCircle(img, data_pnt, 1, CV_RGB(0,255,0), 2, 8, 0);
			}
//			else if(omega->interior[i][j] == INTERNAL_AREA_IN){
//				data_pnt.y = (int)((i+omega->min.y-1+omega->center.y)*(double)SCALE);
//				data_pnt.x = (int)((j+omega->min.x-1+omega->center.x)*(double)SCALE);
//				cvCircle(img, data_pnt, 1, CV_RGB(255,0,0), 2, 8, 0);
//			}
		}
	}
	return;
}

// Calculate Internal Area
void Internal_Area(Curve *omega)
{	
	int i,j,k,sizex,sizey;
	double dtmp,x,y;
	Point p[POINT_NUM], min, max;

	for (k= 0; k < POINT_NUM; k++){
		p[k].x = (int)floor(omega->r[k].x+0.5);
		p[k].y = (int)floor(omega->r[k].y+0.5);
	}
	min.x = (int)floor(omega->min.x+0.5);
	min.y = (int)floor(omega->min.y+0.5);
	max.x = (int)floor(omega->max.x+0.5);
	max.y = (int)floor(omega->max.y+0.5);
	// Initialize
	sizex = (int)(max.x - min.x + 3);
	sizey = (int)(max.y - min.y + 3);
	for (i=0;i < sizey ;i++){
		for (j=0;j < sizex ;j++){
			omega->interior[i][j] = INTERNAL_AREA_UNDEF;
		}
	}

	for (k= 0; k < POINT_NUM -1; k++){
		if(p[k+1].x == p[k].x){
			if(p[k].y < p[k+1].y){
				for (y = p[k].y; y <= p[k+1].y; y++){
					i = (int)(y-min.y+1);
					j = (int)(p[k].x-min.x+1);
					omega->interior[i][j] = INTERNAL_AREA_IN;
//					data_pnt.y = y+omega->center.y;
//					data_pnt.y = height - data_pnt.y;
//					data_pnt.x = p[k].x+omega->center.x;
//					cvCircle(img, data_pnt, 1, CV_RGB(255,0,0), 2, 8, 0);
				}
			}
			else{
				for (y = p[k+1].y; y <= p[k].y; y++){
					i = (int)(y-min.y+1);
					j = (int)(p[k].x-min.x+1);
					omega->interior[i][j] = INTERNAL_AREA_IN;
				}
			}
		}
		else{
			dtmp = (double)(p[k+1].y-p[k].y)/(double)(p[k+1].x-p[k].x);
			if(fabs(dtmp)<1){
				if(p[k].x < p[k+1].x){
					for (x = p[k].x; x <= p[k+1].x; x++){
						i = (int)(dtmp*(x-p[k].x)+p[k].y-min.y+1);
						j = (int)(x-min.x+1);
						omega->interior[i][j] = INTERNAL_AREA_IN;
					}
				}
				else{
					for (x = p[k+1].x; x <= p[k].x; x++){
						i = (int)(dtmp*(x-p[k].x)+p[k].y-min.y+1);
						j = (int)(x-min.x+1);
						omega->interior[i][j] = INTERNAL_AREA_IN;
					}
				}
			}
			else{
				if(p[k].y < p[k+1].y){
					for (y = p[k].y; y <= p[k+1].y; y++){
						i = (int)(y-min.y+1);
						j = (int)(1/dtmp*(y-p[k].y)+p[k].x-min.x+1);
						omega->interior[i][j] = INTERNAL_AREA_IN;
					}
				}
				else{
					for (y = p[k+1].y; y <= p[k].y; y++){
						i = (int)(y-min.y+1);
						j = (int)(1/dtmp*(y-p[k].y)+p[k].x-min.x+1);
						omega->interior[i][j] = INTERNAL_AREA_IN;
					}
				}
			}
		}
	}
	if(p[0].x == p[k].x){
		if(p[k].y < p[0].y){
			for (y = p[k].y; y <= p[0].y; y++){
				i = (int)(y-min.y+1);
				j = (int)(p[k].x-min.x+1);
				omega->interior[i][j] = INTERNAL_AREA_IN;
			}
		}
		else{
			for (y = p[0].y; y <= p[k].y; y++){
				i = (int)(y-min.y+1);
				j = (int)(p[k].x-min.x+1);
				omega->interior[i][j] = INTERNAL_AREA_IN;
			}
		}
	}
	else{
		dtmp = (double)(p[0].y-p[k].y)/(double)(p[0].x-p[k].x);
		if(fabs(dtmp)<1){
			if(p[k].x < p[0].x){
				for (x = p[k].x; x <= p[0].x; x++){
					i = (int)(dtmp*(x-p[k].x)+p[k].y-min.y+1);
					j = (int)(x-min.x+1);
					omega->interior[i][j] = INTERNAL_AREA_IN;
				}
			}
			else{
				for (x = p[0].x; x <= p[k].x; x++){
					i = (int)(dtmp*(x-p[k].x)+p[k].y-min.y+1);
					j = (int)(x-min.x+1);
					omega->interior[i][j] = INTERNAL_AREA_IN;
				}
			}
		}
		else{
			if(p[k].y < p[0].y){
				for (y = p[k].y; y <= p[0].y; y++){
					i = (int)(y-min.y+1);
					j = (int)(1/dtmp*(y-p[k].y)+p[k].x-min.x+1);
					omega->interior[i][j] = INTERNAL_AREA_IN;
				}
			}
			else{
				for (y = p[0].y; y <= p[k].y; y++){
					i = (int)(y-min.y+1);
					j = (int)(1/dtmp*(y-p[k].y)+p[k].x-min.x+1);
					omega->interior[i][j] = INTERNAL_AREA_IN;
				}
			}
		}
	}

	omega->interior[0][0] = INTERNAL_AREA_OUT;
	fill_external(omega,0,0,sizey,sizex);

	for (i=0;i<sizey;i++){
		for (j=0;j<sizex;j++){
			if(omega->interior[i][j] == INTERNAL_AREA_UNDEF){
				omega->interior[i][j] = INTERNAL_AREA_IN;
//				data_pnt.y = (int)((i+min.y-1+omega->center.y)*(double)SCALE);
//				data_pnt.x = (int)((j+min.x-1+omega->center.x)*(double)SCALE);
//				cvCircle(img, data_pnt, 1, CV_RGB(0,255,0), 2, 8, 0);
			}
//			else if(omega->interior[i][j] == INTERNAL_AREA_IN){
//				data_pnt.y = (int)((i+min.y-1+omega->center.y)*(double)SCALE);
//				data_pnt.x = (int)((j+min.x-1+omega->center.x)*(double)SCALE);
//				cvCircle(img, data_pnt, 1, CV_RGB(255,0,0), 2, 8, 0);
//			}
		}
	}
	return;
}

/******************************************************************************/
// check_crossection : check the crossection of two lines g1-g2 and g3-g4
// return : 1(intersected), 0(not intersected)
/******************************************************************************/
int check_crossection(DPoint g1, DPoint g2, DPoint g3, DPoint g4, int neighbor)
{
	double a, b, c, d, e, f, tmp;
	a = g1.x - g2.x;
	b = g3.x - g4.x;
	c = g1.y - g2.y;
	d = g3.y - g4.y;
	e = g3.x - g2.x;
	f = g3.y - g2.y;
	if(neighbor){ // g2 = g3
		if((a*d - b*c) == 0)
			if(a==0)
				if(b==0)
					if((g2.y-g1.y)*(g3.y-g4.y)<=0) 
						return 0;
					else
						return 1;
				else // (a==0)&&(c==0) : g1=g2=g3
					return 0;
			else if(d==0)
				if(c==0)
					if((g2.x-g1.x)*(g3.x-g4.x)<=0) 
						return 0;
					else
						return 1;
				else // (d==0)&&(b==0) : g2=g3=g4
					return 0;
			else
				if((g2.x-g1.x)*(g3.x-g4.x)<=0) 
					return 0;
				else
					return 1;
		else
			return 0;
	}
	else{
		if (  ((g1.x==g4.x)&&(g1.y==g4.y))
			||((g1.x==g3.x)&&(g1.y==g3.y))
			||((g2.x==g4.x)&&(g2.y==g4.y))
			||((g2.x==g3.x)&&(g2.y==g3.y)))
			return 1;
		else{
			if (a==0){	// g1.x = g2.x
				if (c==0)	// g1.y = g2.y
					if (b==0)	// g3.x. = g4.x
						if (d==0)	// g3.y = g4.y
							return 0;
						else // (d!=0)		   g1 = g2, g3.x = g4.x, g3.y != g4.y (vertical line)
							if (e == 0){	// g2.x = g3.x
								tmp = f/d;
								if ((tmp >= 0)&&(tmp <= 1)) // point g1 = g2 is on line g3-g4
									return 1;
								else
									return 0;
							}
							else // (e!= 0) g2.x != g3.x
								return 0;
					else // (b!= 0) g1 = g2, g3.x != g4.x
						if (d==0)	// g3.y = g4.y (horizontal line)
							if (f==0){	// g2.y (= g1.y) = g3.y (= g4.y)
								tmp = e/b;
								if ((tmp >= 0)&&(tmp <= 1)) // point g1 = g2 is on line g3-g4
									return 1; // ok
								else
									return 0;
							}
							else // (f!=0)	   g2.y != g3.y
								return 0;
						else{ // (d!=0)		   g1 = g2, g3 != g4
							tmp = e/b;
							if (f/d == tmp)
								if ((tmp >= 0)&&(tmp <= 1)) // point g1 = g2 is on line g3-g4
									return 1; // ok
								else	// point g1 = g2 is on the same direction of line g3-g4 but not on the line
									return 0;
							else		// point g1 = g2 is not on the line g3-g4
								return 0; // ok
						}
				else{ // (c!=0) g1.y ! = g2.y, line g1-g2 is vertical line
					if (b==0)	// g3.x = g4.x
						if (d==0)	// g3.y = g4.y (g3 = g4)
							if (e==0){	// g2.x (= g1.x) = g3.x (= g4.x)
								tmp = f/c;
								if ((tmp >= 0)&&(tmp <= 1)) // point g3 = g4 is on line g1-g2
									return 1; // ok
								else	// point g3 = g4 is on the same direction of line g1-g2 but not on the line
									return 0; // ok
							}
							else // (e!=0)	g2.x != g3.x point g3 = g4 is not on the line g1-g2
								return 0;
						else // (d!=0)		g3-g4 (vertical line)
							if (e == 0){	// g2.x = g3.x
								if (g2.y > g1.y)	// g1--------g2
									if (g4.y > g3.y) // g3--------g4
										if (g3.y>g2.y) // g1--------g2   g3--------g4
											return 0;
										else
											if (g1.y>g4.y) //  g3--------g4    g1--------g2 
												return 0;
											else	//  g3-----g1---g4-------g2 : line g1-g2 and g3-g4 are overlapped
												return 1;
									else	// g4--------g3
										if (g4.y>g2.y) // g1--------g2   g4--------g3
											return 0;
										else
											if (g1.y>g3.y) //  g4--------g3    g1--------g2 
												return 0;
											else	//  g4-----g1---g3-------g2 : line g1-g2 and g3-g4 are overlapped
												return 1;
								else	// g2--------g1
									if (g4.y > g3.y) // g3--------g4
										if (g3.y>g1.y) // g2--------g1   g3--------g4
											return 0;
										else
											if (g2.y>g4.y) //  g3--------g4    g2--------g1 
												return 0;
											else	//  g3-----g2---g4-------g1 : line g1-g2 and g3-g4 are overlapped
												return 1;
									else	// g4--------g3
										if (g4.y>g1.y) // g2--------g1   g4--------g3
											return 0;
										else
											if (g2.y>g3.y) //  g4--------g3    g2--------g1 
												return 0;
											else	//  g4-----g2---g3-------g1 : line g1-g2 and g3-g4 are overlapped
												return 1;
							}
							else	// (e!=0)  g2.x != g3.x  two lines' direction are same but parallel
								return 0;
					else{	// (b!=0)	g3.x != g4.x
						tmp = (b*f-d*e)/(b*c);
						if ((tmp >= 0)&&(tmp <= 1)){	// cross section point can be on line g1-g2
							tmp = e/b;
							if ((tmp >= 0)&&(tmp <= 1)) // cross section point can be on line g3-g4
								return 1;
							else
								return 0;
						}
						else  // cross section point is not on line g1-g2
							return 0;
					}
				}
			}
			else // (a!=0)  g1.x != g2.x
				if (c==0)	// g1.y = g2.y (g1-g2 horizontal line)
					if (b==0)	// g3.x. = g4.x
						if (d==0)	// g3.y = g4.y (g3 and g4 is same point)
							if (f==0){	// g2.y = g3.y (= g4.x)
								tmp = e/a;
								if ((tmp >= 0)&&(tmp <= 1)) // point g3 = g4 is on line g1-g2
									return 1;
								else
									return 0; // ok
							}
							else // (f!=0)	   point g3 = g4 is not on the same direction of line g1-g2
								return 0; // ok
						else{ // (d!=0)	g3.y != g4.y (vertical line)
							tmp = e/a;
							if ((tmp >= 0)&&(tmp <= 1)){  // cross section point can be on line g1-g2
								tmp = f/d;
								if ((tmp >= 0)&&(tmp <= 1)) // cross section point can be on line g3-g4
									return 1;	// ok
								else
									return 0;				// cross section point can not be on line g3-g4
							}
							else
								return 0;					// cross section point can not be on line g1-g2
						}
					else // (b!= 0) g3.x != g4.x
						if (d==0)	// g3.y = g4.y (horizontal line)
							if (f==0){	// g2.y (= g1.y) = g3.y (= g4.y) (g1-g2 and g3-g4 are same horizontal direction)
								if (g2.x > g1.x)	// g1--------g2
									if (g4.x > g3.x) // g3--------g4
										if (g3.x>g2.x) // g1--------g2   g3--------g4
											return 0;
										else
											if (g1.x>g4.x) //  g3--------g4    g1--------g2 
												return 0;
											else	//  g3-----g1---g4-------g2 : line g1-g2 and g3-g4 are overlapped
												return 1;
									else	// g4--------g3
										if (g4.x>g2.x) // g1--------g2   g4--------g3
											return 0;
										else
											if (g1.x>g3.x) //  g4--------g3    g1--------g2 
												return 0;
											else	//  g4-----g1---g3-------g2 : line g1-g2 and g3-g4 are overlapped
												return 1;
								else	// g2--------g1
									if (g4.x > g3.x) // g3--------g4
										if (g3.x>g1.x) // g2--------g1   g3--------g4
											return 0;
										else
											if (g2.x>g4.x) //  g3--------g4    g2--------g1 
												return 0;
											else	//  g3-----g2---g4-------g1 : line g1-g2 and g3-g4 are overlapped
												return 1;
									else	// g4--------g3
										if (g4.x>g1.x) // g2--------g1   g4--------g3
											return 0;
										else
											if (g2.x>g3.x) //  g4--------g3    g2--------g1 
												return 0;
											else	//  g4-----g2---g3-------g1 : line g1-g2 and g3-g4 are overlapped
												return 1;
							}
							else // (f!=0)	   g2.y != g3.y (g1-g2 and g3-g4 are same horizontal direction but parallel)
								return 0;
						else{	// (d!=0)	g3.y != g4.y
							tmp = (d*e-b*f)/(a*d);
							if ((tmp >= 0)&&(tmp <= 1)){	// cross section point can be on line g1-g2
								tmp = f/d;
								if ((tmp >= 0)&&(tmp <= 1)) // cross section point can be on line g3-g4
									return 1;
								else
									return 0;
							}
							else  // cross section point is not on line g1-g2
								return 0;
						}
				else{ // (c!=0) g1.y ! = g2.y, line g1-g2 has incline (0<m<infinite)
					if (b==0)	// g3.x = g4.x
						if (d==0){	// g3.y = g4.y (g3=g4 same point)
							tmp = e/a;
							if (f/c == tmp)		// g3=g4 is on the direction of g1-g2
								if ((tmp >= 0)&&(tmp <= 1)) // point g3=g4 is on line g1-g2
									return 1; // ok
								else	// point g3 = g4 is on the same direction of line g1-g2 but not on the line
									return 0; // ok
							else		// point g3 = g4 is not on the line g1-g2
								return 0;
						}
						else{	// (d!=0)	g3-g4 (vertical line)
							tmp = e/a;
							if ((tmp >= 0)&&(tmp <= 1)){	// cross section point can be on line g1-g2
								tmp = (a*f-c*e)/(a*d);
								if ((tmp >= 0)&&(tmp <= 1)) // cross section point can be on line g3-g4
									return 1;
								else
									return 0;
							}
							else  // cross section point is not on line g1-g2
								return 0;
						}
					else{	// (b!=0)	g3.x != g4.x
						tmp = a*d-b*c;
						if (tmp == 0){
							if ((a*f-e*c) == 0){	// same direction
								if (g2.x > g1.x)	// g1--------g2
									if (g4.x > g3.x) // g3--------g4
										if (g3.x>g2.x) // g1--------g2   g3--------g4
											return 0; // ok
										else
											if (g1.x>g4.x) //  g3--------g4    g1--------g2 
												return 0;
											else	//  g3-----g1---g4-------g2 : line g1-g2 and g3-g4 are overlapped
												return 1; // ok
									else	// g4--------g3
										if (g4.x>g2.x) // g1--------g2   g4--------g3
											return 0;
										else
											if (g1.x>g3.x) //  g4--------g3    g1--------g2 
												return 0;
											else	//  g4-----g1---g3-------g2 : line g1-g2 and g3-g4 are overlapped
												return 1; // ok
								else	// g2--------g1
									if (g4.x > g3.x) // g3--------g4
										if (g3.x>g1.x) // g2--------g1   g3--------g4
											return 0; // ok
										else
											if (g2.x>g4.x) //  g3--------g4    g2--------g1 
												return 0;
											else	//  g3-----g2---g4-------g1 : line g1-g2 and g3-g4 are overlapped
												return 1; // ok
									else	// g4--------g3
										if (g4.x>g1.x) // g2--------g1   g4--------g3
											return 0;
										else
											if (g2.x>g3.x) //  g4--------g3    g2--------g1 
												return 0;
											else	//  g4-----g2---g3-------g1 : line g1-g2 and g3-g4 are overlapped
												return 1; // ok
							}
							else
								return 0;	// ok
						}
						else{
							tmp = (d*e-b*f)/(a*d-b*c);
							if ((tmp >= 0)&&(tmp <= 1)){	// cross section point can be on line g1-g2
								tmp = (a*f-c*e)/(a*d-b*c);
								if ((tmp >= 0)&&(tmp <= 1)) // cross section point can be on line g3-g4
									return 1; // ok
								else
									return 0; // ok
							}
							else  // cross section point is not on line g1-g2
								return 0; // ok
						}
					}
				}
		}
	}
}

void draw_curve_box(Curve *omega, IplImage *image)
{
	int k;
	CvPoint data_pntA, data_pntB;
	CvScalar blue = CV_RGB(0,0,255);
	CvScalar green = CV_RGB(0,255,0);

	// Draw circle
	for (k = 0;k<POINT_NUM;k++){
		if(k!=0){
			data_pntA.x = (int)((omega->r[k-1].x+omega->center.x)*SCALE);
			data_pntA.y = (int)((omega->r[k-1].y+omega->center.y)*SCALE);
			data_pntB.x = (int)((omega->r[k].x+omega->center.x)*SCALE);
			data_pntB.y = (int)((omega->r[k].y+omega->center.y)*SCALE);
			cvLine(image, data_pntA, data_pntB, green,0,8,0);
		}
	}
	data_pntA.x = (int)((omega->r[k-1].x+omega->center.x)*SCALE);
	data_pntA.y = (int)((omega->r[k-1].y+omega->center.y)*SCALE);
	data_pntB.x = (int)((omega->r[0].x+omega->center.x)*SCALE);
	data_pntB.y = (int)((omega->r[0].y+omega->center.y)*SCALE);
	cvLine(image, data_pntA, data_pntB, green,0,8,0);
	// Draw box
	data_pntA.x = (int)((omega->min.x+omega->center.x)*SCALE);
	data_pntA.y = (int)((omega->min.y+omega->center.y)*SCALE);
	data_pntB.x = (int)((omega->max.x+omega->center.x)*SCALE);
	data_pntB.y = (int)((omega->min.y+omega->center.y)*SCALE);
	cvLine(image, data_pntA, data_pntB, blue,0,8,0);
	data_pntA = data_pntB;
	data_pntB.x = (int)((omega->max.x+omega->center.x)*SCALE);
	data_pntB.y = (int)((omega->max.y+omega->center.y)*SCALE);
	cvLine(image, data_pntA, data_pntB, blue,0,8,0);
	data_pntA = data_pntB;
	data_pntB.x = (int)((omega->min.x+omega->center.x)*SCALE);
	data_pntB.y = (int)((omega->max.y+omega->center.y)*SCALE);
	cvLine(image, data_pntA, data_pntB, blue,0,8,0);
	data_pntA = data_pntB;
	data_pntB.x = (int)((omega->min.x+omega->center.x)*SCALE);
	data_pntB.y = (int)((omega->min.y+omega->center.y)*SCALE);
	cvLine(image, data_pntA, data_pntB, blue,0,8,0);
}


void draw_curve(Curve *omega, IplImage *image, int text_type, CvScalar color)
{
	int k;
	CvPoint data_pntA, data_pntB;
	// for text
	char text[2000];
	double hscale = 0.5;
	double vscale = 0.5;
	double shear = 0.2;
	int line_type = 8;
	CvFont font1;
	CvPoint pt1;
	CvScalar blue = CV_RGB(0,0,255);
	CvScalar red = CV_RGB(255,0,0);
	CvScalar yellow = CV_RGB(255,255,0);
	int thickness = 0;
	int connectivity = 8;

	cvInitFont(&font1,CV_FONT_HERSHEY_DUPLEX,hscale,vscale,shear,thickness,line_type);

	for(k=0;k<POINT_NUM;k++){
		if(k!=0){
			data_pntA.x = (int)((omega->r[k-1].x+omega->center.x)*SCALE);
			data_pntA.y = (int)((omega->r[k-1].y+omega->center.y)*SCALE);
			data_pntB.x = (int)((omega->r[k].x+omega->center.x)*SCALE);
			data_pntB.y = (int)((omega->r[k].y+omega->center.y)*SCALE);
			cvLine(image, data_pntA, data_pntB, color,0,8,0);
		}
	}
	data_pntA.x = (int)((omega->r[k-1].x+omega->center.x)*SCALE);
	data_pntA.y = (int)((omega->r[k-1].y+omega->center.y)*SCALE);
	data_pntB.x = (int)((omega->r[0].x+omega->center.x)*SCALE);
	data_pntB.y = (int)((omega->r[0].y+omega->center.y)*SCALE);
	cvLine(image, data_pntA, data_pntB, color,0,8,0);

	if(text_type != TEXT_NONE){
		switch(text_type)
		{
			case TEXT_SINGLE_E:
				sprintf(text, "%1.2f",omega->single_E);
				break;
			case TEXT_MULTIPLE_E:
				sprintf(text, "%1.2f",omega->multiple_E);
				break;
			case TEXT_TOTAL_E:
				sprintf(text, "%1.2f",omega->single_E + omega->multiple_E);
				break;
			case TEXT_BOTH_E:
				sprintf(text, "(%1.2f, %1.2f)",omega->single_E, omega->multiple_E);
				break;
			case TEXT_E0_E1_E2:
				sprintf(text, "(%1.2f: %1.2f, %1.2f, %1.2f)",omega->single_E,omega->e0, omega->e1, omega->e2);
				break;
			case TEXT_NUM_E0_E1_E2:
				hscale = 0.3;
				vscale = 0.3;
				cvInitFont(&font1,CV_FONT_HERSHEY_DUPLEX,hscale,vscale,shear,thickness,line_type);
				sprintf(text, "(%d: %1.0f, %1.0f, %1.0f)",omega->num,omega->e0, omega->e1, omega->e2);
				break;
			case TEXT_NUMBER:
			default:
				sprintf(text, "%d",omega->num);
				break;
		}
		pt1 = cvPoint((int)((omega->center.x+4)*SCALE), (int)(omega->center.y*SCALE));
		cvPutText(image,text,pt1,&font1,yellow);
	}
}

void draw_curve2(DPoint *new_r, DPoint *old_r, IplImage *image, CvScalar color)
{
	int k;
	CvPoint data_pntA, data_pntB;
	CvScalar red = CV_RGB(255,0,0);
	CvScalar green = CV_RGB(0,255,0);

	for(k=0;k<POINT_NUM;k++){
		if(k!=0){
			data_pntA.x = (int)((new_r[k-1].x*(double)SCALE));
			data_pntA.y = (int)((new_r[k-1].y*(double)SCALE));
			data_pntB.x = (int)((new_r[k].x  *(double)SCALE));
			data_pntB.y = (int)((new_r[k].y  *(double)SCALE));
			cvLine(image, data_pntA, data_pntB, color,0,8,0);
		}
		data_pntA.x = (int)((old_r[k].x  *(double)SCALE));
		data_pntA.y = (int)((old_r[k].y  *(double)SCALE));
		data_pntB.x = (int)((new_r[k].x  *(double)SCALE));
		data_pntB.y = (int)((new_r[k].y  *(double)SCALE));
		cvLine(image, data_pntA, data_pntB, green,0,8,0);
	}
	data_pntA.x = (int)((new_r[k-1].x*(double)SCALE));
	data_pntA.y = (int)((new_r[k-1].y*(double)SCALE));
	data_pntB.x = (int)((new_r[0].x  *(double)SCALE));
	data_pntB.y = (int)((new_r[0].y  *(double)SCALE));
	cvLine(image, data_pntA, data_pntB, color,0,8,0);
}

void draw_fill_curve(Curve *omega, IplImage *image)
{
	int i,j;
	CvPoint data_pntA;
	CvScalar red = CV_RGB(255,0,0);

	for (i=0;i < omega->max.y - omega->min.y + 3;i++){
		for (j=0;j < omega->max.x - omega->min.x + 3;j++){
			if(omega->interior[i][j] == INTERNAL_AREA_IN){
				data_pntA.y = (int)(i+omega->min.y-1+omega->center.y)*SCALE;
				data_pntA.x = (int)(j+omega->min.x-1+omega->center.x)*SCALE;
				cvCircle(image, data_pntA, 1*SCALE, red, 2, 8, 0);
			}
		}
	}
/*
	int k;
	CvPoint data_pnts[1][POINT_NUM];
    const CvPoint* ppt[1] = { data_pnts[0] };
	int npts[1];
	CvScalar red = CV_RGB(255,0,0);

	for (k= 0; k < POINT_NUM ; k++){
		data_pnts[0][k].x = (omega->r[k].x + omega->center.x);
		data_pnts[0][k].y = (omega->r[k].y + omega->center.y);
	}
	npts[0] = POINT_NUM;
	cvFillPoly( image, ppt, npts, 1, red,8,0);
*/

}


int check_twistofcurve(DPoint *new_r)
{
	int k, m, succeed=1;

	for (k = 0; k <= POINT_NUM-2 ; k++){
		if (k==0){
			succeed *= !check_crossection(new_r[0],new_r[1],new_r[1],new_r[2],1);
			if((succeed) == 0) return 0;
			for (m = 2; m <= POINT_NUM-2 ; m++){
				succeed *= !check_crossection(new_r[0],new_r[1],new_r[m],new_r[m+1],0);
			}
			succeed *= !check_crossection(new_r[POINT_NUM-1],new_r[0],new_r[0],new_r[1],1);
		}
		else if (k == POINT_NUM-2){
			succeed *= !check_crossection(new_r[POINT_NUM-2],new_r[POINT_NUM-1],new_r[POINT_NUM-1],new_r[0],1);
		}
		else{
			succeed *= !check_crossection(new_r[k],new_r[k+1],new_r[k+1],new_r[k+2],1);
			for (m = k+2; m <= POINT_NUM-1 ; m++){
				if (m == POINT_NUM-1){
					succeed *= !check_crossection(new_r[k],new_r[k+1],new_r[POINT_NUM-1],new_r[0],0);
				}
				else{
					succeed *= !check_crossection(new_r[k],new_r[k+1],new_r[m],new_r[m+1],0);
				}
			}
		}
	}
	return succeed;
}

/******************************************************************************/
// Qbirth : 
// 
/******************************************************************************/
void Qbirth(Curve *omega, DPoint center, double delta_t, int rmin, int rmax, int cols, int rows)
{
	int k;
	double radius;
	DPoint min, max, r;

	radius = (rmin + ((double)rand()/RAND_MAX)*(rmax-rmin));
	omega->state = STATE_NEW_BORN;
	omega->center = center;
	min.x = (double)cols;
	max.x = -(double)cols;
	min.y = (double)rows;
	max.y = -(double)rows;
	for (k = 0;k<POINT_NUM;k++){
//		r.x = (5*(delta_t*k-M_PI)*(delta_t*k-M_PI)/(M_PI*M_PI)+radius-5)*cos(delta_t*k); // to test smoothness constraint
		r.x = radius*cos(delta_t*k);
		if (r.x + center.x < 0.) r.x = - center.x;
		else if (r.x + center.x >= (double)cols-2) r.x = cols-2-center.x;
		if(r.x > (double)CURVE_MAX) r.x = (double)CURVE_MAX;
		if (r.x > max.x) max.x = r.x;
		if (r.x < min.x) min.x = r.x;
//		r.y = (5*(delta_t*k-M_PI)*(delta_t*k-M_PI)/(M_PI*M_PI)+radius-5)*sin(delta_t*k);
		r.y = radius*sin(delta_t*k);
		if (r.y + center.y < 0.) r.y = -center.y;
		else if (r.y + center.y >= (double)rows-2) r.y = (double)rows-2-center.y;
		if(r.y > (double)CURVE_MAX) r.y = (double)CURVE_MAX;
		if (r.y > max.y) max.y = r.y;
		if (r.y < min.y) min.y = r.y;
		omega->r[k] = r;
	}
	omega->min = min;
	omega->max = max;
}
/******************************************************************************/
// Single_Object : 
// 
/******************************************************************************/
int Single_Object(unsigned char **y, Curve *omega, double GD_step, double lambda_g, double lambda_G,
				  double GDepsilon,  double delta_t, double *mean, double *vari, int cols, int rows,
				  IplImage *image, const char* win_name, CvVideoWriter *video)
{
	int k,m,n,iter,succeed = 1;
	double d2r_d2tx, d2r_d2ty, Laplacian_I, dtmp;
	double G_r, G_hat_r, normalx, normaly, delta_Ex[POINT_NUM], delta_Ey[POINT_NUM];
	double delta_r = 100*GDepsilon, delta_r_old;
	DPoint old_r[POINT_NUM], new_r[POINT_NUM], normal_r[POINT_NUM];
	double dxp1, dxm1, dyp1, dym1, sum1, sum2, sum3;
	int i, j, sizex, sizey;
	double gradientx, gradienty;
	double delta_t2 = delta_t*delta_t;
	CvScalar red = CV_RGB(255,0,0);
	CvScalar green = CV_RGB(0,255,0);
	CvScalar blue = CV_RGB(0,0,255);
	CvScalar yellow = CV_RGB(255,255,0);
	CvScalar cyan = CV_RGB(0,255,255);
	CvScalar white = CV_RGB(255,255,255);
//	unsigned char **laplacian;

	// Draw Laplacian
//	laplacian = (unsigned char **)get_img(cols, rows, sizeof(unsigned char));
//	for (i = 1; i < rows-1; i++)
//		for (j = 1; j < cols-1; j++){
//			tmp = y[i+1][j]+y[i-1][j]+y[i][j+1]+y[i][j-1]-4*y[i][j];
//			if(tmp==0)
//				laplacian[i][j] = 0;
//			else
//				laplacian[i][j] = tmp+68;
//		}

	iter = 0;
	for (k= 0; k < POINT_NUM ; k++){
		new_r[k].x = (omega->r[k].x + omega->center.x);
		new_r[k].y = (omega->r[k].y + omega->center.y);
	}
	// Do Gradient Descent Optimization to find minimum E
	do {
		delta_r_old = delta_r;
		for (k= 0; k < POINT_NUM ; k++){
			old_r[k].x = new_r[k].x;
			old_r[k].y = new_r[k].y;
		}
		for (k = 0; k < POINT_NUM ; k++){
			if (k==0){
				d2r_d2tx = old_r[POINT_NUM-1].x;
				d2r_d2tx = (old_r[k+1].x-2*old_r[k].x+old_r[POINT_NUM-1].x)/delta_t2;
				d2r_d2ty = (old_r[k+1].y-2*old_r[k].y+old_r[POINT_NUM-1].y)/delta_t2;
				normalx =  (old_r[k+1].y-old_r[k].y);
				normaly = -(old_r[k+1].x-old_r[k].x);
			}
			else if(k==POINT_NUM-1){
				d2r_d2tx = (old_r[0].x-2*old_r[k].x+old_r[k-1].x)/delta_t2;
				d2r_d2ty = (old_r[0].y-2*old_r[k].y+old_r[k-1].y)/delta_t2;
				normalx =  (old_r[0].y-old_r[k].y);
				normaly = -(old_r[0].x-old_r[k].x);
			}
			else{
				d2r_d2tx = (old_r[k+1].x-2*old_r[k].x+old_r[k-1].x)/delta_t2;
				d2r_d2ty = (old_r[k+1].y-2*old_r[k].y+old_r[k-1].y)/delta_t2;
				normalx =  (old_r[k+1].y-old_r[k].y);
				normaly = -(old_r[k+1].x-old_r[k].x);
			}
//			dtmp = sqrt(normalx*normalx+normaly*normaly);
//			normalx = normalx/dtmp;
//			normaly = normaly/dtmp;
			dxp1 = old_r[k].x + 1.;
			if(dxp1 > (double)(cols-2)) dxp1 = (double)(cols-2);
			dxm1 = old_r[k].x - 1.;
			if(dxm1 < 0) dxm1 = 0.;
			dyp1 = old_r[k].y + 1.;
			if(dyp1 > (double)(rows-2)) dyp1 = (double)(rows-2);
			dym1 = old_r[k].y - 1.;
			if(dym1 < 0) dym1 = 0.;
			Laplacian_I = real_coord(y, old_r[k].y, dxp1) + real_coord(y, old_r[k].y, dxm1)
						+ real_coord(y, dyp1, old_r[k].x) + real_coord(y, dym1, old_r[k].x)
						- 4*real_coord(y, old_r[k].y, old_r[k].x);
			//Laplacian_I = (double)y[(int)(old_r[k].y+0.5)][(int)(dxp1+0.5)] + (double)y[(int)(old_r[k].y+0.5)][(int)(dxm1+0.5)] 
			//			+ (double)y[(int)(dyp1+0.5)][(int)(old_r[k].x+0.5)] + (double)y[(int)(dym1+0.5)][(int)(old_r[k].x+0.5)]
			//			- 4*(double)y[(int)(old_r[k].y+0.5)][(int)(old_r[k].x+0.5)];
			G_r =  real_coord(y, old_r[k].y, old_r[k].x) - mean[1];
			G_r = G_r*G_r;
			G_r = G_r/(2.*vari[1]);
			G_hat_r = real_coord(y, old_r[k].y, old_r[k].x) - mean[0];
			G_hat_r = G_hat_r*G_hat_r;
			G_hat_r = G_hat_r/(2.*vari[0]);
			dtmp = lambda_g*Laplacian_I+lambda_G*(G_r-G_hat_r);
			delta_Ex[k] = -2.*d2r_d2tx+normalx*dtmp;
			delta_Ey[k] = -2.*d2r_d2ty+normaly*dtmp;
			new_r[k].x = old_r[k].x - (GD_step*delta_Ex[k]);
			if(new_r[k].x < 0) new_r[k].x = 0.;
			else if(new_r[k].x > (double)(cols-2)) new_r[k].x = (double)(cols-2);
			new_r[k].y = old_r[k].y - (GD_step*delta_Ey[k]);
			if(new_r[k].y < 0) new_r[k].y = 0.;
			else if(new_r[k].y > (double)(rows-2)) new_r[k].y = (double)(rows-2);
			
		}

		if(SHOW_WINDOW*LEVEL4){
			// Clear and load image
			LoadImageFromMemory(image, y);
			draw_curve2(new_r, old_r, image, red);
			// show
			cvShowImage(win_name, image);
			if (iter%10 ==0) cvWriteFrame( video, image );
//			printf("iter = %d\n",iter);
			// wait for a key
			cvWaitKey((int)(1));
			//cvWaitKey((int)(0));
		}

	// Check cross-section. If happen fail
		succeed = check_twistofcurve(new_r);

		// || r(k+1) - rk) ||
		delta_r = 0.;
		for (k= 0; k < POINT_NUM ; k++){
			dtmp = (new_r[k].x - old_r[k].x)*(new_r[k].x - old_r[k].x)
					+(new_r[k].y - old_r[k].y)*(new_r[k].y - old_r[k].y);
			delta_r += sqrt(dtmp);
		}
		iter++;
	} while((fabs(delta_r-delta_r_old) > GDepsilon)&&(iter<GD_ITER_MAX));

//	succeed = check_twistofcurve(new_r);

	if (succeed){
		// update curve
		sum1 = 0.;
		sum2 = 0.;
		for (k= 0; k < POINT_NUM ; k++){
			sum1 += new_r[k].x;
			sum2 += new_r[k].y;
		}
		sum1 /= (double)POINT_NUM;
		sum2 /= (double)POINT_NUM;
		omega->center.x = sum1;
		omega->center.y = sum2;
		for (k= 0; k < POINT_NUM ; k++){
			omega->r[k].x = (new_r[k].x - omega->center.x);
			omega->r[k].y = (new_r[k].y - omega->center.y);
			if(omega->min.x > omega->r[k].x)
				omega->min.x = omega->r[k].x;
			if(omega->max.x < omega->r[k].x)
				omega->max.x = omega->r[k].x;
			if(omega->min.y > omega->r[k].y)
				omega->min.y = omega->r[k].y;
			if(omega->max.y < omega->r[k].y)
				omega->max.y = omega->r[k].y;
		}

		//  calculate Energy
		Internal_Area(omega);

		if(SHOW_WINDOW*LEVEL5){
			draw_fill_curve(omega, image);		
			cvShowImage(win_name, image);
			printf("pause\n");
			// wait for a key
			cvWaitKey(0);
		}
		
		sum1 = 0.;
		sum2 = 0.;
		sum3 = 0.;
		for (k = 0; k < POINT_NUM ; k++){
			sum1 += dist(omega->center,new_r[k]);
		}
		sum1 = sum1/POINT_NUM;
		for (k = 0; k < POINT_NUM ; k++){
			normal_r[k].x = (new_r[k].x - omega->center.x)/sum1;
			normal_r[k].y = (new_r[k].y - omega->center.y)/sum1;
		}
		sum1 = sum1/POINT_NUM;

		sum1 = 0.;
		sum2 = 0.;
		sum3 = 0.;
		for (k = 0; k < POINT_NUM ; k++){
			if(k==POINT_NUM-1){
				d2r_d2tx = (normal_r[0].x-normal_r[k].x)/delta_t;
				d2r_d2ty = (normal_r[0].y-normal_r[k].y)/delta_t;
			}
			else{
				d2r_d2tx = (normal_r[k+1].x-normal_r[k].x)/delta_t;
				d2r_d2ty = (normal_r[k+1].y-normal_r[k].y)/delta_t;
			}
			normalx =  d2r_d2ty;
			normaly = -d2r_d2tx;
			dtmp = sqrt(normalx*normalx+normaly*normaly);
			normalx = normalx/dtmp;
			normaly = normaly/dtmp;
			// E_curve
			sum1 += (d2r_d2tx*d2r_d2tx+d2r_d2ty*d2r_d2ty)*delta_t;
			// E_grad
			dxp1 = new_r[k].x + 1.;
			if(dxp1 > (double)(cols-2)) dxp1 = (double)(cols-2);
			dyp1 = new_r[k].y + 1.;
			if(dyp1 > (double)(rows-2)) dyp1 = (double)(rows-2);
			gradientx = real_coord(y, new_r[k].y, dxp1) - real_coord(y, new_r[k].y, new_r[k].x);
			gradienty = real_coord(y, dyp1, new_r[k].x) - real_coord(y, new_r[k].y, new_r[k].x);
			sum2 += (lambda_g*(normalx*gradientx+normaly*gradienty))*delta_t;
		}
		
		sizex = (int)(omega->max.x - omega->min.x + 3);
		sizey = (int)(omega->max.y - omega->min.y + 3);
		
		// E_gauss
		/*
		for (k = 0; k < POINT_NUM ; k++){
			dxp1 = omega->r[k].x + omega->center.x;
			dyp1 = omega->r[k].y + omega->center.y;
			G_r =  real_coord(y, dyp1, dxp1) - mean[1];
			G_r = G_r*G_r;
			G_r = G_r/(2.*vari[1]);
			G_hat_r = real_coord(y, dyp1, dxp1) - mean[0];
			G_hat_r = G_hat_r*G_hat_r;
			G_hat_r = G_hat_r/(2.*vari[0]);
			sum1 += lambda_G*(G_r-G_hat_r);
		}
		*/
		n = 0;
		for (k = 0; k < sizey ; k++){
			for (m = 0; m < sizex ; m++){
				if (omega->interior[k][m] == INTERNAL_AREA_IN){
					i = k+(int)(omega->center.y + omega->min.y) - 1;
					j = m+(int)(omega->center.x + omega->min.x) - 1;
					G_r =  real_coord(y, i, j) - mean[1];
					G_r = G_r*G_r;
					G_r = G_r/(2.*vari[1]);
					G_hat_r = real_coord(y, i, j) - mean[0];
					G_hat_r = G_hat_r*G_hat_r;
					G_hat_r = G_hat_r/(2.*vari[0]);
					sum3 += lambda_G*(G_r-G_hat_r);
					n++;
				}
			}
		}
		sum3 = sum3/(double)n;
		omega->e0 = sum1;
		omega->e1 = sum2;
		omega->e2 = sum3;
//		printf("energy2 = %1.3f\n",sum2);
		sum2 = sum2;///40.;
		sum3 = sum3;///20.;
		omega->single_E = (sum1+sum2+sum3)/1000;//50.;
/* Another method to calculate Energy of curve ???
		sum1 = 0.;
		for (k = 0; k < POINT_NUM ; k++){
			if (k==POINT_NUM-1){
				sum1 = (new_r[0].x - new_r[k].x)*delta_Ex[k]+(new_r[0].y - new_r[k].y)*delta_Ey[k];
			}
			else{
				sum1 = (new_r[k+1].x - new_r[k].x)*delta_Ex[k]+(new_r[k+1].y - new_r[k].y)*delta_Ey[k];
			}
		}
		omega->single_E = sum1;
*/
	}
	else{
		omega->single_E = 0;
	}
	if(SHOW_WINDOW*LEVEL3){
		draw_curve2(new_r, old_r, image, blue);
		cvShowImage(win_name, image);
		cvWriteFrame( video, image );
		// wait for a key
		cvWaitKey(1);
	}

//	free_img((void **)laplacian);

	return succeed;
}


/******************************************************************************/
// test_single_object : 
/******************************************************************************/
void test_single_object(unsigned char **yimg, double GD_step,  
						double lambda_g, double lambda_G, double GDepsilon,
						double *mean, double *vari, int cols, int rows, 
						IplImage *image, const char* win_name, CvVideoWriter *video)
{
	double test_pointx[15] = {6,-39,-41,-33,-23,-23,-57,-49,-28,1,0,49,50,13,60};
	double test_pointy[15] = {-59,-38,-15,-23,-23,10,14,51,36,48,77,82,22,-13,-29};
	double test_centerx = 63;
	double test_centery = 97;
	DPoint center;
//	CvPoint data_pntA, data_pntB;
	int i;
	int succeed;
	Curve omega;
	double delta_t = DELTA_T,r;
	char ch;
	CvScalar red = CV_RGB(255,0,0);
	CvSize frame_size;

#ifdef TEST_POLY
	center.x = test_centerx;
	center.y = test_centery;
#else
	center.x = 115; 
	center.y = 70;	
#endif
	r = 2;
	omega.state = STATE_NEW_BORN;
	omega.single_E = 0;
	omega.e0 = 0;
	omega.e1 = 0;
	omega.e2 = 0;
	Qbirth(&omega, center, delta_t, r, r, cols, rows);
	if(SHOW_WINDOW*LEVEL1){
		LoadImageFromMemory(image, yimg);
		draw_curve(&omega, image, TEXT_NONE, red);
		/* write frame to video file */
		cvWriteFrame( video, image );
		// create a window
		cvShowImage(win_name, image);
		cvWaitKey(0);
	}
	printf("waiting key..... \n");
	while((ch=getch())!='q'){   
		putchar(ch);
		switch(ch)
		{
			case '8':
				center.y--;
				printf("\ny=%d\n",(int)center.y);
				Qbirth(&omega, center, delta_t, r, r, cols, rows);
				break;
			case '2':
				center.y++;
				printf("\ny=%d\n",(int)center.y);
				Qbirth(&omega, center, delta_t, r, r, cols, rows);
				break;
			case '4':
				center.x--;
				printf("\nx=%d\n",(int)center.x);
				Qbirth(&omega, center, delta_t, r, r, cols, rows);
				break;
			case '6':
				center.x++;
				printf("\nx=%d\n",(int)center.x);
				Qbirth(&omega, center, delta_t, r, r, cols, rows);
				break;
			case '7':
				r--;
				printf("\nr=%1.2f\n",r);
				Qbirth(&omega, center, delta_t, r, r, cols, rows);
				break;
			case '9':
				r++;
				printf("\nr=%1.2f\n",r);
				Qbirth(&omega, center, delta_t, r, r, cols, rows);
				break;
			case '5':
				printf("\n(%d,%d,%1.2f)\n", (int)center.x, (int)center.y, r);
				succeed = Single_Object(yimg, &omega, GD_step, lambda_g,lambda_G, GDepsilon,
								delta_t, mean, vari, cols, rows, image, win_name, video );
				if(!succeed) printf("GD Fail!\n");
				break;
			case '0':
				Qbirth(&omega, center, delta_t, r, r, cols, rows);
				break;
			default:
				break;
		}
		if(SHOW_WINDOW*LEVEL1){
			LoadImageFromMemory(image, yimg);
			draw_curve(&omega, image, TEXT_E0_E1_E2, red);
			cvWriteFrame( video, image );
			// create a window
			cvShowImage(win_name, image);
			cvWaitKey((int)(1));
		}
	}
	/* close video file writer */
	printf("GD..... \n");
	// GD optimization
#if 0
	printf("(C,%1.5f,%1.5f) ",omega.center.x,omega.center.y);
	printf("(MN,%1.5f,%1.5f) ",omega.min.x,omega.min.y);
	printf("(MX,%1.5f,%1.5f)\n",omega.max.x,omega.max.y);
	for(i=0;i<POINT_NUM;i++){
		printf("(%d,%1.5f,%1.5f) ",i,omega.r[i].x,omega.r[i].y);
	}
	printf("\n");
#endif
#if 0
	printf("(C,%1.5f,%1.5f) ",omega.center.x,omega.center.y);
	printf("(MN,%1.5f,%1.5f) ",omega.min.x,omega.min.y);
	printf("(MX,%1.5f,%1.5f)\n",omega.max.x,omega.max.y);
	for(i=0;i<POINT_NUM;i++){
		printf("(%d,%1.5f,%1.5f) ",i,omega.r[i].x,omega.r[i].y);
	}
	printf("\n");
#endif

	// Calculate internal Area test
	/*
	if(succeed){
		Internal_Area(&omega);
		draw_fill_curve(&omega, image);
		if(SHOW_WINDOW*LEVEL1){
			// create a window
			cvShowImage(win_name, image);
			printf("pause\n");
			// wait for a key
			cvWaitKey((int)((100-i)/5.));
		}
	}
	*/
}

#define TEST_OBJ_NUM	26
/******************************************************************************/
// test_single_object : 
/******************************************************************************/
void test_single_object2(unsigned char **yimg, double GD_step,  
						double lambda_g, double lambda_G, double GDepsilon,
						double *mean, double *vari, int cols, int rows, 
						IplImage *image, const char* win_name, CvVideoWriter *video)
{
	double test_circle[TEST_OBJ_NUM][4] = {	{89,  39,  9, -4},	// 0o
											{115, 90, 23, -4},	// 1o
											{140, 70,  4, -4},	// 2o
											{149, 54,  9, -4},	// 3o
											{159, 76,  9, -4},	// 4o
											{ 67, 28,  5, -4},	// 5o
											{ 42, 46,  4, -4},	// 6o
											{ 32, 38,  4, -4},	// 7o
											{ 19,106, 10, -4},	// 8o
											{ 36,100,  6, -4},	// 9o
											{ 65, 75,  6, -4},	// 10o
											{ 77, 85,  6, -4},	// 11o
											{ 63,102, 11, -4},	// 12o
											{113, 68,  9, -3},	// 13x falsly disconnected
											{121,104, 22, -3},	// 14x 
											{ 71, 79, 11, -2},	// 15xx merge 2 objects
											{ 37, 41,  7, -2},	// 16xx
											{142, 59, 14, -2},	// 17xx
											{ 85, 29, 13, -2},	// 18xx
											{ 24, 41,  9, -2},	// 19xx
											{ 21,105, 15, -1},	// 20xxx merge 3 objects
											{ 81, 31, 19, -1},	// 21xxx
											{148, 65, 19, -1},	// 22xxx
											{113, 67, 12, -1},	// 23xxx stop during GD
											{138, 79, 22, -1},	// 24xxxx merge many objects
											{ 63, 79, 22,  0}};	// 25xxxx
	DPoint center;
//	CvPoint data_pntA, data_pntB;
	int i;
	int succeed;
	Curve omega[TEST_OBJ_NUM];
	double delta_t = DELTA_T,r;
	char ch;
	CvScalar red = CV_RGB(255,0,0);
	CvSize frame_size;
	char outfileName[100];
	FILE *fp;

	for(i=0;i<TEST_OBJ_NUM;i++){
		omega[i].state = STATE_NEW_BORN;
		omega[i].num = i;
		omega[i].single_E = 0;
		omega[i].e0 = 0;
		omega[i].e1 = 0;
		omega[i].e2 = 0;
		center.x	= test_circle[i][0];
		center.y	= test_circle[i][1];
		r			= test_circle[i][2];

		Qbirth(&(omega[i]), center, delta_t, r, r, cols, rows);
	}
	if(SHOW_WINDOW*LEVEL1){
		LoadImageFromMemory(image, yimg);
		draw_all_curves(omega, TEST_OBJ_NUM, image, TEXT_NONE);
		/* write frame to video file */
		cvWriteFrame( video, image );
		// create a window
		cvShowImage(win_name, image);
		cvWaitKey(0);
	}
	for(i=0;i<TEST_OBJ_NUM;i++){
		succeed = Single_Object(yimg, &(omega[i]), GD_step, lambda_g,lambda_G, GDepsilon,
						delta_t, mean, vari, cols, rows, image, win_name, video );
		if(!succeed) printf("%d item GD Fail!\n",i);
	}
	if(SHOW_WINDOW*LEVEL1){
		LoadImageFromMemory(image, yimg);
		draw_all_curves(omega, TEST_OBJ_NUM, image, TEXT_NUM_E0_E1_E2);
		cvWriteFrame( video, image );
		// create a window
		cvShowImage(win_name, image);
		cvWaitKey((int)(0));
	}
	sprintf(outfileName, "single_obj_test.txt");
	if ((fp = fopen(outfileName, "wb")) == NULL ) {
		printf("Cannot open file %s\n", outfileName);
		exit(1);
	}
	fprintf(fp,"smoothness, gradient, gaussian, s_E hope\r\n");
	for(i=0;i<TEST_OBJ_NUM;i++){
		fprintf(fp,"%1.4f, %1.4f, %1.4f, %1.4f\r\n", omega[i].e0, omega[i].e1, omega[i].e2, test_circle[i][3]);
	}
	fclose(fp);
}

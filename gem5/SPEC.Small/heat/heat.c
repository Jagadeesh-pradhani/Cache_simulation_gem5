
#include<stdio.h>
#include<math.h>
/*
	 Aryabartta Sahu, Dept of Computer Sc and Engg 
	 asahu@cse.iitd.ac.in, IIT Delhi 

 **/



#define f(x,y)     (sin(x)*sin(y))
#define randa(x,t) (0.0)
#define randb(x,t) (exp(-2*(t))*sin(x))
#define randc(y,t) (0.0)
#define randd(y,t) (exp(-2*(t))*sin(y))
#define solu(x,y,t) (exp(-2*(t))*sin(x)*sin(y))

#define XSIZE 32 
#define YSIZE 16

int nx, ny, nt;
float xu, xo, yu, yo, tu, to;
float dx, dy, dt;

float dtdxsq, dtdysq;
float t;

float old[XSIZE][YSIZE];
float new[XSIZE][YSIZE];

int leafmaxcol;

/*****************   Initialization of grid partition  ********************/

void initgrid(int lb, int ub){

				int a, b, llb, lub;
				float v1,v2;

				llb = (lb == 0) ? 1 : lb;
				lub = (ub == nx) ? nx - 1 : ub;
				b=0;
				for (a=llb; a < lub; a++)	
								old[a][b] = randa(xu + a * dx, 0);

				b=ny-1;
				for (a=llb; a < lub; a++)
								old[a][b] = randb(xu + a * dx, 0);

				if (lb == 0) {
								for (a=0, b=0; b < ny; b++)
												old[a][b] = randc(yu + b * dy, 0);
				}
				if (ub == nx) {
								b=0;
								for (a=nx-1; b < ny; b++)
												old[a][b] = randd(yu + b * dy, 0);
				}

				for (a=llb; a < lub; a++) {	
								for (b=1; b < ny-1; b++) {
												v1=xu + a * dx;
												v2=yu + b * dy;
												old[a][b] = f(v1,v2);
								}
				} 
}


/***************** Five-Point-Stencil Computation ********************/

void compstripeO2N(int lb, int ub)
{

				int a, b, llb, lub;

				llb = (lb == 0) ? 1 : lb;
				lub = (ub == nx) ? nx - 1 : ub;

				for (a=llb; a < lub; a++) {
								for (b=1; b < ny-1; b++) {

												new[a][b] =   dtdxsq * (old[a+1][b] - 2 * old[a][b] + old[a-1][b])
																+ dtdysq * (old[a][b+1] - 2 * old[a][b] + old[a][b-1])
																+ old[a][b];
								}
				}
				b=ny-1;
				for (a=llb; a < lub; a++)
								new[a][b] = randb(xu + a * dx, t);
				b=0;
				for (a=llb; a < lub; a++)
								new[a][b] = randa(xu + a * dx, t);

				if (lb == 0) {
								a=0;
								for (b=0; b < ny; b++)
												new[a][b] = randc(yu + b * dy, t);
				}
			
				if (ub == nx) {
								a=nx-1;
								for (b=0; b < ny; b++)
												new[a][b] = randd(yu + b * dy, t);
				} 
}


void compstripeN2O(int lb, int ub)
{

				int a, b, llb, lub;

				llb = (lb == 0) ? 1 : lb;
				lub = (ub == nx) ? nx - 1 : ub;

				for (a=llb; a < lub; a++) {
								for (b=1; b < ny-1; b++) {
												old[a][b] =   dtdxsq * (new[a+1][b] - 2 * new[a][b] + new[a-1][b])
																+ dtdysq * (new[a][b+1] - 2 * new[a][b] + new[a][b-1])
																+ new[a][b];
								}
				}
				b=ny-1;
				for (a=llb; a < lub; a++)
								old[a][b] = randb(xu + a * dx, t);
				b=0;
				for (a=llb; a < lub; a++)
								old[a][b] = randa(xu + a * dx, t);

				if (lb == 0) {
								a=0;
								for (b=0; b < ny; b++)
												old[a][b] = randc(yu + b * dy, t);
				}
				if (ub == nx) {
								a=nx-1;
								for (b=0; b < ny; b++)
												old[a][b] = randd(yu + b * dy, t);
				}
}



/***************** Decomposition of 2D grids in stripes ********************/

void heat() {
				int  c;
				initgrid( 0, nx); 
				for (c = 1; c <= nt; c++) {
								t = tu + c * dt;
								if ( c % 2 )
												compstripeO2N(0, nx);
								else 	compstripeN2O(0, nx);
				}
}


void init(){

				nx = XSIZE;
				ny = YSIZE;
				nt = 5;
				xu = 0.0;
				xo = 1.5707;
				yu = 0.0;
				yo = 1.570;
				tu = 0.0;
				to = 0.0001;
				leafmaxcol = 10;

				dx = (xo - xu) / (nx - 1);
				dy = (yo - yu) / (ny - 1);
				dt = (to - tu) / nt;
				dtdxsq = dt / (dx * dx);

}


int main(int argc, char *argv[])
{

				init();	
				heat();

				return 0;
}

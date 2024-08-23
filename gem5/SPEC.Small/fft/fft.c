
#include <stdio.h>
#include <math.h>

/* Sahu : Incsease SIZE to increase Execution time : Power of SMALLSIZE*/
/*This Program can run upto maxsize 256, Beause of the reverse() written in that way */
#define SIZE  2048
#define SMALLSIZE 16  /*it can be any number*/



#define PI 3.14159265358979323846
/*
double fabs(double x) {
				return x < 0 ? -x : x;
}
double sin(double rad) {
				double app, diff;
				int inc = 1;
				while( rad > 2*PI ) rad -= 2*PI;
				while( rad < -2*PI ) rad += 2*PI;
				app = diff = rad;
				diff = (diff * (-(rad*rad))) / ((2.0 * inc) * (2.0 * inc + 1.0));
				app = app + diff;
				inc++;
				while( fabs(diff) >= 0.00001 ) {
								diff = (diff * (-(rad*rad))) / ((2.0 * inc) * (2.0 * inc + 1.0));
								app = app + diff;
								inc++;
				}
				return app;
}

double cos(double rad) {
				double app,diff;
				int inc = 1;
				rad += PI/2.0;
				while( rad > 2*PI ) rad -= 2*PI;
				while( rad < -2*PI ) rad += 2*PI;
				app = diff = rad;
				diff = (diff * (-(rad*rad))) / ((2.0 * inc) * (2.0 * inc + 1.0));
				app = app + diff;
				inc++;
				while( fabs(diff) >= 0.00001 ) {
								diff = (diff * (-(rad*rad))) / ((2.0 * inc) * (2.0 * inc + 1.0));
								app = app + diff;
								inc++;
				}
				return app;
}
*/

typedef struct complex {
				float r,i;			/* Real and imaginery parts */
} COMPLEX;

COMPLEX TF[SIZE];
COMPLEX input[SIZE]; 
int new_index[SIZE];


int reverse(int input, int n) {
				int output=0;
				int i;
				int ip=input;
				for(i=0;i<n;i++){
								output=output*2+ip%2;
								ip=ip/2;
				}
				return output;
}

void par_k (int k,int M){

				int i1,i2,u;
				COMPLEX A,B,W2M,Half,tmp1, tmp2, tmp3, tmp4,tmp5;
				Half.r = 0.5;
				Half.i = 0.0;
				i1 = k*2*M;        
				i2 = (k*2+1)*M;    
				for (u=0; u<M; u++) {

								W2M.r = cos(PI*u/M);
								W2M.i = -sin(PI*u/M);

								tmp1=TF[i2+u];

								tmp2.r = tmp1.r*W2M.r - tmp1.i*W2M.i;
								tmp2.i = tmp1.r*W2M.i + W2M.r*tmp1.i;

								tmp3=TF[i1+u];
								tmp4.r=tmp3.r+tmp2.r; tmp4.i=tmp3.i+tmp2.i;
								tmp5.r=tmp3.r-tmp2.r; tmp5.i=tmp3.i-tmp2.i;

								A.r = Half.r*tmp4.r ;
								A.i = Half.r*tmp4.i ;
								B.r = Half.r*tmp5.r ;
								B.i = Half.r*tmp5.i ;

								TF[i1+u] = A;
								TF[i1+u+M] = B;
				}

}



/********************************************
 ** FFT calculates the 1D-FFT of a complex  **
 ** array passed into the function. The ar- **
 ** ray is expected to be of size N.        **
 ********************************************/
void FFT(COMPLEX F[], int N) {

				int i, M, j, k ,n=0, Ni;
				int i1,i2,u;
				COMPLEX A,B,W2M,Half,tmp1, tmp2, tmp3, tmp4,tmp5;
				Half.r = 0.5;
				Half.i = 0.0;
				/**Assumed N is 2^n log N calculation as follow**/	
				Ni=N;
				while(Ni!=0) {
								Ni=Ni/2;
								n++;
				}
				/*n=(int) log(N) / log(2);*/

				Half.r = 0.5;	Half.i = 0.0;

				for (i=0; i<N; i++) {
								new_index[i] = reverse(i,n);
				}
				/*****************************************
				 ** Rearrange array elements
				 *****************************************/
				for (i=0; i<N; i++) {
								TF[new_index[i]] = F[i];
				}
				M = 1;			/* The initial length of subgroups  */
				j = N/2;			/* The number of pairs of subgroups */
				/** Begin successive merges for n levels */

				for (i=0; i<n; i++) {
								if(M>SMALLSIZE) break;
								for (k=0; k<j; k++) {
												i1 = k*2*M;        
												i2 = (k*2+1)*M;    
												for (u=0; u<M; u++) {

																W2M.r = cos(PI*u/M);
																W2M.i = -sin(PI*u/M); 


																tmp1=TF[i2+u];
																tmp2.r = tmp1.r*W2M.r - tmp1.i*W2M.i;
																tmp2.i = tmp1.r*W2M.i + W2M.r*tmp1.i;
																tmp3=TF[i1+u];
																tmp4.r=tmp3.r+tmp2.r; tmp4.i=tmp3.i+tmp2.i;
																tmp5.r=tmp3.r-tmp2.r; tmp5.i=tmp3.i-tmp2.i;
																A.r = Half.r*tmp4.r ;
																A.i = Half.r*tmp4.i ;
																B.r = Half.r*tmp5.r ;
																B.i = Half.r*tmp5.i ;
																TF[i1+u] = A;
																TF[i1+u+M] = B;
												}
								}
								M = 2*M;
								j = j/2;
				}

				for(i=i; i<n; i++) {
								/*printf("\nPAR:i=%d M=%d j=%d\n", i,M,j);*/
								for(k=0; k<j; k++) {
												par_k(k,M);
								}

								M = 2*M;
								j = j/2;
				}
}


int main(){
				/*Init();*/
				FFT(input, SIZE);
				/*Output is stored in TF[SIZE]*/

				return 0;
}


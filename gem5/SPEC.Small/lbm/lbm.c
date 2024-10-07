/******************Poiseuille flow using LBM***********************
  Based on procedures as explained in 'Lattice Gas Cellular Automata and Lattice  Boltzmann
  Models'by Wolf Gladrow	
  */
#include <stdio.h>
#include <string.h>

float OMEGA=0.2; /*Relaxation factor*/
float TAU=0.5;

/* Sahu: Increase XMAX and YMAX to increase execution time*/

#define XMAX 64 /*Mesh size in x-direction  64*/
#define YMAX 64 /*Mesh size in y-direction  64*/
#define ZMAX 10 /*Not a-direction */
#define STEP 8 
#define MAXiter 200

float FORCING=(1.024)/YMAX*YMAX*YMAX;/*Forcing term in Poiseuille equation*/
float KVISC = ((1/0.2)-0.5)/3;
float UX;
int cx[9]={1, 0, -1, 0, 1, -1, -1, 1,0}; /*components of lattice velocities in x-direction*/
int cy[9]={0, 1, 0, -1, 1, 1, -1, -1, 0}; /*components of lattice velocities in y-direction*/
float jx[XMAX][YMAX];
float jy[XMAX][YMAX];
float rho[XMAX][YMAX];
float f[XMAX][YMAX][XMAX];
float fprop[XMAX][YMAX][XMAX];
float feq[XMAX][YMAX][XMAX];
float u[XMAX][YMAX],v[XMAX][YMAX];
float force;
float uprofile[YMAX];

void  JxJyRhoFFprop(int s, int e){
		int i,j;
		for(i=s;i<e;i++)
				for(j=0;j<YMAX;j++){
						jx[i][j]=0; jy[i][j]=0; rho[i][j]=1; f[0][i][j]=0;fprop[0][i][j]=0;
				}
}
void initFValues(int s, int e){
		int x,y;
		float rhxy,uxy,vxy;
		for(x=s;x<e;x++) {
				for(y=0;y<YMAX;y++){
						u[x][y]=jx[x][y]/rho[x][y];
						v[x][y]=jy[x][y]/rho[x][y];
						rhxy=rho[x][y]; uxy=u[x][y]; vxy=v[x][y];
						f[1][x][y]= (rhxy/9.0)*(1+3.0*uxy+ 4.5*uxy*uxy-1.5*(uxy*uxy+vxy*vxy));
						f[2][x][y]= (rhxy/9.0)*(1+3.0*uxy+ 4.5*uxy*uxy-1.5*(uxy*uxy+vxy*vxy));
						f[3][x][y]= (rhxy/9.0)*(1-3.0*uxy+ 4.5*uxy*uxy-1.5*(uxy*uxy+vxy*vxy));
						f[4][x][y]= (rhxy/9.0)*(1-3.0*uxy+ 4.5*uxy*uxy-1.5*(uxy*uxy+vxy*vxy));
						f[5][x][y]=(rhxy/36.0)*(1+3.0*( uxy+vxy)+4.5*( uxy+uxy)*(uxy+uxy)-1.5*(uxy*uxy+vxy*vxy));
						f[6][x][y]=(rhxy/36.0)*(1+3.0*(-uxy+vxy)+4.5*(-uxy+uxy)*(uxy+uxy)-1.5*(uxy*uxy+vxy*vxy));
						f[7][x][y]=(rhxy/36.0)*(1-3.0*( uxy+vxy)+4.5*( uxy+uxy)*(uxy+uxy)-1.5*(uxy*uxy+vxy*vxy));
						f[8][x][y]=(rhxy/36.0)*(1+3.0*( uxy-vxy)+4.5*( uxy-uxy)*(uxy+uxy)-1.5*(uxy*uxy+vxy*vxy));
						f[9][x][y]= (rhxy*4.0/9.0)*(1-1.5*(uxy*uxy+vxy*vxy));
				}
		}
}
void CopyFpropAndFeq(int s, int e){
		int x,y,z;
		for(x=s;x<e;x++){
				for(y=0;y<YMAX;y++){
						for(z=0;z<9;z++)
								fprop[z][x][y]=f[z][x][y];
						feq[z][x][y]=0;
				}
		}
}
void Init(){
		int i,x;
		UX=(0.5*FORCING*YMAX*YMAX)/KVISC;
		for(i=0;i<XMAX;i=i+STEP){
				JxJyRhoFFprop(i,i+STEP);
		}
		for(x=0;x<XMAX;x=x+STEP){
				initFValues(x,x+STEP);
		}
		/* Assign first fprop to equilibrium distributions.*/
		for(x=0;x<XMAX;x=x+STEP){
				CopyFpropAndFeq(x,x+STEP);
		}
}
void ValueUpdateMain(int s, int e){
		int x,y,z;
		float rhxy,uxy,vxy;
		for(x=s;x<e;x++){
				for(y=2;y<YMAX-1;y++){
						rhxy=rho[x][y]; uxy=u[x][y];vxy=v[x][y];
						feq[1][x][y]=(rhxy/9.0) * (1+ 3*uxy +4.5*uxy*uxy -1.5*(uxy*uxy+vxy*vxy) );
						feq[2][x][y]=(rhxy/9.0) * (1+ 3*vxy +4.5*vxy*vxy -1.5*(uxy*uxy+vxy*vxy) );
						feq[3][x][y]=(rhxy/9.0) * (1- 3*uxy +4.5*uxy*uxy -1.5*(uxy*uxy+vxy*vxy) );
						feq[4][x][y]=(rhxy/9.0) * (1- 3*vxy +4.5*vxy*vxy -1.5*(uxy*uxy+vxy*vxy) );
						feq[5][x][y]=(rhxy/36.0) * (1+3*(uxy+vxy)+4.5*(uxy+vxy)*(uxy+vxy)-1.5*(uxy*uxy+vxy*vxy) );
						feq[6][x][y]=(rhxy/36.0) * (1+3*(-uxy+vxy)+4.5*(-uxy+vxy)*(-uxy+vxy)-1.5*(uxy*uxy+vxy*vxy));
						feq[7][x][y]=(rhxy/36.0) * (1-3*(uxy+vxy)+4.5*(uxy+vxy)*(uxy+vxy)-1.5*(uxy*uxy+vxy*vxy) );
						feq[8][x][y]=(rhxy/36.0) * (1-3*(uxy-vxy)+4.5*(uxy-vxy)*(uxy-vxy)-1.5*(uxy*uxy+vxy*vxy) );
						feq[9][x][y]=(rhxy*4.0/36.0) * (1-1.5*(uxy*uxy+vxy*vxy));
						for(z=0;z<9;z++)
								fprop[z][x][y]= (1.0-OMEGA)*f[z][x][y] + OMEGA * feq[z][x][y];
						fprop[1][x][y]+=force;
						fprop[3][x][y]-=force;
						fprop[5][x][y]+=force;
						fprop[6][x][y]-=force;
						fprop[7][x][y]-=force;
						fprop[8][x][y]+=force;
				}
				for(z=0;z<ZMAX;z++) {
						fprop[z][x][0]=f[z][x][0];
						fprop[z][x][YMAX-1]=f[z][x][YMAX-1];
				}
		}
}
void Propagate(int x_min, int x_max, int y_min, int y_max, int xp, int yp, int plane){
		int x,y;
		for(x=x_min;x<x_max;x++){
				for(y=y_min;y<y_max;y++){
						f[plane][x][y]=fprop[plane][x+xp][y+yp];
				}
		}
}
void BounceBackBoundaryCondition(int ymin, int ymax, int lx, int rx, int plane){
		int y;
		for(y=0;y<YMAX;y++)
				f[plane][lx][y]=fprop[plane][rx][y];
}
void BounceBackBoundaryConditionX(int xmin, int xmax, int y, int plane1, int plane2){
		float temp;
		int x;
		for(x=xmin;x<xmax;x++){
				temp=f[plane1][x][y]; f[plane1][x][y]=f[plane2][x][y]; f[plane2][x][y]=temp;
		}
}
void UpdateRHO(){
		int x,y,wt;
		for (x=0;x<XMAX;x++){
				for(y=2;y<YMAX-1;y++){
						rho[x][y]=0.0;
						for(wt=1;wt<=9;wt++)
								rho[x][y]+=f[wt][x][y];
				}
		}
}
void UpdateJXJY(){
		int x,y;
		for(x=0;x<XMAX;x++){
				for(y=2;y<YMAX-1;y++){
						jx[x][y]=f[1][x][y]-f[3][x][y]+f[5][x][y]-f[6][x][y]-f[7][x][y]+f[8][x][y];
						jx[x][y]=f[1][x][y]-f[3][x][y]+f[5][x][y]-f[6][x][y]-f[7][x][y]+f[8][x][y];
				}
		}
}
void CalculateUxyVxy(int s, int e){
		int x,y;
		for(x=s;x<e;x++){
				for(y=0;y<YMAX;y++){
						u[x][y]=jx[x][y]/rho[x][y];
						v[x][y]=jy[x][y]/rho[x][y];
				}
		}
}
void UpdateFinalUprofile(int s, int e){
		int x,y;
		float upf;
		for(y=s;y<e;y++){
				upf=0.0; //Critical Operation
				for(x=0;x<XMAX;x++){
						upf +=jx[x][y]/rho[x][y];
				}	
				upf /=XMAX ;
				uprofile[y] = upf;
		}
}
void LBM(){
		int iter;
		int x,y;
		FILE *fp;
		fp=fopen("lbmout.txt","w");
		force = FORCING/6.0;
		for (iter=0; iter < MAXiter;iter++) {

				for(x=0;x<XMAX;x=x+STEP){
						CalculateUxyVxy(x,x+STEP);
				}
				for (x=0;x<XMAX;x=x+STEP){
						ValueUpdateMain(x,x+STEP);
				}

				/* %Propagating C1   f(p,y,1) = fprop(p-1,y,1);*/
				Propagate(1,XMAX,0,YMAX,-1,0,1);
				/* %Propagating C2 : f(x,q,2) = fprop(x,q-1,2);*/
				Propagate(0,XMAX,1,YMAX,0,-1,2);
				/* %Propagating C3: f(r,y,3) = fprop(r+1,y,3);*/
				Propagate(0,XMAX-1,0,YMAX,1,0,3);
				/*%Propagating C4: f(x,s,4) = fprop(x,s+1,4);*/
				Propagate(0,XMAX,0,YMAX-1,0,1,4);
				/* %Propagating C5  : f(p,q,5) = fprop(p-1,q-1,5);*/
				Propagate(1,XMAX,1,YMAX, -1,-1,5);
				/* %Propagating C6: f(r,p,6) = fprop(r+1,p-1,6);*/
				Propagate(0,XMAX-1,0,YMAX, 1,-1,6);
				/* %Propagating C7: f(r,s,7) = fprop(r+1,s+1,7); */
				Propagate(0,XMAX-1,1,YMAX-1, 1,1,7);
				/* %Propagating C8 : f(p,s,8) = fprop(p-1,s+1,8);*/ 
				Propagate(1,XMAX,0,YMAX-1,-1,1,8);
				/* %Propagating C9 : f(x,y,9) = fprop(x,y,9); */
				Propagate(0,XMAX,0,YMAX,0,0,9);


				BounceBackBoundaryCondition(0,YMAX,0,XMAX-1,1);
				BounceBackBoundaryCondition(0,YMAX,XMAX-1,0,3);	
				BounceBackBoundaryCondition(1,YMAX-1,0,XMAX-1,5);
				BounceBackBoundaryCondition(1,YMAX-1,XMAX-1,0,6);
				BounceBackBoundaryCondition(0,YMAX-1,XMAX-1,0,7);
				BounceBackBoundaryCondition(0,YMAX-1,0,XMAX-1,8);

				BounceBackBoundaryConditionX(0,XMAX,0,2,4);
				BounceBackBoundaryConditionX(0,XMAX,YMAX-1,2,4);
				BounceBackBoundaryConditionX(0,XMAX,0,5,7);
				BounceBackBoundaryConditionX(0,XMAX,YMAX-1,7,5);
				BounceBackBoundaryConditionX(0,XMAX,0,6,8);
				BounceBackBoundaryConditionX(0,XMAX,YMAX-1,8,6);	

				UpdateRHO();
				UpdateJXJY();
				// File operation need to be done in serial
				if (iter>0) fprintf(fp,"%d ", iter);
				for (y=0;y<YMAX; y+=STEP) {
						UpdateFinalUprofile(y, y+STEP);	
						if (iter>0) fprintf(fp," %5.5f",  uprofile[y]);
				}
				fprintf(fp,"\n");
				
		} /*for each iteration*/
		fclose(fp);
}


int main(){

		Init(); /* Initialize LBM Structure*/
		LBM();  /* Simulate LBM*/

		return 0;

}


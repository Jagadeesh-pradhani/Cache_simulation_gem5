
/******************Poiseuille flow using LBM***********************
  Based on procedures as explained in 'Lattice Gas Cellular Automata and Lattice  Boltzmann
  Models'by Wolf Gladrow  
  */

#include <stdio.h>
#include <string.h>
#include <stdlib.h>

float OMEGA=0.2; /*Relaxation factor*/
float TAU=0.5;

/* Sahu: Increase XMAX and YMAX to increase execution time*/

#define XMAX 64 /*Mesh size in x-direction */
#define YMAX 64 /*Mesh size in y-direction */
#define ZMAX 10 /*Not a-direction */
#define STEP 8 
#define MAXiter 10

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

void JxJyRhoFFprop(int s, int e){
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
                        f[7][x][y]=(rhxy/36.0)*(1+3.0*( uxy-vxy)+4.5*(uxy-vxy)*(uxy-vxy)-1.5*(uxy*uxy+vxy*vxy));
                        f[8][x][y]=(rhxy/36.0)*(1+3.0*(-uxy-vxy)+4.5*(-uxy-vxy)*(-uxy-vxy)-1.5*(uxy*uxy+vxy*vxy));
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

void UpdateUV(){
        int x,y;
        for(x=0;x<XMAX;x++){
                for(y=2;y<YMAX-1;y++){
                        u[x][y]=jx[x][y]/rho[x][y];
                        v[x][y]=jy[x][y]/rho[x][y];
                }
        }
}

void UpdateUProfile(){
        int y;
        for(y=0;y<YMAX;y++)
                uprofile[y]=u[2][y];
}

/* Function to read input file */
void read_input_file(char *filename) {
    FILE *file = fopen(filename, "r");
    if (file == NULL) {
        printf("Error: Could not open file %s \n", filename);
        exit(1);
    }
    fscanf(file, "%f", &OMEGA);
    fscanf(file, "%f", &TAU);
    fclose(file);
}

int main(int argc, char *argv[]){
    if (argc < 2) {
        printf("Usage: %s <input_file> \n", argv[0]);
        return 1;
    }
    read_input_file(argv[1]);
    printf("OMEGA: %f, TAU: %f \n", OMEGA, TAU);
    return 0;
}

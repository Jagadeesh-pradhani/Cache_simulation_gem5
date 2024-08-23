#include <string.h>
#include <stdio.h>
#include<stdlib.h>
#include <math.h>

#define N 8

float C[N][N];
float Ct[N][N];
int Quantum[N][N];

void Initialize()
{
    int i;
    int j;
    float pi = atan( 1.0 ) * 4.0;
    for ( i = 0 ; i < N ; i++ )
        for ( j = 0 ; j < N ; j++ )
            Quantum[ i ][ j ] = 1 + ( ( 1 + i + j )  );
    for ( j = 0 ; j < N ; j++ ) {
        C[ 0 ][ j ] = 1.0 / sqrt( (float) N );
        Ct[ j ][ 0 ] = C[ 0 ][ j ];
    }
    for ( i = 1 ; i < N ; i++ ) {
        for ( j = 0 ; j < N ; j++ ) {
            C[ i ][ j ] = sqrt( 2.0 / N ) *
                          cos( pi * ( 2 * j + 1 ) * i / ( 2.0 * N ) );
            Ct[ j ][ i ] = C[ i ][ j ];
        }
    }
}

void DCT( int input[N][N], int output[ N ][ N ])
{
    float temp[ N ][ N ];
    float temp1;
    int i;
    int j;
    int k;

/*  MatrixMultiply( temp, input, Ct ); */
    for ( i = 0 ; i < N ; i++ ) {
        for ( j = 0 ; j < N ; j++ ) {
            temp[ i ][ j ] = 0.0;
            for ( k = 0 ; k < N ; k++ )
                 temp[ i ][ j ] += ( (int) input[ i ][ k ] - 128 ) *
                                   Ct[ k ][ j ];
        }
    }

/*  MatrixMultiply( output, C, temp ); */
    for ( i = 0 ; i < N ; i++ ) {
        for ( j = 0 ; j < N ; j++ ) {
            temp1 = 0.0;
            for ( k = 0 ; k < N ; k++ )
                temp1 += C[ i ][ k ] * temp[ k ][ j ];
            output[ i ][ j ] = (int) temp1 ;
        }
    }
}


void IDCT( int input[ N ][ N ], int output[ N ] [N] )
{
    float temp[ N ][ N ];
    float temp1;
    int i;
    int j;
    int k;

//  MatrixMultiply( temp, input, C ); 
    for ( i = 0 ; i < N ; i++ ) {
        for ( j = 0 ; j < N ; j++ ) {
            temp[ i ][ j ] = 0.0;
            for ( k = 0 ; k < N ; k++ )
                temp[ i ][ j ] += input[ i ][ k ] * C[ k ][ j ];
        }
    }

//  MatrixMultiply( output, Ct, temp ); 
    for ( i = 0 ; i < N ; i++ ) {
        for ( j = 0 ; j < N ; j++ ) {
            temp1 = 0.0;
            for ( k = 0 ; k < N ; k++ )
                temp1 += Ct[ i ][ k ] * temp[ k ][ j ];
            temp1 += 128.0;
            if ( temp1 < 0 )
                output[ i ][ j ] = 0;
            else if ( temp1 > 255 )
                output[ i ][ j ] = 255;
            else
                output[ i ][ j ] = (int) temp1;
        }
    }
}


/////////////////////////////////
typedef struct BMP{
  
  unsigned short bType;           /* Magic number for file */
  unsigned int   bSize;           /* Size of file */
  unsigned short bReserved1;      /* Reserved */
  unsigned short bReserved2;      /* ... */
  unsigned int   bOffBits;        /* Offset to bitmap data */

  unsigned int  bISize;           /* Size of info header */
  unsigned int  bWidth;          /* Width of image */
  unsigned int  bHeight;         /* Height of image */
  unsigned short bPlanes;         /* Number of color planes */
  unsigned short bBitCount;       /* Number of bits per pixel */
  unsigned int  bCompression;    /* Type of compression to use */
  unsigned int  bSizeImage;      /* Size of image data */
  int           bXPelsPerMeter;  /* X pixels per meter */
  int      	bYPelsPerMeter;  /* Y pixels per meter */
  unsigned int  bClrUsed;        /* Number of colors used */
  unsigned int  bClrImportant;   /* Number of important colors */
}BMP;

int Y[8][8],U[8][8],V[8][8],R[8][8],G[8][8],B[8][8];
int YZ[64],UZ[64],VZ[64];
unsigned char *RGB;
unsigned char *RGB1;
BMP bmp;
BMP *b;
int h,w,H,W,current_blk =0;
static int output_writtensofar=0;

/*--------------------------------------------------------*/


int Read_BMP_Header(char *filename,int *h,int *w,BMP *bmp) 
{
   FILE *f;
   int *p;
      f=fopen(filename,"rb");
      /*printf("\nReading BMP Header ");*/
	fread(&bmp->bType,sizeof(unsigned short),1,f);
	printf(" bmp->Type=%d",bmp->bType); 
	
	
      p=(int *)bmp;
	  fread(p+1,sizeof(BMP)-4,1,f);
      if (bmp->bType != 19778) {
/*      printf("Error, not a BMP file!\n");*/
      return 0;
      } 
      *w = bmp->bWidth;
      *h = bmp->bHeight;
	printf(" bmp->bHeight=%d",bmp->bHeight); 
	printf(" bmp->bWidth=%d",bmp->bWidth); 
	printf(" bmp->bOffBits=%d\n",bmp->bOffBits); 
  return 1;
}

void Read_BMP_Data(char *filename,int *h,int *w,BMP *bmp){
		
      int i,j,i1,H,W,Wp,PAD;
      FILE *f;
/*      printf("\nReading BMP Data ");*/
      f=fopen(filename,"rb");
      fseek(f, 0, SEEK_SET);
      fseek(f, bmp->bOffBits, SEEK_SET);
      W = bmp->bWidth;
      H = bmp->bHeight;
/*      printf("\nheight = %d width= %d\n ",H,W);*/
      PAD = (3 * W) % 4 ? 4 - (3 * W) % 4 : 0;
      Wp = 3 * W + PAD;
/*/printf("%d" ,Wp*H);*/
      RGB = (unsigned char *)malloc(Wp* H * sizeof(unsigned char));
      RGB1 = (unsigned char *)malloc(3*Wp* H * sizeof(unsigned char));
      fread(RGB, sizeof(unsigned char), Wp * H, f);
   fclose(f);
  }

void  getBlock(int current_blk) {

	unsigned int H,W,num_blk,count,placey,placex,i,j;
	int i1,Wp,PAD;
	W = bmp.bWidth;
        H = bmp.bHeight;
        PAD = (3 * W) % 4 ? 4 - (3 * W) % 4 : 0;
        Wp = 3 * W + PAD;
 	if(bmp.bHeight%8)  H= bmp.bHeight + (8 - bmp.bHeight%8);
	   else  H= bmp.bHeight; 
	if(bmp.bWidth%8)  W= bmp.bWidth +(8-bmp.bWidth%8);
    	    else   W= bmp.bWidth ;
   	num_blk=H*W/64;
	          count=current_blk;
	  /*/Actual 8X8 Matrix formulation and Sending Data and YUV Conversion.*/
	     if(count <num_blk){
			  placey=(count/(W/8))*8;
			  placex=(count%(W/8))*8;
			  for(i=0;i<8;i++){
			  	  for(j=0;j<8;j++){
				  i1=(placey+i)*Wp + (placex+j)*3;
				  if(((placex+j)<bmp.bWidth)&&((placey+i)<bmp.bHeight)){
				  		  R[i][j]=(unsigned char)RGB[i1];
				  		  G[i][j]=(unsigned char)RGB[i1+1];
				  		  B[i][j]=(unsigned char)RGB[i1+2];
					  }
					  else{
				  		  B[i][j]=0;
				  		  G[i][j]=0;
				  		  R[i][j]=0;
					  }
				  }
			  }
		}
}

void MCT(){
	int i,j;
	int r,g,b,y,u,v;
	for(i=0;i<8;i++){
  		for(j=0;j<8;j++){
	  		r=R[i][j];g=G[i][j];b=B[i][j];
	  		y=(r+2*g+b)/4;
  	  		u=b-g;
          		v=r-g;
	  		R[i][j]=y;G[i][j]=u;B[i][j]=v;	
	  	}
	}
}


struct zigzag {
    int row;
    int col;
} ZigZag[ N * N ] =
{
    {0, 0},
    {0, 1}, {1, 0},
    {2, 0}, {1, 1}, {0, 2},
    {0, 3}, {1, 2}, {2, 1}, {3, 0},
    {4, 0}, {3, 1}, {2, 2}, {1, 3}, {0, 4},
    {0, 5}, {1, 4}, {2, 3}, {3, 2}, {4, 1}, {5, 0},
    {6, 0}, {5, 1}, {4, 2}, {3, 3}, {2, 4}, {1, 5}, {0, 6},
    {0, 7}, {1, 6}, {2, 5}, {3, 4}, {4, 3}, {5, 2}, {6, 1}, {7, 0},
    {7, 1}, {6, 2}, {5, 3}, {4, 4}, {3, 5}, {2, 6}, {1, 7},
    {2, 7}, {3, 6}, {4, 5}, {5, 4}, {6, 3}, {7, 2},
    {7, 3}, {6, 4}, {5, 5}, {4, 6}, {3, 7},
    {4, 7}, {5, 6}, {6, 5}, {7, 4},
    {7, 5}, {6, 6}, {5, 7},
    {6, 7}, {7, 6},
    {7, 7}
};
void zigzag(int YZ[64], int  Y[8][8]){
	int a;
  for ( a = 0; a < 64; a++ )
    	YZ[a]=Y[ZigZag[a].row][ZigZag[a].col];
}


void  RL(int  zdata[64]){
	unsigned char count=0, current,a,reg;
	/*FILE *xx=fopen("enc.dump","a");
	fprintf(xx,"\n RL ........\n");
	*/
		current=zdata[0];
	     	for ( a = 0; a < 64; a++ )
	      	{
		     	if(current==(unsigned char)zdata[a]){
			   	count++;
			   	if(a==63)
			   	{
			    		reg=count;
			    		RGB1[output_writtensofar++]=reg;
		//			fprintf(xx," %d ",reg);
					
			    		reg=current;
					RGB1[output_writtensofar++]=reg;
		//			fprintf(xx," %d ",reg);
			    		
			     	}
			}
		     	else
		     	{
		     		reg=count;
		     		RGB1[output_writtensofar++]=reg;
		//			fprintf(xx," %d ",reg);
		     		reg=current;
		     		RGB1[output_writtensofar++]=reg;
		//			fprintf(xx," %d ",reg);
			   	count=1;
			   	current=(unsigned char)zdata[a];
			   	if(a==63)
			   	{
			    		reg=count;
			    		RGB1[output_writtensofar++]=reg;
		//			fprintf(xx," %d ",reg);
			    		reg=current;
					RGB1[output_writtensofar++]=reg;
		//			fprintf(xx," %d ",reg);
				}
			}
		}
/*/printf("JPEG value of out%d\n",output_writtensofar);*/
//fclose(xx);
}


void quant(int dir){
     int u,v;
    if(dir==0){
    for ( u = 0; u < 8; u++ ){
     for ( v = 0; v < 8; v++ )  {
         Y[u][v]=Y[u][v]*Quantum[u][v];
         U[u][v]=U[u][v]*Quantum[u][v];
         V[u][v]=V[u][v]*Quantum[u][v];       
	}
    }
}else{
   for ( u = 0; u < 8; u++ ){
     for ( v = 0; v < 8; v++ )  {
         Y[u][v]=Y[u][v]/Quantum[u][v];
         U[u][v]=U[u][v]/Quantum[u][v];
         V[u][v]=V[u][v]/Quantum[u][v];       
	}
    }
}
}

void writeJPG(char *filename){
	FILE *fp;
	int *p;
/*/printf("JPEG value of out%d\n",output_writtensofar);*/
	fp=fopen(filename,"wb");  
	p=(int *)b;
	fwrite(&(b->bType),sizeof(unsigned short),1,fp);
	fwrite(p+1,sizeof(BMP)-2,1,fp);
	fwrite(RGB1,sizeof(short),output_writtensofar,fp);
	fclose(fp);
}

int main(){
	int len_y,len_u,len_v,num_blk,h,w;
	b=&bmp;
	Read_BMP_Header("test.bmp",&h,&w,b);
	Read_BMP_Data("test.bmp",&h,&w,b);
	output_writtensofar=0;

	if(bmp.bHeight%8)  H= bmp.bHeight + (8 - bmp.bHeight%8);
	   else  H= bmp.bHeight; 
	if(bmp.bWidth%8)  W= bmp.bWidth +(8-bmp.bWidth%8);
    	    else   W= bmp.bWidth ;

   	num_blk=H*W/64;
	Initialize();
	for(current_blk=0;current_blk<num_blk;current_blk++){
		getBlock(current_blk);
		MCT();
		DCT(R,Y);DCT(G,U);DCT(B,V);
		quant(1);// 1  for forward
		zigzag(YZ,Y);zigzag(UZ,U);zigzag(VZ,V);
		RL(YZ);RL(UZ);RL(VZ);	
	}
	writeJPG("myfile.jpg");
	return 0;
}

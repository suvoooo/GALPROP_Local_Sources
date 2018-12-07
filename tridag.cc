
//**.****|****.****|****.****|****.****|****.****|****.****|****.****|****.****|
// * tridag.cc *                                   galprop package * 4/14/2000 
//**"****!****"****!****"****!****"****!****"****!****"****!****"****!****"****|

/*****************************************************************************\
*                 *** I.Moskalenko (MPE,Garching) *** version of 04.08.99 *** *
* PURPOSE:                                                                    *
*   tridag(a,b,c,r,u,n) - solves for vector u[] the tridiagonal linear set:   *
*                 | b1 c1  0 ...                     | |  u1  | |  r1  |      *
*                 | a2 b2 c2 ...                     | |  u2  | |  r2  |      *
*                 |          ...                     |*|  ... |=|  ... |      *
*                 |          ... a(n-1) b(n-1) c(n-1)| |u(n-1)| |r(n-1)|      *
*                 |          ...   0    a(n)   b(n)  | | u(n) | | r(n) |      *
* INPUT parameters:                                                           *
*   (float ) a[],b[],c[] are the input vectors (j=1..n),                      *
*   (int)    n - the length of the vectors.                                   *
* OUTPUT:                                                                     *
*   returns 0 when solution exists, in this case                              *
*      (float ) u[] is a vector of solutions (j=1..n);                        *
*   returns 1 for a bad defined matrix;                                       * 
*   returns 2 when solution doesn't exist.                                    *
* NOTES:                                                                      *
* - The algorithm  is adopted from the "Numerical Recipes" by W.H.Press et al.*
\*****************************************************************************/
// 19991124  replaced dynamic gam with static                                 
// 19991203  float version (double version copied to tridag_double.cc)        
using namespace std;//AWS20050624
#include <cstdio>

//**.****|****.****|****.****|****.****|****.****|****.****|****.****|****.****|

int tridag(float  a[], float  b[], float  c[], float  r[], float  u[], int n) 
{
#define NMAX 1000
   static float gam[NMAX];
   int j;
   float bet;

   if(n>NMAX) { printf("\ntridag ... n>NMAX\n\n"); return 3; }
   if(b[0] == 0) { printf("\ntridag ... rewrite equations\n\n"); return (1); }

   bet = b[0];
   u[0]= r[0]/bet;
   for(j=1; j<n; j++)
   {
      gam[j] = c[j-1]/bet;
      bet    = b[j]-a[j]*gam[j];
      if(bet == 0.) { printf("\ntridag ... tridag failed\n\n"); return (2); }
      u[j] = (r[j]-a[j]*u[j-1])/bet;
   }
   for(j=n-2; j>=0; u[j]=u[j]-gam[j+1]*u[j+1],j--);
   return (0);
}

//**.****|****.****|****.****|****.****|****.****|****.****|****.****|****.****|

int tridag(double a[], double b[], double c[], double r[], double u[], int n) 
{
   int j,key;
   double bet;
	 //Gulli20070810
//#define NMAX 1000
   //static double gam[NMAX];
	 double *gam;
	 gam = new double[n];
   //if(n>NMAX) {printf("\ntridag ... n>NMAX\n\n");return 3;}
   if(b[0] == 0) {  printf("\ntridag ... rewrite equations\n\n"); return (1); }

   bet = b[0];
   u[0]= r[0]/bet;
   for(j=1; j<n; j++) 
   {
      gam[j] = c[j-1]/bet;
      bet    = b[j]-a[j]*gam[j];
      if(bet == 0.) { printf("\ntridag ... tridag failed\n\n"); return (2); }
      u[j] = (r[j]-a[j]*u[j-1])/bet;
   }
   for(j=n-2; j>=0; u[j]=u[j]-gam[j+1]*u[j+1],j--);
	 delete[] gam;
   return (0);
}

//**.****|****.****|****.****|****.****|****.****|****.****|****.****|****.****|

/************** TESTS **********************************
int tridag( float*, float*, float*, float*, float*,int);
int tridag(double*,double*,double*,double*,double*,int);
main() 
{
   int i,n = 5, n_sym;
   float a[n],b[n],c[n],r[n],u[n],a_sym[n],b_sym[n],c_sym[n],r_sym[n];

   for(i=0, printf("\n"); i<n; i++)
   {
      a[i] =-1.;
      b[i] = 1.;
      c[i] = 1.;
      r[i] = 1.;
//printf("%4d%12.2e%12.2e%12.2e  %12.2e\n",i,a[i],b[i],c[i],r[0]);
   }
   a[n/2] = c[n/2];
   for(i=0, printf("\n"); i<n/2; i++)
   {
      a[i] = c[n-i-1];
      b[i] = b[n-i-1];
      c[i] = a[n-i-1];
      r[i] = r[n-i-1];
//printf("%4d%12.2e%12.2e%12.2e  %12.2e\n",i,a[i],b[i],c[i],r[0]);
   }
/////////////////////////////
      n=3; n_sym=2*n-1;
      for(i=0; i<n; i++)
      {
         a[i] = i+1;
         b[i] = -2*a[i];
         c[i] = a[i]+2;
         r[i] = i+4;
      }
      a[0]=c[0]; // essential for tridag_sym and creating a_sym
      for(i=0; i<n; i++)
      {
         a_sym[i+n-1] = a[i];
         b_sym[i+n-1] = b[i];
         c_sym[i+n-1] = c[i];
         r_sym[i+n-1] = r[i];
      }
      for(i=0; i<n-1; i++)
      {
         a_sym[i  ] = c_sym[n_sym-i-1];
         b_sym[i  ] = b_sym[n_sym-i-1];
         c_sym[i  ] = a_sym[n_sym-i-1];
         r_sym[i  ] = r_sym[n_sym-i-1];
      }
      for(i=0; i<n_sym; i++)
      {
         a[i] = a_sym[i];
         b[i] = b_sym[i];
         c[i] = c_sym[i];
         r[i] = r_sym[i];
      }
      n=n_sym;
/////////////////////////////
   printf("\n%12.2e%12.2e%12.2e  %12.2e\n", 0., b[0],c[0],r[0]);
   printf(  "%12.2e%12.2e%12.2e  %12.2e\n",a[1],b[1],c[1],r[1]);
   printf(  "%12.2e%12.2e%12.2e  %12.2e\n",a[2],b[2],c[2],r[2]);
   printf(  "%12.2e%12.2e%12.2e  %12.2e\n",a[3],b[3],c[3],r[3]);
   printf(  "%12.2e%12.2e%12.2e  %12.2e\n",a[4],b[4], 0. ,r[4]);

   tridag(a,b,c,r,u,n);

   for(i=0,printf("\n"); i<n; printf(" %e ",u[i]),i++); printf("\n");
   return (0);
}
****************************************************************/


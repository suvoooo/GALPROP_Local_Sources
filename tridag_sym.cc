
//**.****|****.****|****.****|****.****|****.****|****.****|****.****|****.****|
// * tridag_sym.cc *                               galprop package * 4/05/2002 
//**"****!****"****!****"****!****"****!****"****!****"****!****"****!****"****|
/////////////////////////////////////// revised version, symmetric about low end
using namespace std;//AWS20050624
#include<cstdio>

//**.****|****.****|****.****|****.****|****.****|****.****|****.****|****.****|

int tridag_sym(float  a[], float  b[], float  c[], float  r[], float  u[], int n) 
{
#define NMAX 1000
   static float gam[NMAX], u_sym[NMAX];
   int j,k;
   float bet;
   int n_sym=2*n-1; // dimension of expanded matrix
   int m=n-1;       // index of central bin of expanded matrix

   if(n_sym>NMAX) { printf("\ntridag_sym ... n>NMAX\n\n"); return 3; }
   if(b[0] == 0)  { printf("\ntridag_sym ... rewrite equations\n\n"); return (1); }

   bet     = b[m];
   u_sym[0]= r[m]/bet;
   for(j=1; j<=m; j++)
   {
      k=n-j-1;
      gam[j] = a[k+1]/bet;
      bet    = b[k]-c[k]*gam[j];
      if(bet == 0.) { printf("\ntridag_sym ... tridag failed\n\n"); return (2); }
      u_sym[j] = (r[k]-c[k]*u_sym[j-1])/bet;
   }
   for(j=m+1; j<n_sym; j++) 
   {
      k=j-m;
      gam[j] = c[k-1]/bet;
      bet    = b[k]-a[k]*gam[j];      
      if(bet == 0.) { printf("\ntridag_sym ... tridag_sym failed\n\n"); return (2); }
      u_sym[j] = (r[k]-a[k]*u_sym[j-1])/bet;
   }
   for(j=n_sym-2; j>=m; u_sym[j]=u_sym[j]-gam[j+1]*u_sym[j+1],j--);
   for(j=0; j<n; j++) u[j] = u_sym[j+m];
   return (0);
}

//**.****|****.****|****.****|****.****|****.****|****.****|****.****|****.****|

int tridag_sym(double a[], double b[], double c[], double r[], double u[], int n)  //IMOS20030217
{
#define NMAX 1000
   static double gam[NMAX], u_sym[NMAX];
   int j,k;
   double bet;
   int n_sym=2*n-1; // dimension of expanded matrix
   int m=n-1;       // index of central bin of expanded matrix

   if(n_sym>NMAX) { printf("\ntridag_sym ... n>NMAX\n\n"); return 3; }
   if(b[0] == 0)  { printf("\ntridag_sym ... rewrite equations\n\n"); return (1); }

   bet     = b[m];
   u_sym[0]= r[m]/bet;
   for(j=1; j<=m; j++)
   {
      k=n-j-1;
      gam[j] = a[k+1]/bet;
      bet    = b[k]-c[k]*gam[j];
      if(bet == 0.) { printf("\ntridag_sym ... tridag failed\n\n"); return (2); }
      u_sym[j] = (r[k]-c[k]*u_sym[j-1])/bet;
   }
   for(j=m+1; j<n_sym; j++) 
   {
      k=j-m;
      gam[j] = c[k-1]/bet;
      bet    = b[k]-a[k]*gam[j];      
      if(bet == 0.) { printf("\ntridag_sym ... tridag_sym failed\n\n"); return (2); }
      u_sym[j] = (r[k]-a[k]*u_sym[j-1])/bet;
   }
   for(j=n_sym-2; j>=m; u_sym[j]=u_sym[j]-gam[j+1]*u_sym[j+1],j--);
   for(j=0; j<n; j++) u[j] = u_sym[j+m];
   return (0);
}

/************** TESTS **********************************

main() 
{
   int i,n = 7, n1=n/2+1,n_sym;
   float a[n],b[n],c[n],r[n],u[n],a_sym[n],b_sym[n],c_sym[n],r_sym[n];

   for(i=0; i<n; i++)
   {
      a[i] =-1.;
      b[i] = 1.;
      c[i] = 1.;
      r[i] = 1.;
   }
   a[0]=c[0];
/////////////////////////////
      n=n1; n_sym=2*n-1;
      for(i=0; i<n; i++)
      {
         a[i] = i+1;
         b[i] = -2*a[i];
         c[i] = a[i]+2;
         r[i] = i+4;
      }
      a[0]=c[0]; // essential for tridag_sym and creating a_sym
/////////////////////////////
                      printf("\n%12.2e%12.2e%12.2e  %12.2e\n", 0., b[0],c[0],r[0]);
   for(i=1;i<n-1;i++) printf(  "%12.2e%12.2e%12.2e  %12.2e\n",a[i],b[i],c[i],r[i]);
                      printf(  "%12.2e%12.2e%12.2e  %12.2e\n",a[n-1],b[n-1], 0. ,r[n-1]);
   tridag_sym(a,b,c,r,u,n1);

   for(i=n1-1,printf("\n"); i>0; printf(" %e ",u[i]),i--);
   for(i=0; i<n1; printf(" %e ",u[i]),i++); printf("\n");
   return (0);
}
*/




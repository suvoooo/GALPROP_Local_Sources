
//**.****|****.****|****.****|****.****|****.****|****.****|****.****|****.****|
// * tridag_ext.cc *                               galprop package * 3/29/2002
//**"****!****"****!****"****!****"****!****"****!****"****!****"****!****"****|
using namespace std;//AWS20050624
#include <cstdio>
#include <iostream>

int tridag_ext_init=0;
//Gulli20070810 Changed from float to double
int tridag_ext(double  a[], double  b[], double  c[], double  r[], double  u[], int n,int nk) 
{
   int j,key,k,kk;
	 //Gulli20070810 Why on earth are they static???
   //static int NMAX;
   //static double  *gam;
   //static double  *bet;
	 double *gam;
	 double *bet;

	 gam = new double[n*nk];
	 bet = new double[nk];

	 /* //Gulli20070810
   if(tridag_ext_init==0)
   {
      tridag_ext_init=1;
      NMAX=n*nk;
      cout<<"  tridag_ext: initializing arrays: n="<<n<<" nk= "<<nk<<endl;
      cout<<"  tridag_ext: initializing arrays to dimension ="<<NMAX<<endl;
      gam  =new double[NMAX];
      bet  =new double[NMAX];
   }
//   cout<<" tridag_ext: : n="<<n<<" nk= "<<nk<<" NMAX="<<NMAX<<endl;

// the dimensions in x,y,z,p are different so may need to reallocate: AWS20010123
   if(n*nk>NMAX)
   { 
      cout<<" tridag_ext ... n*nk>NMAX"<<endl;
      cout<<" tridag_ext: : n="<<n<<" nk= "<<nk<<" NMAX="<<NMAX<<endl;
      NMAX=n*nk;
      cout<<"  tridag_sym_ext: re-initializing arrays: n="<<n<<" nk= "<<nk<<endl;
      cout<<"  tridag_sym_ext: re-initializing arrays to dimension ="<<NMAX<<endl;
      delete[] gam;
      delete[] bet;
      gam  =new double[NMAX];
      bet  =new double[NMAX];
   }
	 */
   if(b[0] == 0) { printf("\ntridag ... rewrite equations\n\n"); return (1); }

#pragma vdir nodep   
   for(k=0,kk=0;k<nk;k++,kk+=n) bet[k]=b[kk];
#pragma vdir nodep
   for(k=0,kk=0;k<nk;k++,kk+=n) u[kk]= r[kk]/bet[k];
   for(j=1; j<n; j++) 
   {
#pragma vdir nodep
      for(k=0,kk=0; k<nk; k++,kk+=n) gam[j+kk] = c[j-1+kk]/bet[k];

#pragma vdir nodep
      for(k=0,kk=0;k<nk;k++,kk+=n) bet[k]    = b[j+kk]-a[j+kk]*gam[j+kk];
//      if(bet == 0.) { printf("\ntridag ... tridag failed\n\n"); return (2); }

#pragma vdir nodep
     for(k=0,kk=0;k<nk;k++,kk+=n)      
        u[j+kk] = (r[j+kk]-a[j+kk]*u[j-1+kk])/bet[k];
   }

   for(j=n-2; j>=0;j--)
#pragma vdir nodep
      for(k=0,kk=0;k<nk;k++,kk+=n)
         u[j+kk]=u[j+kk]-gam[j+1+kk]*u[j+1+kk];

	 delete[] gam;
	 delete[] bet;
   return (0);
}

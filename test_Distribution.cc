
//**.****|****.****|****.****|****.****|****.****|****.****|****.****|****.****|
// * test_Distribution.cc *                        galprop package * 4/14/2000 
//**"****!****"****!****"****!****"****!****"****!****"****!****"****!****"****|
using namespace std;//AWS20050624

#include"galprop_classes.h"
#include"galprop_internal.h"



int test_Distribution(){

cout<<">>>>test_Distribution"<<endl;
int n_xgrid=50;
int n_ygrid=100;
int n_zgrid=100;
int n_pgrid=10;

 cout<<"enter integer"<<endl;
int q;
 cin>>q;

Distribution dist(n_xgrid,n_ygrid,n_zgrid,n_pgrid);
int size=n_xgrid*n_ygrid*n_zgrid*n_pgrid*8     ;
 cout<<"size in MB="<<size/1.0e6<<endl;
 cout<<"enter integer"<<endl;

 cin>>q;
 cout<<"creating double a[size/8]"<<endl;
double *a;
a=new double[size/8];
for(int i=0;i<size/8;i++)a[i]=i;
 cout<<"enter integer"<<endl;

 cin>>q;
 cout<<"creating double b[size/8]"<<endl;
double *b;
b=new double[size/8];
for(int i=0;i<size/8;i++)b[i]=i;
 cout<<"enter integer"<<endl;
 cin>>q;
cout<<"<<<<test_Distribution"<<endl;

return 0;
}

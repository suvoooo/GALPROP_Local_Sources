
//**.****|****.****|****.****|****.****|****.****|****.****|****.****|****.****|
// * test_float_accuracy.cc *                      galprop package * 4/14/2000 
//**"****!****"****!****"****!****"****!****"****!****"****!****"****!****"****|
using namespace std;//AWS20050624
#include<iostream>
int test_float_accuracy(){
  cout<<">>>>test_float_accuracy"<<endl;
  cout<<">>>>float:"<<endl;
float x=1.0;
float r=0.1;
 for (int i=0;i<100;i++){x*=r;cout<<x<<" ";}
			  cout<<endl;
  cout<<">>>>double:"<<endl;
double xx=1.0;
double rr=0.1;
 for (int i=0;i<700;i++){xx*=rr;cout<<xx<<" ";}
			  cout<<endl;
  cout<<"<<<<test_float_accuracy"<<endl;
return 0;
}

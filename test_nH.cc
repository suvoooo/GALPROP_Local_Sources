
//**.****|****.****|****.****|****.****|****.****|****.****|****.****|****.****|
// * test_nH.cc *                                  galprop package * 4/14/2000 
//**"****!****"****!****"****!****"****!****"****!****"****!****"****!****"****|
using namespace std;//AWS20050624
#include"galprop_classes.h"
#include"galprop_internal.h"
int test_nH(){

  cout<<">>>>test_nH\n";
double R,z;
double dR=3,dz=.1;
double Rmax=30,zmin=-5,zmax=+5;

 cout<<"\n======================== nHI\n";
 cout<<"z/         R ";
 for(double R=0;   R<=Rmax;R+=dR){cout<<R<<"       ";};cout<<endl;
 for(double z=zmin;z<=zmax;z+=dz){
   cout<<z<<"       ";
   for(double R=0;   R<=Rmax;R+=dR){
   cout<<nHI(R,z)<<" ";
}
   cout<<endl;
}

 dR=1;dz=.1;
 Rmax=10;zmin=-3;zmax=+3;
 cout<<"\n======================== nH2\n";
 cout<<"z/         R ";
 for(double R=0;   R<=Rmax;R+=dR){cout<<R<<"       ";};cout<<endl;
 for(double z=zmin;z<=zmax;z+=dz){
   cout<<z<<"       ";
   for(double R=0;   R<=Rmax;R+=dR){
   cout<<nH2(R,z)<<" ";
}
   cout<<endl;
}

 dR=2;dz=.1;
 Rmax=20;zmin=-10;zmax=+10;
 cout<<"\n======================== nHII\n";
 cout<<"z/         R ";
 for(double R=0;   R<=Rmax;R+=dR){cout<<R<<"       ";};cout<<endl;
 for(double z=zmin;z<=zmax;z+=dz){
   cout<<z<<"       ";
   for(double R=0;   R<=Rmax;R+=dR){
   cout<<nHII(R,z)<<" ";
}
   cout<<endl;
}


 ////////////////// averages ////////////////////
 double xmin=-21;
 double xmax=+21;
 double dx=3;
 double y=0;
 double dzz=.01;
 zmin=-1;zmax=+1;dz=.1;
 
 cout<<"\n======================== nHI_av\n";
 cout<<"z/         x ";
 for(double x=xmin;x<=xmax;x+=dx){cout<<x<<"       ";};cout<<endl;
 for(double z=zmin;z<=zmax;z+=dz){
   cout<<z<<"       ";
   for(double x=xmin;   x<=xmax;x+=dx){
   cout<<nH_av(x,y,z,dz,dzz,nHI)<<" "; //IMOS20080114
   //cout<<nHI_av(x,y,z,dz,dzz)<<" ";
}
   cout<<endl;
}
  
xmin=-10;xmax=10;dx=1;

 cout<<"\n======================== nH2_av\n";
 cout<<"z/         x ";
 for(double x=xmin;x<=xmax;x+=dx){cout<<x<<"       ";};cout<<endl;
 for(double z=zmin;z<=zmax;z+=dz){
   cout<<z<<"       ";
   for(double x=xmin;   x<=xmax;x+=dx){
      cout<<nH_av(x,y,z,dz,dzz,nH2)<<" "; //IMOS20080114
     //cout<<nH2_av(x,y,z,dz,dzz)<<" ";
}
   cout<<endl;
}

zmin=-10;zmax=10;
 cout<<"\n======================== nHII_av\n";
 cout<<"z/         x ";
 for(double x=xmin;x<=xmax;x+=dx){cout<<x<<"       ";};cout<<endl;
 for(double z=zmin;z<=zmax;z+=dz){
   cout<<z<<"       ";
   for(double x=xmin;   x<=xmax;x+=dx){
        cout<<nH_av(x,y,z,dz,dzz,nHII)<<" "; //IMOS20080114
     //   cout<<nHII_av(x,y,z,dz,dzz)<<" ";
}
   cout<<endl;
}



  cout<<"<<<<test_nH\n";
return 0;}

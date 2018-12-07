
//**.****|****.****|****.****|****.****|****.****|****.****|****.****|****.****|
// * sigma_boron_dec_heinbach_simon.cc *           galprop package * 4/14/2000 
//**"****!****"****!****"****!****"****!****"****!****"****!****"****!****"****|
using namespace std;//AWS20050624
#include<iostream>

double  sigma_boron_dec_heinbach_simon(int IZ,int IA,int JZ,int JA,double EJ){

/*
  Heinbach and Simon 1995 ApJ 441, 209 Table 1
  
  This gives 12C,16O -> 10C+10B, 11C+11B from a fit to experimental data
  In this subroutine the cross-sections are not split into 10C,10B,11C,11B,
  but the full decayed cross section is returned
  (cf routine sigma_boron_heinbach_simon where this split is made)

  EJ= kinetic energy per nucleon of 12C of 16O
*/


#define  npoints 17 

                            
  //      sig_12_10             12C->10C+10B
  //      sig_12_11             12C->11C+11B
  //      sig_16_10             16O->10C+10B
  //      sig_16_11             16O->11C+11B

double QJ;
int   iuse,ie;

//cout<<  ">>>>    sigma_boron_dec_heinbach_simon"<<endl;



iuse=0;
 
// NB here was an error in f90 version

if(IA==12 &&    IZ == 6 &&    JA == 10 &&    JZ == 5)iuse=1;//12C->10B
if(IA==12 &&    IZ == 6 &&    JA == 11 &&    JZ == 5)iuse=1;//12C->11B

 

if(IA==16 &&    IZ == 8 &&    JA == 10 &&    JZ == 5)iuse=1;//16O->10B
if(IA==16 &&    IZ == 8 &&    JA == 11 &&    JZ == 5)iuse=1;//16O->11B


if (iuse == 0 ){    
  QJ=-99999.;
  return QJ;
}    

float e[]       ={  100.,200.,300.,400.,500.,600.,700.,800.,900.,1000.,1200.,1400.,1600.,1800.,2000.,2400.,2800.};

float sig_12_10[]={  21.1,16.8,16.1,16.6,22.5,24.6,24.3,23.7,23.2,22.8, 22.0, 21.5, 21.3, 21.1, 20.8, 20.8, 20.8  };
float sig_12_11[]={  67.4,59.2,62.7,65.9,67.5,64.3,61.2,57.9,53.9,51.3, 52.7, 55.4, 57.1, 56.8, 56.6, 56.4, 56.1  };
float sig_16_10[]={  13.1, 8.2, 9.1, 9.8, 9.8, 9.9, 9.9, 9.9, 9.9, 9.9, 10.0, 10.0, 10.0, 10.0, 10.0, 10.0, 10.0  };
float sig_16_11[]={  32.0,17.0,22.3,26.0,26.0,26.0,26.0,26.0,26.0,26.0, 26.0, 26.0, 26.0, 26.0, 26.0, 26.0, 26.0  };

ie=0;
for(ie=0;ie<npoints-1;ie++){
if(EJ > e[ie] &&  EJ <= e[ie+1]  ) break;
}     
//cout<< EJ<<" "<<ie<<" "<<e[ie]<<endl;

float *sig_buff;

if (IA ==12 &&    JA == 10)sig_buff=sig_12_10;
if (IA ==12 &&    JA == 11)sig_buff=sig_12_11;
if (IA ==16 &&    JA == 10)sig_buff=sig_16_10;
if (IA ==16 &&    JA == 11)sig_buff=sig_16_11;



QJ=sig_buff[ie]+(EJ-e[ie])* ( sig_buff[ie+1]-sig_buff[ie] ) /(e[ie+1]-e[ie]);

if(EJ <= e[0]      )   QJ= sig_buff[0];
if(EJ >= e[npoints-1]) QJ= sig_buff[npoints-1];

 

//cout<<  IZ<<" "<<IA<<" "<<JZ<<" "<<JA<<" "<<EJ<<" "<<QJ<<endl;

//cout<<  "<<<<<< sigma_boron_dec_heinbach_simon"<<endl;

return QJ;
}   

int test_sigma_boron_dec_heinbach_simon(){
cout<<  ">>>> test_boron_dec_heinbach_simon"<<endl;
int IZ; 
int IA;
int JZ;
int JA;
double EJ;
IA=12;IZ=6; 
JA=10;JZ=5;
 for(EJ=50.;EJ<1.e5;EJ*=1.5){
   cout<<"IZ IA JZ JA EJ sigma "<< IZ<<" "<<IA<<" "<<JZ<<" "<<JA<<" "<<EJ<<" "<<
   sigma_boron_dec_heinbach_simon(IZ,IA,JZ,JA ,EJ)<<endl;
 }
 cout<<endl;
IA=12;IZ=6; 
JA=11;JZ=5;
 for(EJ=50.;EJ<1.e5;EJ*=1.5){
 cout<<"IZ IA JZ JA EJ sigma "<< IZ<<" "<<IA<<" "<<JZ<<" "<<JA<<" "<<EJ<<" "<<
   sigma_boron_dec_heinbach_simon(IZ,IA,JZ,JA ,EJ)<<endl;
 }
 cout<<endl;
IA=16;IZ=8; 
JA=10;JZ=5;
 for(EJ=50.;EJ<1.e5;EJ*=1.5){
 cout<<"IZ IA JZ JA EJ sigma "<< IZ<<" "<<IA<<" "<<JZ<<" "<<JA<<" "<<EJ<<" "<<
   sigma_boron_dec_heinbach_simon(IZ,IA,JZ,JA ,EJ)<<endl;
 }
 cout<<endl;
IA=16;IZ=8; 
JA=11;JZ=5;
 for(EJ=50.;EJ<1.e5;EJ*=1.5){
 cout<<"IZ IA JZ JA EJ sigma "<< IZ<<" "<<IA<<" "<<JZ<<" "<<JA<<" "<<EJ<<" "<<
   sigma_boron_dec_heinbach_simon(IZ,IA,JZ,JA ,EJ)<<endl;
 }
 cout<<endl;
 JA=9; // test case where cross section undefined
 cout<<" undefined case:IZ IA JZ JA EJ sigma "<< IZ<<" "<<IA<<" "<<JZ<<" "<<JA<<" "<<EJ<<" "<<
   sigma_boron_dec_heinbach_simon(IZ,IA,JZ,JA ,EJ)<<endl;
cout<<  "<<<< test_boron_dec_heinbach_simon"<<endl;
return 0;
}

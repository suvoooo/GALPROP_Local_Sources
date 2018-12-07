#include <iostream>   //AWS20050919
#include <cmath>       //AWS20080620
#include <cstring>
#include <cstdlib>
#include <cstdio>

using namespace std;

#include <ErrorLogger.h>

#include"GCR_data.h"

GCR_data::GCR_data() : 
  n(0), 
  entry(0),
  status(0),
  reference(0),
  experiment(0),
  Y_name(0),
  E_low_input(0),
  E_mean_input(0),
  value_input(0),
  err_minus_input(0),
  err_plus_input(0),
  err_type(0),
  Z_numerator(0),
  A_numerator(0),
  Z_denominator(0),
  A_denominator(0),
  comment(0),
  X_units(0),
  Y_units(0),
  E_low(0),
  E_high(0),
  E_mean(0),
  value(0),
  err_minus(0),
  err_plus(0),
  color(0),
  style(0),
  size(0) {

}

GCR_data::~GCR_data() {

  if (n) {
    //Delete entries
    for (int i = 0; i <= N; ++i) {
      delete[] entry[i];
    }
    
    delete[] entry;
    delete[] status;
    delete[] E_low_input;
    delete[] E_high_input;
    delete[] E_mean_input;
    delete[] value_input;
    delete[] err_minus_input;
    delete[] err_plus_input;
    delete[] E_low;
    delete[] E_high;
    delete[] E_mean;
    delete[] value;
    delete[] err_minus;
    delete[] err_plus;
    delete[] color;
    delete[] style;
    delete[] size;
    delete[] err_type;
    
    for (int i=0;i<N;++i) {
	
      delete[] reference[i];
      delete[] experiment[i];
      delete[] X_units[i];
      delete[] Y_units[i];
      delete[] Y_name[i];
      delete[] Z_numerator[i];
      delete[] A_numerator[i];
      delete[] Z_denominator[i];
      delete[] A_denominator[i];
    
    }
    
    for(int i=0;i<n;++i)
      delete comment[i];//free(comment[i]); //WTF????
      
    delete[] reference;
    delete[] experiment;
    delete[] X_units;
    delete[] Y_units;
    delete[] Y_name;
    delete[] Z_numerator;
    delete[] A_numerator;
    delete[] Z_denominator;
    delete[] A_denominator;
    delete[] comment;
   
  }

}

int 
GCR_data::
read(const char* database_file_, 
     char *area_units_, 
     char *energy_units_) {
  
  INFO("Entry");
  //cout<<">>GCR_data::read"<<endl;
  
  int i,j,start;
  char* flag;
  double area_conversion_factor,energy_conversion_factor;
  
  strcpy(database_file,database_file_);
  strcpy(area_units,     area_units_);cout<<area_units<<endl;
  strcpy(energy_units, energy_units_);

  ostringstream buf;
  buf << "Reading from " << database_file;
  INFO(buf.str());
  
  cout<<"area_units   ="<<area_units   <<endl;
  FILE* ft = fopen(database_file,"r"); 
  
  if (0 == ft) {
    
    cout<<"no database file called "<<database_file<<endl; return -1;
      
  }

  n=10000;
  entry=new char*[n];
  entry[0]=new char[1000];
  start=0;
  while(start==0) {

    flag=fgets(entry[0],1000,ft); // read string until newline (Schildt p.222)
    flag=strstr(entry[0],"====");
    
    if(flag!=NULL){cout<<"start of data found"<<endl;start=1;}
    
  }
  delete[] entry[0];
  
  for(i=0;i<n;i++) {

    //cout << i << endl;

    entry[i]=new char[1000];
    memset(entry[i],0,1000);
    
    flag=fgets(entry[i],1000,ft); // read string until newline (Schildt p.222)
    if(flag==NULL){cout<<"end of data after "<<i<<" entries"<<endl;break;}
    //cout<<entry[i]<<endl;
  }
  
  N=i;
  n=i;
  
  fclose(ft);
  
  status           =new char  [n];
  reference        =new char *[n];
  experiment       =new char *[n];
  E_low_input      =new double[n];
  E_high_input     =new double[n];
  E_mean_input     =new double[n];
  value_input      =new double[n];
  err_minus_input  =new double[n];
  err_plus_input   =new double[n];
  err_type         =new char  [n];
  
  Z_numerator  =new int*[n];
  A_numerator  =new int*[n];
  Z_denominator=new int*[n];
  A_denominator=new int*[n];
  
  comment=new char*[n];
  X_units=new char*[n];
  Y_units=new char*[n];
  Y_name =new char*[n];
  
  E_low       =new double[n];
  E_high      =new double[n];
  E_mean      =new double[n];
  value       =new double[n];
  err_minus   =new double[n];
  err_plus    =new double[n];
  
  n_ZA_numerator   =3;
  n_ZA_denominator =3;
  
  color       =new    int[n]; //AWS20060621
  style       =new    int[n]; //AWS20060621
  size        =new double[n]; //AWS20060621
  
  for(i=0;i<n;i++) {
      
    reference    [i] =new char[20];
    experiment   [i]= new char[20];
    X_units      [i]= new char[20];
    Y_units      [i]= new char[20];
    Y_name       [i]= new char[20];
    Z_numerator  [i]= new int [n_ZA_numerator];
    A_numerator  [i]= new int [n_ZA_numerator];
    Z_denominator[i]= new int [n_ZA_denominator];
    A_denominator[i]= new int [n_ZA_denominator];
    //comment      [i] =new char[40 ];
    for (j = 0; j < n_ZA_numerator; ++j){
      Z_numerator  [i][j] = 0;
      A_numerator  [i][j] = 0;
      Z_denominator[i][j] = 0;
      A_denominator[i][j] = 0;
    }
  }
  
  i=0;
  for(j=0;j<n;j++) {

    //cout<<i<<endl;
    //We must handle comment lines separately to avoid buffer overflow.
    sscanf(entry[j],"%c",&status[i]);
    if (status[i] != '#') {
      int comment_start;
      int nent = sscanf(entry[j], "%s%s%s%s%s%le%le%le%le%le%le %c%d.%d%d.%d%d.%d |%d.%d%d.%d%d.%d%n",
			reference[i],experiment[i],X_units[i],Y_units[i],Y_name[i],
			//       &E_low_input[i],&E_high_input[i],&E_mean_input[i],   //AWS20060320
			&E_mean_input[i],  &E_low_input[i],&E_high_input[i], //AWS20060320 to correspond to data format used
			&value_input[i],&err_minus_input[i],&err_plus_input[i],&err_type[i],
			& Z_numerator  [i][0],&A_numerator  [i][0],
			& Z_numerator  [i][1],&A_numerator  [i][1],
			& Z_numerator  [i][2],&A_numerator  [i][2],
			& Z_denominator[i][0],&A_denominator[i][0],
			& Z_denominator[i][1],&A_denominator[i][1],
			& Z_denominator[i][2],&A_denominator[i][2],
			&comment_start
			);
      if (nent < 24 ) {
	cout<<"Format of line broken "<<entry[j]<<endl;
      } else {
	comment[i] = strdup(&entry[j][comment_start]);
	comment[i][strlen(comment[i])-1] = char(0);
	//if (i < j) strcpy(entry[i],entry[j]);               // OK since j>=i;
	++i;
      }
    }
    //cout<<"Comment? "<<i<<":"<<status[i]<<":"<<entry[i]<<std::endl;
    
  }//i
  
  
  n=i;
  
  // ---------- units conversion
  
  for(i=0;i<n;i++)
    {
      area_conversion_factor=-1.;
      if(                                  strstr(Y_units[i],"none")!=NULL)                                  area_conversion_factor=1.0;
      if(strcmp(area_units,"cm2")==NULL && strstr(Y_units[i],"cm2" )!=NULL)                                  area_conversion_factor=1.0;
      if(strcmp(area_units,"cm2")==NULL && strstr(Y_units[i],"m2"  )!=NULL && strstr(Y_units[i],"cm2")==NULL)area_conversion_factor=1.0e-4;
      if(strcmp(area_units,"m2") ==NULL && strstr(Y_units[i],"cm2" )!=NULL)                                  area_conversion_factor=1.0e+4;
      if(strcmp(area_units,"m2") ==NULL && strstr(Y_units[i],"m2"  )!=NULL && strstr(Y_units[i],"cm2")==NULL)area_conversion_factor=1.0;
      
      //cout<<Y_units[i]<<"-> "<<area_units<<" area_conversion_factor="<<area_conversion_factor<<endl;
      if(area_conversion_factor<0.0){cout<<"GCR_data::read: error in converting units"<<endl; exit(0);}
      
      energy_conversion_factor=-1.;
      if(                                    strstr(Y_units[i],"none")!=NULL)   energy_conversion_factor=1.0; 
      
      if(strcmp(energy_units,"MeV")==NULL && strstr(Y_units[i],"/MeV" )!=NULL)  energy_conversion_factor=1.0;
      if(strcmp(energy_units,"MeV")==NULL && strstr(Y_units[i],"/GeV" )!=NULL)  energy_conversion_factor=1.0e-3;
      if(strcmp(energy_units,"MeV")==NULL && strstr(Y_units[i],"MeV/" )!=NULL)  energy_conversion_factor=1.0*pow(E_mean_input[i],-2);
      if(strcmp(energy_units,"MeV")==NULL && strstr(Y_units[i],"GeV/" )!=NULL)  energy_conversion_factor=1.0*pow(E_mean_input[i],-2)*1.e-3;
      
      
      if(strcmp(energy_units,"GeV")==NULL && strstr(Y_units[i],"/GeV" )!=NULL)  energy_conversion_factor=1.0;
      if(strcmp(energy_units,"GeV")==NULL && strstr(Y_units[i],"/MeV" )!=NULL)  energy_conversion_factor=1.0e+3;
      if(strcmp(energy_units,"GeV")==NULL && strstr(Y_units[i],"GeV/" )!=NULL)  energy_conversion_factor=1.0*pow(E_mean_input[i],-2);
      if(strcmp(energy_units,"GeV")==NULL && strstr(Y_units[i],"MeV/" )!=NULL)  energy_conversion_factor=1.0*pow(E_mean_input[i],-2)*1.e+3;
      
      
      // cout<<Y_units[i]<<"-> "<<energy_units<<" energy_conversion_factor="<<energy_conversion_factor<<endl;
      if(energy_conversion_factor<0.0){cout<<"GCR_data::read: error in converting units"<<endl; exit(0);}
      
      value    [i]=    value_input[i]*area_conversion_factor*energy_conversion_factor;
      
      
      //   cout<<"GCR_data:"<<i<<" "<<reference[i]<<" err_type="<<err_type[i]<<" strstr: "<<strncmp(&err_type[i],"a",1)<<" " <<strncmp(&err_type[i],"r",1) <<endl;//AWS20060411
      
      //  if(strstr(&err_type[i],"a")!=NULL)// absolute errors  //AWS20060411
      
      // have to do it this way because err_type is array of characters
      if(err_type[i] == 'a' )// absolute errors //AWS20060411
	{
	  //      cout<<"GCR_data:"<<i<<"  absolute errors: err_type="<<err_type[i]<<endl;//AWS20060411
	  err_minus[i]=err_minus_input[i]*area_conversion_factor*energy_conversion_factor;
	  err_plus [i]=err_plus_input [i]*area_conversion_factor*energy_conversion_factor;
	}
      
      // if(strstr(&err_type[i],"r")!=NULL)// relative errors //AWS20060411
      if(err_type[i] == 'r' )// relative errors   //AWS20060411
	{
	  //     cout<<"GCR_data:"<<i<<" relative errors: err_type="<<err_type[i]<<endl;//AWS20060411
	  err_minus[i]=err_minus_input[i]*value[i];
	  err_plus [i]=err_plus_input [i]*value[i];
	}
      
      
      if(err_plus [i]==0.0)err_plus [i]=err_minus[i]; // if only one value then assume symmetric errors
      if(err_minus[i]==0.0)err_minus[i]=err_plus [i]; 
      
      energy_conversion_factor=-1.;
      if(                                    strstr(X_units[i],"none")!=NULL)  energy_conversion_factor=1.0; 
      if(strcmp(energy_units,"MeV")==NULL && strstr(X_units[i],"MeV" )!=NULL)  energy_conversion_factor=1.0;
      if(strcmp(energy_units,"MeV")==NULL && strstr(X_units[i],"GeV" )!=NULL)  energy_conversion_factor=1.0e+3;
      if(strcmp(energy_units,"GeV")==NULL && strstr(X_units[i],"GeV" )!=NULL)  energy_conversion_factor=1.0;
      if(strcmp(energy_units,"GeV")==NULL && strstr(X_units[i],"MeV" )!=NULL)  energy_conversion_factor=1.0e-3;
      // cout<<X_units[i]<<"-> "<<energy_units<<" energy_conversion_factor="<<energy_conversion_factor<<endl;
      if(energy_conversion_factor<0.0){cout<<"GCR_data::read: error in converting units"<<endl; exit(0);}
            
      E_low [i]=E_low_input [i]*energy_conversion_factor;
      E_high[i]=E_high_input[i]*energy_conversion_factor;
      E_mean[i]=E_mean_input[i]*energy_conversion_factor;
	
    }//i
    
  INFO("Exit");

  //cout<<"<<GCR_data::read"<<endl;
  return 0;
}
//////////////////////////////////////////////
int 
GCR_data::
read(const char* database_file_, 
     char *area_units_, 
     char *energy_units_, 
     char *Y_name_,
     int n_ZA_numerator_select,  
     int *Z_numerator_select,  
     int* A_numerator_select,
     int n_ZA_denominator_select,
     int *Z_denominator_select,
     int* A_denominator_select) {

  INFO("Entry");

  //cout<<">>GCR_data::read(..Z,A ratio) "<<endl;

  int i,j,i_ZA,selected;

  GCR_data::read(database_file_, area_units_, energy_units_); 
  GCR_data::set_plotting();//AWS20060621

  i=0;

  for(j=0;j<n;j++) {

    selected=1;
    
    if(strstr(Y_name[j],Y_name_)!=NULL && n_ZA_numerator_select==n_ZA_numerator    && n_ZA_denominator_select==n_ZA_denominator   )
       {
	 
     for(i_ZA=0;i_ZA< n_ZA_numerator_select;i_ZA++)
       {
	 //cout<<"candidate..Z"<<Z_numerator_select[i_ZA]<<" "<<Z_numerator[j][i_ZA]<<" "<< (Z_numerator_select[i_ZA]==Z_numerator[j][i_ZA])<<endl;
        selected*=        (Z_numerator_select[i_ZA]==Z_numerator[j][i_ZA]);
       }

     for(i_ZA=0;i_ZA< n_ZA_numerator_select;i_ZA++)
       {
	 //cout<<"candidate..A"<<A_numerator_select[i_ZA]<<" "<<A_numerator[j][i_ZA]<<" "<< (A_numerator_select[i_ZA]==A_numerator[j][i_ZA])<<endl;
         selected*=        (A_numerator_select[i_ZA]==A_numerator[j][i_ZA]);
       }
       

    for(i_ZA=0;i_ZA< n_ZA_denominator_select;i_ZA++)
       {
	 // cout<<"candidate..Z"<<Z_denominator_select[i_ZA]<<" "<<Z_denominator[j][i_ZA]<<" "<< (Z_denominator_select[i_ZA]==Z_denominator[j][i_ZA])<<endl;
        selected*=        (Z_denominator_select[i_ZA]==Z_denominator[j][i_ZA]);
       }


    for(i_ZA=0;i_ZA< n_ZA_denominator_select;i_ZA++)
       {
	 //cout<<"candidate..A"<<A_denominator_select[i_ZA]<<" "<<A_denominator[j][i_ZA]<<" "<< (A_denominator_select[i_ZA]==A_denominator[j][i_ZA])<<endl;
        selected*=        (A_denominator_select[i_ZA]==A_denominator[j][i_ZA]);
       }


       }//if
     else {selected=0;}
     //cout<<"j= "<<j<<" selected="<<selected<<endl;

     if(selected==1)
       {
         strcpy(entry [i],    entry[j]);
	 status[i]=status[j];
	 strcpy(reference[i],reference[j]);
	 strcpy(experiment[i],experiment[j]);
	 strcpy(X_units[i],X_units[j]);
	 strcpy(Y_units[i],Y_units[j]);
	 strcpy(Y_name[i],Y_name[j]);

         E_low_input      [i]=     E_low_input[j];
         E_high_input     [i]=    E_high_input[j];
         E_mean_input     [i]=    E_mean_input[j];
         value_input      [i]=     value_input[j];
         err_minus_input  [i]=err_minus_input [j];
         err_plus_input   [i]=err_plus_input  [j];
         err_type         [i]=err_type        [j];

         for (i_ZA=0;i_ZA<n_ZA_numerator;i_ZA++)
         {
          Z_numerator    [i][i_ZA] = Z_numerator    [j][i_ZA];
          A_numerator    [i][i_ZA] = A_numerator    [j][i_ZA];
         }

         for (i_ZA=0;i_ZA<n_ZA_denominator;i_ZA++)
         {
          Z_denominator  [i][i_ZA] = Z_denominator  [j][i_ZA];
          A_denominator  [i][i_ZA] = A_denominator  [j][i_ZA];
         }

         strcpy(comment[i],comment[j]);

         E_low      [i]=     E_low[j];
         E_high     [i]=    E_high[j];
         E_mean     [i]=    E_mean[j];
         value      [i]=     value[j];
         err_minus  [i]=err_minus [j];
         err_plus   [i]=err_plus  [j];


         color      [i]=color     [j]; //AWS20060621
         style      [i]=style     [j]; //AWS20060621
         size       [i]=size      [j]; //AWS20060621

	 i++;
         
       }//if n_ZA_numberator_select

   
   }//j

  n=i;
  return 0;
}
//////////////////////////////////////////////
int GCR_data::set_plotting()
{
  cout<<">>GCR_data::set_plotting"<<endl;

 int i;

 // Here we use the ROOT definitions, but this class is kept independent of ROOT
 // so define the k.. values explicitly.

  //Gtypes.h:enum EColor { kWhite, kBlack, kRed, kGreen, kBlue, kYellow, kMagenta, kCyan };

     // TAttMarker
     //*-*-*-*-*-*-*-*-*-*-*-*Marker Attributes class*-*-*-*-*-*-*-*-*-*-*-*-*-*
     //*                      =======================
     //*  Marker attributes are:
     //*    Marker Color
     //*    Marker style
     //*    Marker Size
     //*
     //*  This class is used (in general by secondary inheritance)
     //*  by many other classes (graphics, histograms).
     //*
     //*  List of the currently supported markers (screen and PostScript)
     //*  ===============================================================
     //*      1 : dot                     kDot
     //*      2 : +                       kPlus
     //*      3 : *                       kStar
     //*      4 : o                       kCircle
     //*      5 : x                       kMultiply
     //*      6 : small scalable dot      kFullDotSmall
     //*      7 : medium scalable dot     kFullDotMedium
     //*      8 : large scalable dot      kFullDotLarge
     //*      9 -->15 : dot
     //*     16 : open triangle down      kOpenTriangleDown
     //*     18 : full cross              kFullCross
     //*     20 : full circle             kFullCircle
     //*     21 : full square             kFullSquare
     //*     22 : full triangle up        kFullTriangleUp
     //*     23 : full triangle down      kFullTriangleDown
     //*     24 : open circle             kOpenCircle
     //*     25 : open square             kOpenSquare
     //*     26 : open triangle up        kOpenTriangleUp
     //*     27 : open diamond            kOpenDiamond
     //*     28 : open cross              kOpenCross
     //*     29 : open star               kOpenStar
     //*     30 : full star               kFullStar
     //*

   int  kWhite=0, kBlack=1, kRed=2, kGreen=3, kBlue=4, kYellow=5, kMagenta=6, kCyan=7;

   int  kDot              =1,
        kPlus             =2,
        kStar             =3,
        kCircle           =4,
        kMultiply         =5,
        kFullDotSmall     =6,
        kFullDotMedium    =7,
        kFullDotLarge     =8,
     //      9 -->15 : dot
        kOpenTriangleDown=16,
        kFullCross       =18,
        kFullCircle      =20,
        kFullSquare      =21,
        kFullTriangleUp  =22,
        kFullTriangleDown=23,
        kOpenCircle      =24,
        kOpenSquare      =25,
        kOpenTriangleUp  =26,
        kOpenDiamond     =27,
        kOpenCross       =28,
        kOpenStar        =29,
        kFullStar        =30;

  

  


  // all experiments in GCR_data_1.dat
  for(i=0;i<n;i++)
  {


   // default marker
  color[i]=kBlack;                                          
  style[i]=kOpenCircle;                                   
  size [i]=1.0;

  if(strncmp(experiment[i],"AMS01"   ,5)==0){ color[i]=kRed ;   style[i]=kOpenCircle;  } 
  if(strncmp(experiment[i],"JACEE"   ,5)==0){ color[i]=kRed ;   style[i]=kOpenSquare;  }                  
  if(strncmp(experiment[i],"BESS"    ,4)==0){ color[i]=kRed ;   style[i]=kOpenTriangleUp  ;}            
  if(strncmp(experiment[i],"IMP"     ,3)==0){ color[i]=kRed ;   style[i]=kFullDotLarge;}//kOpenTriangleDown: produces dot
  if(strncmp(experiment[i],"CAPRICE" ,7)==0){ color[i]=kRed  ;  style[i]=kOpenStar  ;  } 
  if(strncmp(experiment[i],"RUNJOB"  ,6)==0){ color[i]=kRed ;   style[i]=kOpenCross ;  } 

  if(strncmp(experiment[i],"HEAT"    ,4)==0){ color[i]=kBlue ;  style[i]=kFullCircle;  } 
  if(strncmp(experiment[i],"HEAO3"   ,5)==0){ color[i]=kBlue;   style[i]=kFullSquare;  }
  if(strncmp(experiment[i],"SOKOL"   ,5)==0){ color[i]=kBlue;   style[i]=kFullTriangleUp  ; }
  if(strncmp(experiment[i],"SANRIKU" ,7)==0){ color[i]=kBlue;   style[i]=kFullTriangleDown;}
  if(strncmp(experiment[i],"IMAX"    ,4)==0){ color[i]=kBlue;   style[i]=kFullStar ;  } 
  if(strncmp(experiment[i],"CRN"     ,3)==0){ color[i]=kBlue;   style[i]=kOpenDiamond; } // kFullCross:no symbol plotted
  if(strncmp(experiment[i],"ATIC-2"  ,6)==0){ color[i]=kCyan;   style[i]=kOpenTriangleUp  ; }
  if(strncmp(experiment[i],"ATIC-1-2",8)==0){ color[i]=kCyan;   style[i]=kOpenDiamond   ;   } //AWS20081124
  if(strncmp(experiment[i],"BETS"    ,4)==0){ color[i]=kRed ;   style[i]=kMultiply      ;   } //AWS20081128
  if(strncmp(experiment[i],"PPB-BETS",7)==0){ color[i]=kRed ;   style[i]=kMultiply      ;   } //AWS20081124
  if(strncmp(experiment[i],"HESS"    ,4)==0){ color[i]=kGreen;  style[i]=kOpenSquare    ;   } //AWS20081126

  if(strncmp(experiment[i],"ACE"     ,3)==0){ color[i]=kCyan;   style[i]=kPlus      ;  }
  if(strncmp(experiment[i],"Voyager" ,7)==0){ color[i]=kMagenta;style[i]=kMultiply  ;  } 
  if(strncmp(experiment[i],"MUBEE"   ,5)==0){ color[i]=kGreen;  style[i]=kOpenDiamond     ; } 

  cout<<" GCR_set_plotting:  experiment: "<<  experiment[i]<<" colour= "<<color[i]<<" style="<<style[i]
      <<" E_mean="<<E_mean[i]<<" plot value="<<  value [i]*pow(E_mean[i],2) <<endl;        

  }



  return 0;
}
//////////////////////////////////////////////
int GCR_data::print()
{
  cout<<">>GCR_data::print"<<endl;

 int i,i_ZA;

 cout<<"database_file="<<database_file<<endl; 
 cout<<"area_units   ="<<area_units   <<endl;
 cout<<"energy_units ="<<energy_units <<endl;

 cout<<"number of entries ="<<n<<"   (# commented eliminated)"<<endl;

 cout<<"selected original database entries  :"; cout<<endl;
 for (i=0;i<n;i++)cout<<entry[i]<<endl; cout<<endl;

 cout<<"selected values   :"; cout<<endl;

 for (i=0;i<n;i++)
{
  cout<<status[i]<<" "<<reference[i]<<" "<<experiment[i]<<" "<<X_units[i]<<" "<<Y_units[i]<<" "<<Y_name[i];
  cout<<" E_low_input="<<E_low_input[i]<<" ";
  cout<<" E_high_input="<<E_high_input[i]<<" ";
  cout<<" E_mean_input="<<E_mean_input[i]<<" ";
  cout<<"  value_input="<< value_input[i]<<" ";
  cout<<"err_minus_input="<< err_minus_input[i]<<" ";
  cout<<"err_plus_input="<< err_plus_input[i]<<" ";
  cout<<"err_type="      << err_type      [i]<<" ";
  cout<<"Z.A numerator =";
  for (i_ZA=0;i_ZA<n_ZA_numerator;i_ZA++)
    {
  cout   <<Z_numerator  [i][i_ZA] <<"."
         <<A_numerator  [i][i_ZA] <<" ";
    }

  cout<<"Z.A denominator =";
  for (i_ZA=0;i_ZA<n_ZA_denominator;i_ZA++)
    {
  cout   <<Z_denominator[i][i_ZA] <<"."
         <<A_denominator[i][i_ZA] <<" ";
    }


  cout<<"    E_low="<<     E_low[i]<<" ";
  cout<<"   E_high="<<    E_high[i]<<" ";
  cout<<"   E_mean="<<    E_mean[i]<<" ";
  cout<<"    value="<<     value[i]<<" ";
  cout<<"err_minus="<< err_minus[i]<<" ";
  cout<<" err_plus="<<  err_plus[i]<<" ";
  cout<<endl;
  cout<<"Comment:"<<comment[i]<<endl;

  //cout<<"........"<<endl;

} 

 cout<<"----------------"<<endl;
 cout<<"<<GCR_data::print"<<endl;
return 0;
}

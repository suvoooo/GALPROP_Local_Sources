
//**.****|****.****|****.****|****.****|****.****|****.****|****.****|****.****|
// * print_BC.cc *                                 galprop package * 4/14/2000 
//**"****!****"****!****"****!****"****!****"****!****"****!****"****!****"****|
using namespace std;//AWS20050624
#include<iostream>


#include"galprop_classes.h"

int Galprop::print_BC()
{
int stat;

cout<<" >>>> print_BC"<<endl;
stat=0;

int iC12=-1,iB10=-1,iB11=-1;
  
  
  for (int i=0; i<n_species; i++){


 

    if(gcr[i].A==12 && gcr[i].Z==6) iC12=i;
    if(gcr[i].A==10 && gcr[i].Z==5) iB10=i;
    if(gcr[i].A==11 && gcr[i].Z==5) iB11=i;
  }

  // run may not include C and B
  if(iC12==-1 || iB10==-1 || iB11==-1){cout<<"no C or no B\n<<<print_BC\n";return stat;}

    /*
 if(galdef.n_spatial_dimensions==2){

   
   for(int ir=0;ir<gcr[i].n_rgrid;ir++){
   for(int iz=0;iz<gcr[i].n_zgrid;iz++){
   for(int ip=0;ip<gcr[i].n_pgrid;ip++){
 
     // if(galdef.source_specification==2 && iz==gcr[i].n_zgrid/2        )
       gcr[i].primary_source_function.d2[ir][iz].s[ip]= gcr[i].primary_abundance  ;
   }
   }
   }
 }

 if(galdef.n_spatial_dimensions==3){

 
   for(int ix=0;ix<gcr[i].n_xgrid;ix++){
   for(int iy=0;iy<gcr[i].n_ygrid;iy++){
   for(int iz=0;iz<gcr[i].n_zgrid;iz++){
   for(int ip=0;ip<gcr[i].n_pgrid;ip++){
  

 
   }
   }
   }
   }
 }

    */

 

  cout<<"iC12 iB10 iB11 "<< iC12<<" " <<iB10<<" " <<iB11<<endl;
   gcr[iC12].print();
   gcr[iB10].print();
   gcr[iB11].print();

    Distribution BC;
    if(galdef.n_spatial_dimensions==2) BC.init(gcr[0].n_rgrid,gcr[0].n_zgrid,gcr[0].n_pgrid);
    if(galdef.n_spatial_dimensions==3) BC.init(gcr[0].n_xgrid,gcr[0].n_ygrid,gcr[0].n_zgrid,gcr[0].n_pgrid);

       BC=(gcr[iB10].cr_density+gcr[iB11].cr_density)/ gcr[iC12].cr_density ;
  
       /*
BC+=gcr[iB10].cr_density;
 
BC+=gcr[iB11].cr_density;
 
BC/=gcr[iC12].cr_density;
       */
    cout<<"BC.print"<<endl;
    BC.print();

		BC.delete_array();  //Gulli20070810



cout<<" <<<< print_BC"<<endl;
return stat;
}

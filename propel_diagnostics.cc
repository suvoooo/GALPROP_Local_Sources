#include"ErrorLogger.h"
#include"galprop_classes.h"
#include"galprop_internal.h"

#include<sstream>

// variables maintained over calls of propel_diagnostics  AWS20101231
namespace propel_diagnostics
{
  int number_of_calls = 0;
  int    total_steps  = 0;
  double total_time   = 0.;
  double last_max     = 0;
  Distribution last_cr_density,diff_cr_density;
}
//**.****|****.****|****.****|****.****|****.****|****.****|****.****|****.****|
// * propel_diagnostics.cc *                       galprop package * 4/14/2000 
//**"****!****"****!****"****!****"****!****"****!****"****!****"****!****"****|
using namespace std;//AWS20050624
using namespace propel_diagnostics; //AWS20101231


int Galprop::propel_diagnostics()   //AWS20110113
{
  INFO("Entry");
  INFO("initializing variables for new species");
  
  number_of_calls = 0;
    total_steps   = 0;
    total_time    = 0.;

    return 0;

}


///////////////////////////////////////////////////////////////////////////////////
int Galprop::propel_diagnostics
(Particle &particle,
Distribution &alpha1_r,
Distribution &alpha1_z,
Distribution &alpha1_p,
Distribution &alpha2_r,
Distribution &alpha2_z,
Distribution &alpha2_p,
Distribution &alpha3_r,
Distribution &alpha3_z,
Distribution &alpha3_p,
Distribution &total_source_function,
double dt
)
{




int ir,iz,ip;


 if(galdef.control_diagnostics>0) //AWS20101227
 { 
  INFO("Entry");
  INFO("generating diagnostics for 2D");
 }
 
 if(number_of_calls==0)
 {
  last_cr_density=particle.cr_density;   //AWS20101231
  diff_cr_density=particle.cr_density;
//last_cr_density.print(); 
 }

 if(number_of_calls>0)
 {
 diff_cr_density = particle.cr_density;
 last_cr_density*= -1.0; // since no "-" operator available
 diff_cr_density+=     last_cr_density; 
 diff_cr_density/= particle.cr_density;
 if(galdef.control_diagnostics==4)diff_cr_density.print();
 }

 number_of_calls += 1;                                         //AWS20101231
 total_steps     += galdef.timestep_diagnostics;               //AWS20101231 
 total_time      += galdef.timestep_diagnostics * dt/year2sec; //AWS20101231 
 

Distribution dphidt         (particle.n_rgrid, particle. n_zgrid,particle. n_pgrid);
 
Distribution timescale      (particle.n_rgrid, particle. n_zgrid,particle. n_pgrid);

  

   // r propagation
 if(galdef.prop_r==1){

 //Gulli20070821
#pragma omp parallel for schedule(dynamic) default(shared) private(iz,ir,ip) 
 for(iz=0;iz<particle.n_zgrid;iz++){
 for(ip=0;ip<particle.n_pgrid;ip++){

   for(ir=1;ir<particle.n_rgrid-1;ir++) {

dphidt.d2[ir][iz].s[ip]+=alpha1_r.d2[ir][iz].s[ip] *particle.cr_density.d2[ir-1][iz  ].s[ip  ];
dphidt.d2[ir][iz].s[ip]-=alpha2_r.d2[ir][iz].s[ip] *particle.cr_density.d2[ir  ][iz  ].s[ip  ];
dphidt.d2[ir][iz].s[ip]+=alpha3_r.d2[ir][iz].s[ip] *particle.cr_density.d2[ir+1][iz  ].s[ip  ];

   }                  
    
 }//ip
 }//iz
 }//prop_r


 
   // z propagation
 if(galdef.prop_z==1){

 //Gulli20070821
#pragma omp parallel for schedule(dynamic) default(shared) private(iz,ir,ip) 
 for(ir=0;ir<particle.n_rgrid;ir++){
 for(ip=0;ip<particle.n_pgrid;ip++){

   for(iz=1;iz<particle.n_zgrid-1;iz++) {

dphidt.d2[ir][iz].s[ip]+=alpha1_z.d2[ir][iz].s[ip] *particle.cr_density.d2[ir  ][iz-1].s[ip  ];
dphidt.d2[ir][iz].s[ip]-=alpha2_z.d2[ir][iz].s[ip] *particle.cr_density.d2[ir  ][iz  ].s[ip  ];
dphidt.d2[ir][iz].s[ip]+=alpha3_z.d2[ir][iz].s[ip] *particle.cr_density.d2[ir  ][iz+1].s[ip  ];

   } 
 

 }//ip
 }//ir
 }//prop_z



  // p propagation

 if(galdef.prop_p==1){ 

 //Gulli20070821
#pragma omp parallel for schedule(dynamic) default(shared) private(iz,ir,ip) 
 for(ir=0;ir<particle.n_rgrid;ir++){
 for(iz=0;iz<particle.n_zgrid;iz++){

   for(ip=1;ip<particle.n_pgrid-1;ip++) {
dphidt.d2[ir][iz].s[ip]+=alpha1_p.d2[ir][iz].s[ip] *particle.cr_density.d2[ir  ][iz  ].s[ip-1];
dphidt.d2[ir][iz].s[ip]-=alpha2_p.d2[ir][iz].s[ip] *particle.cr_density.d2[ir  ][iz  ].s[ip  ];
dphidt.d2[ir][iz].s[ip]+=alpha3_p.d2[ir][iz].s[ip] *particle.cr_density.d2[ir  ][iz  ].s[ip+1];
  
   } 
 
 }//ir
 }//ix
 }//prop_p

 
dphidt/=dt; // since alpha is passed as alpha*dt

 if (galdef.control_diagnostics>=3){  //AWS20101221             
    INFO("dphidt without source function");
dphidt.print();

INFO("               source function");
total_source_function.print();  
 }

 
dphidt+=   total_source_function;     
      
 if (galdef.control_diagnostics>=2){
    INFO("dphidt with    source function");
dphidt.print();
 
 }

timescale=particle.cr_density;
timescale/=dphidt;
timescale/=year2sec;

 //Gulli20070821
#pragma omp parallel for schedule(dynamic) default(shared) private(iz,ir,ip) 
for(ir=0;ir<particle.n_rgrid;ir++)
for(iz=0;iz<particle.n_zgrid;iz++)
for(ip=0;ip<particle.n_pgrid;ip++) timescale.d2[ir][iz].s[ip]=fabs(timescale.d2[ir][iz].s[ip]);

double timescale_min=1.e30;
double timescale_max=0.0;

 for(ir=1;ir<particle.n_rgrid-1;ir++){  // edges not included in search for min,max
 
 for(iz=1;iz<particle.n_zgrid-1;iz++){

 for(ip=1;ip<particle.n_pgrid-1;ip++){

if(timescale.d2[ir]    [iz].s[ip]<timescale_min)timescale_min=timescale.d2[ir]    [iz].s[ip];
if(timescale.d2[ir]    [iz].s[ip]>timescale_max)timescale_max=timescale.d2[ir]    [iz].s[ip];

 }
 }
 }
 


double units=1.0e6;
ostringstream buf;

/* AWS20101221 since print(units) gives wrong output
 buf<<"timescale  ("<<units<<" years)            ";
 INFO(buf.str());
if(galdef.control_diagnostics>=1)timescale.print(units);
*/




 if(galdef.control_diagnostics>=1)                   //AWS20101227
 { 
  buf<<"timescale  (years)            ";             //AWS20101221
  INFO(buf.str());                                   //AWS20101221 
  timescale.print();//AWS20101221
 }

buf.str("");

 buf<<"call #"<<number_of_calls<<" dt="<<dt/year2sec<<" total steps="<<total_steps<<" total time="<<total_time ; //AWS20101231
 buf<<" alpha timescale min max ="<<timescale_min<<"  "<<timescale_max;
 //buf << " yrs CR max = "<< particle.cr_density.max(); //AWS20101229
 // buf<< " diff from last max= "<< particle.cr_density.max()-last_max; //AWS20101231
 buf<< " max frac diff wrt last CR="<<diff_cr_density.max();        //AWS20101231
 buf<< " gives min timescale="<< dt/year2sec*galdef.timestep_diagnostics / diff_cr_density.max(); //AWS20110101
 buf<< " yrs";

INFO(buf.str());

last_max=particle.cr_density.max();
last_cr_density=particle.cr_density; //AWS20101231

dphidt.delete_array();
timescale.delete_array();

 return 0;}

/////////////////////////////////////////////////////////////
/////////////////////   3D    ///////////////////////////////
///////////////////////////////////////////////////////////// 

int Galprop::propel_diagnostics
(Particle &particle,
Distribution &alpha1_x,
Distribution &alpha1_y,
Distribution &alpha1_z,
Distribution &alpha1_p,
Distribution &alpha2_x,
Distribution &alpha2_y,
Distribution &alpha2_z,
Distribution &alpha2_p,
Distribution &alpha3_x,
Distribution &alpha3_y,
Distribution &alpha3_z,
Distribution &alpha3_p,
Distribution &total_source_function,
double dt
)
{

  if(galdef.control_diagnostics>0)//AWS20101227
  {
   INFO("Entry");
   INFO("generating diagnostics  for 3D");
  }

 if(number_of_calls==0)
 {
  last_cr_density=particle.cr_density;   //AWS20101231
  diff_cr_density=particle.cr_density;
//last_cr_density.print(); 
 }

 if(number_of_calls>0)
 {
 diff_cr_density = particle.cr_density;
 last_cr_density*= -1.0; // since no "-" operator available
 diff_cr_density+=     last_cr_density; 
 diff_cr_density/= particle.cr_density;
 if(galdef.control_diagnostics==4)diff_cr_density.print();
 }

 number_of_calls += 1;                                         //AWS20101231
 total_steps     += galdef.timestep_diagnostics;               //AWS20101231 
 total_time      += galdef.timestep_diagnostics * dt/year2sec; //AWS20101231 


Distribution dphidt         (particle.n_xgrid,particle.n_ygrid,particle. n_zgrid,particle. n_pgrid);
 
Distribution timescale      (particle.n_xgrid,particle.n_ygrid,particle. n_zgrid,particle. n_pgrid);

 
  int ix,iy,iz,ip;

   // x propagation

 if(galdef.prop_x==1){

  //Gulli20070821
#pragma omp parallel for schedule(dynamic) default(shared) private(iz,ix,iy,ip) 
 for(ip=0;ip<particle.n_pgrid;ip++){
 for(iz=0;iz<particle.n_zgrid;iz++){
 for(iy=0;iy<particle.n_ygrid;iy++){

   for(ix=1;ix<particle.n_xgrid-1;ix++) {

dphidt.d3[ix][iy][iz].s[ip]+=alpha1_x.d3[ix][iy][iz].s[ip] *particle.cr_density.d3[ix-1][iy][iz  ].s[ip  ];
dphidt.d3[ix][iy][iz].s[ip]-=alpha2_x.d3[ix][iy][iz].s[ip] *particle.cr_density.d3[ix  ][iy][iz  ].s[ip  ];
dphidt.d3[ix][iy][iz].s[ip]+=alpha3_x.d3[ix][iy][iz].s[ip] *particle.cr_density.d3[ix+1][iy][iz  ].s[ip  ];
   
   }//ix   
 }//iy
 }//iz
 }//ip
 }// prop_x

   // y propagation

 if(galdef.prop_y==1){

  //Gulli20070821
#pragma omp parallel for schedule(dynamic) default(shared) private(iz,ix,iy,ip) 
 for(ip=0;ip<particle.n_pgrid;ip++){
 for(iz=0;iz<particle.n_zgrid;iz++){
 for(ix=0;ix<particle.n_xgrid;ix++){
 
   for(iy=1;iy<particle.n_ygrid-1;iy++) {

dphidt.d3[ix][iy][iz].s[ip]+=alpha1_y.d3[ix][iy][iz].s[ip] *particle.cr_density.d3[ix  ][iy-1][iz  ].s[ip  ];
dphidt.d3[ix][iy][iz].s[ip]-=alpha2_y.d3[ix][iy][iz].s[ip] *particle.cr_density.d3[ix  ][iy  ][iz  ].s[ip  ];
dphidt.d3[ix][iy][iz].s[ip]+=alpha3_y.d3[ix][iy][iz].s[ip] *particle.cr_density.d3[ix  ][iy+1][iz  ].s[ip  ];
   
   }//iy
   

 }//iy
 }//iz
 }//ip
 }// prop_y



   // z propagation

 if(galdef.prop_z==1){

  //Gulli20070821
#pragma omp parallel for schedule(dynamic) default(shared) private(iz,ix,iy,ip) 
 for(ip=0;ip<particle.n_pgrid;ip++){
 for(ix=0;ix<particle.n_xgrid;ix++){
 for(iy=0;iy<particle.n_ygrid;iy++){

   for(iz=1;iz<particle.n_zgrid-1;iz++) {

dphidt.d3[ix][iy][iz].s[ip]+=alpha1_z.d3[ix][iy][iz].s[ip] *particle.cr_density.d3[ix  ][iy  ][iz-1].s[ip  ];
dphidt.d3[ix][iy][iz].s[ip]-=alpha2_z.d3[ix][iy][iz].s[ip] *particle.cr_density.d3[ix  ][iy  ][iz  ].s[ip  ];
dphidt.d3[ix][iy][iz].s[ip]+=alpha3_z.d3[ix][iy][iz].s[ip] *particle.cr_density.d3[ix  ][iy  ][iz+1].s[ip  ];
      
   } 
      

 }//iy
 }//iz
 }//ip
 }//prop_z

  // p propagation

 if(galdef.prop_p==1){ 

  //Gulli20070821
#pragma omp parallel for schedule(dynamic) default(shared) private(iz,ix,iy,ip) 
 for(ix=0;ix<particle.n_xgrid;ix++){
 for(iy=0;iy<particle.n_ygrid;iy++){
 for(iz=0;iz<particle.n_zgrid;iz++){

   for(ip=1;ip<particle.n_pgrid-1;ip++) {

dphidt.d3[ix][iy][iz].s[ip]+=alpha1_p.d3[ix][iy][iz].s[ip] *particle.cr_density.d3[ix  ][iy  ][iz  ].s[ip-1];
dphidt.d3[ix][iy][iz].s[ip]-=alpha2_p.d3[ix][iy][iz].s[ip] *particle.cr_density.d3[ix  ][iy  ][iz  ].s[ip  ];
dphidt.d3[ix][iy][iz].s[ip]+=alpha3_p.d3[ix][iy][iz].s[ip] *particle.cr_density.d3[ix  ][iy  ][iz  ].s[ip+1];
      
   
   } 
 

 }//iz
 }//iy
 }//ix
 }//prop_p




dphidt/=dt; // since alpha is passed as alpha*dt

 if (galdef.control_diagnostics>=3){            //AWS20101221
    INFO("dphidt without source function");
dphidt.print();

INFO("               source function");
total_source_function.print();  
 }

 
dphidt+=   total_source_function;     
      
 if (galdef.control_diagnostics>=2){            //AWS20101221
    INFO("dphidt with    source function");
dphidt.print();
 
 }

timescale=particle.cr_density;
timescale/=dphidt;
timescale/=year2sec;

  //Gulli20070821
#pragma omp parallel for schedule(dynamic) default(shared) private(iz,ix,iy,ip) 
 for(ix=0;ix<particle.n_xgrid  ;ix++)   
 for(iy=0;iy<particle.n_ygrid  ;iy++)
 for(iz=0;iz<particle.n_zgrid  ;iz++)
 for(ip=0;ip<particle.n_pgrid  ;ip++) timescale.d3[ix][iy][iz].s[ip]=fabs(timescale.d3[ix][iy][iz].s[ip]);
 
double timescale_min=1.e30;
double timescale_max=0.0;



 for(ix=1;ix<particle.n_xgrid-1;ix++){  // edges not included in search for min,max
 for(iy=1;iy<particle.n_ygrid-1;iy++){
 for(iz=1;iz<particle.n_zgrid-1;iz++){
 for(ip=1;ip<particle.n_pgrid-1;ip++) {

  if(timescale.d3[ix][iy][iz].s[ip]<timescale_min)timescale_min=timescale.d3[ix][iy][iz].s[ip];
  if(timescale.d3[ix][iy][iz].s[ip]>timescale_max)timescale_max=timescale.d3[ix][iy][iz].s[ip];

   }
 }
 }
 }


double units=1.0e6;
ostringstream buf;
 /* since print(units) gives wrong output
 buf<<"timescale  ("<<units<<" years)            ";
 INFO(buf.str());
if(galdef.control_diagnostics>=1)timescale.print(units);
 */


 if(galdef.control_diagnostics>=1)                   //AWS20101227
 {
  buf<<"timescale  (years)            ";             //AWS20101221
  INFO(buf.str());                                   //AWS20101221 
  timescale.print();//AWS20101221
 }

buf.str("");


 buf<<"call #"<<number_of_calls<<" dt = "<<dt/year2sec<<" total steps="<<total_steps<<" total time="<<total_time ; //AWS20101231
 buf<<"  alpha timescale min max  = "<<timescale_min<<"  "<<timescale_max << " yrs. CR max = "<< particle.cr_density.max(); //AWS20101229
 // buf<< " diff from last max= "<< particle.cr_density.max()-last_max; //AWS20101231
 buf<< " max frac diff from last CR="<<diff_cr_density.max();        //AWS20101231
 buf<< " corresponding min timescale="<< dt/year2sec*galdef.timestep_diagnostics / diff_cr_density.max(); //AWS20110101
 buf<< " yrs";

INFO(buf.str());


dphidt.delete_array();
timescale.delete_array();

last_max=particle.cr_density.max();
last_cr_density=particle.cr_density; //AWS20101231

//cout<<"\n<<<<propel_diagnostics"<<endl;

return 0;
}







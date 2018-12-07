void B_field_3D_model
(char* name, double *parameters,
 double x, double y, double z,
 int options,
 double &Breg,  double &Bregx, double &Bregy, double &Bregz,
 double &Bran,  double &Branx, double &Brany, double &Branz,
 int debug=0);

double Bperp(double x,double y,double z, double Bx,double By,double Bz, double x0,double y0,double z0);




double B_field_3D_model_tot
(char* name, double *parameters,
 double x, double y, double z,
 int debug=0 );



double B_field_3D_model_tot
(char* name, double *parameters,
 double R, double z,
 int debug=0 );

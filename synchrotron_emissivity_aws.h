double synchrotron_emissivity_aws(double gamma, double nu, double x, double y, double z,  char* name, double *parameters, int debug=0);

double synchrotron_emissivity_aws(double gamma, double nu, double R, double z, char* name, double *parameters, int debug=0);


double synchrotron_emissivity_aws(double gamma, double nu, double x, double y, double z,  char* name, double *parameters, double &I, double &Q, double &U, int debug=0);//AWS20100707

double synchrotron_emissivity_aws(double gamma, double nu, double R, double z, char* name, double *parameters, double &I, double &Q, double &U, int debug=0);//AWS20100707


double synchrotron_emissivity_aws(double gamma, double nu, double x, double y, double z,  int galdef_B_field_model,int debug=0);

double synchrotron_emissivity_aws(double gamma, double nu, double R, double z, int galdef_B_field_model, int debug=0);

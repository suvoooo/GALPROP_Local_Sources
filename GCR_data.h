#ifndef GCR_DATA_H
#define GCR_DATA_H
class GCR_data
{

 public:

  char database_file[400];
  int n; // number of entries
  int N; // number of entries read from database, including comments
  
  // as read
  char **entry;     // complete input string
  char  *status;    // # etc use, dont use etc.
  char **reference;
  char **experiment;
  char **Y_name;   // e.g. flux
  double *E_low_input;    
  double *E_high_input; 
  double *E_mean_input;
  double *value_input;
  double *err_minus_input;
  double *err_plus_input; 

  char   *err_type; // a=abs r=relative

  int **Z_numerator;
  int **A_numerator;
  int **Z_denominator;
  int **A_denominator;

  char **comment;
  char **X_units;
  char **Y_units;

  // converted and derived quantities

  char energy_units[4];// "MeV"|"GeV"
  char area_units  [4];  // "cm2"|"m2"

  double *E_low;    // in energy units
  double *E_high;   
  double *E_mean;

  double *value;
  double *err_minus;
  double *err_plus; 


  int n_ZA_numerator;
  int n_ZA_denominator;

  
  // information for plotting: user-defined      // AWS20060621
  int    *color;                                 // AWS20060621
  int    *style;                                 // AWS20060621
  double *size ;                                 // AWS20060621

  GCR_data();
  ~GCR_data();

  int read(const char* database_file_, char *area_units, char *energy_units_); 

  int read(const char *database_file_, char *area_units_, char *energy_units_,
           char *Y_name_,
	   int n_ZA_numerator_select,  int *Z_numerator_select,  int* A_numerator_select,
           int n_ZA_denominator_select,int *Z_denominator_select,int* A_denominator_select);

  int set_plotting();                            // AWS20060621

  int print();
};
#endif

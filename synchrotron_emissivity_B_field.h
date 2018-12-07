double synchrotron_emissivity_B_field( double gamma, double nu,
                                       char *name, double *parameters,
                                       double x, double y, double z,
                                       int options,
                                       double x0, double y0, double z0,
                                       double     &synch_emissivity_reg, double &synch_emissivity_par, double &synch_emissivity_perp, double &synch_emissivity_random,
                                       int debug );



double synchrotron_emissivity_B_field( double gamma, double nu,
                                       char *name, double *parameters,
                                       double x, double y, double z,
                                       int options,
                                       double x0, double y0, double z0,
                                       double     &synch_emissivity_reg, double &synch_emissivity_par, double &synch_emissivity_perp, double &synch_emissivity_random,
                                       double     &I, double     &Q, double &U,
                                       int debug ); //AWS20100706


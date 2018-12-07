
//**.****|****.****|****.****|****.****|****.****|****.****|****.****|****.****|
// * nuc_package.cc *                              galprop package * 2001/08/16
//**"****!****"****!****"****!****"****!****"****!****"****!****"****!****"****|

//**.****|****.****|****.****|****.****|****.****|****.****|****.****|****.****|
// This file contains 4 routines on c++ which read the nuclear data from files,
// evaluate the isotopic cross sections and take care of the unstable nuclei
// decay channels. A sample driver program is in the bottom of this file.
// The isotope_cs routine uses Webber's and Tsao & Silberberg fortran codes,
// wsigma and yieldx.
//                            ### Igor Moskalenko, NASA/GSFC ### 2001/08/16 ###
//
// Description of routines
// ^^^^^^^^^^^^^^^^^^^^^^^
//1 void read_nucdata() - reads 3D data arrays from files, it is organized
//1 as following:
//1 a constant N_DATA_FILES defines the total number of data-files to read;
//1 an array of pointers data_filename[N_DATA_FILES] specifies the file names;
//1 the dimensions of arrays in files are read into n_data[3][N_DATA_FILES];
//1 the correspondence between data files and routines is established by  
//1 file_no[N_DATA_FILES] = {0,1,2}, where 0=isotope_cs.dat, 1=nucdata.dat etc.;
//1 finally an array of pointers data_file[N_DATA_FILES] indicates the areas
//1 allocated for the data arrays which store the data from files.
//1 The DATA FILES to read by the read_nucdata are organized as following:
//1 in the beginning there are comments with basic description (any line number),
//1 which follows by 3 integer numbers - array dimensions, which follows by
//1 the data lines with some lines commented.
//1 The comments should contain any simbol in the 1st column (not space !), 
//1 space in the 1st column will mean the data line.
//
//2 float nucdata
//2   (int ksp,int iz,int ia,int K_electron,int izf,int iaf,int* izl,int* ial,float* To) -
//2 for any given primary nucleus (iz,ia) and any stable (or long-lived) 
//2 final nucleus (izf,iaf) calculates the decay chain and returns branching ratio.
//2    input:
//2 ksp =0 - network of decays compiled from [NDS]; =1 from [GM87];
//2 iz,ia - primary charge and atomic number (iz<29, ia<65);
//2 K_electron =0 for the naked (iz,ia) nucleus; =1 for H-like atoms.
//2 izf,iaf - final charge and atomic number;
//2 uses file_no[1] as indicator of the data array.
//2    output:
//2 izl,ial - return charge and atomic number of an INTERMEDIATE long-lived 
//2           isotope (if any); otherwise (izl,ial)=(0,0);
//2 To  - if (izl,ial)#0, returns half-lifetime of (izl,ial) in sec; 
//2       if (izl,ial)=0, but primary (iz,ia) is a long-lived isotope,
//2          To gives its half-lifetime;
//2       =0 otherwise;
//2 branching ratio is returned by nucdata: 
//2 if (iz,ia)=(izf,iaf), B =1 for stables and long-lived nuclei; =0 otherwise.
//2    comments:
//2 For each Z and A boundary nuclei are chosen (on the left and right), so that
//2 outside of this area all decay treated as proton- or neutron-emission.
//2 Inside this area all decays treated as beta(+/-) decays unless there are
//2 special decay channels - see file nucdata.dat. The program allows up
//2 to 4 generations in any of the special decay channels (up to 3^4 nuclei
//2 in the final state !), beta,p,n-decays - no limits, and returns the final
//2 branching for specified nucleus and To for intermediate long-lived isotope 
//2 if exists.
//2 The nucdata.dat file contains zi.ai,zf1.af1,br1,zf2.af2,br2,zf3.af3,br3:
//2 decay channels with their branchings in order of increasing ai; and zi 
//2 for the same ai.
//2 Rule - secondary zf.af can appear ONLY BELOW the line where it was primary !
//2 Any character in the 1st column can be used to comment the line.
//2 Data on the lifetime T1/2 of long-live isotopes in cosmic rays include
//2 reanalysis of their decay probabilities. This is particular important for
//2 secondary K-electron capture isotopes which is suppressed in CR where 
//2 the main channel is to be beta(+/-) decay.
//2    references:
//2 [B76]  Berenyi, D. et al. 1976, NPA 256, 87
//2 [F99]  Fisker,J.L., Martinez-Pinedo,G., Langamke,K. 1999, Eur.Phys.J. A5, 229
//2 [GM87] Garcia-Munoz M. et al. 1987, ApJS 64, 269 (non-beta decays & branch.)
//2 [M98]  Martinez-Pinedo G., Vogel P. 1998, PRL 81, 281
//2 [NDS]  Nuclear Data Sheets
//2 [ToI]  Table of Isotopes, 8th Ed., by R.B.Firestone etal.(J.Wiley & Sons,Inc.),1996. 
//2 [W98]  Wuosmaa A.H. et al. 1998, PRL 80, 2085 (54Mn[b+ & b-] half lifetime)  
//
//3 float isotope_cs(float emev,int iz,int ia,int izf,int iaf,int kopt,int* info)
//3 For any given primary nucleus (iz,ia) and any final nucleus (izf,iaf) 
//3 returns the cross section (mb) of the reaction p+(zi,ai) -> (zf,af) +X.
//3    input:
//3 emev  - energy of the primary nucleus in MeV/nucleon;
//3 iz,ia - primary charge and atomic number;
//3 izf,iaf - final charge and atomic number/ izf=0 is alloyed, gives fits of
//3           isobaric cross sections (only with kopt = 11, 21);
//3 kopt  =0 - uses best alghorithm described in comments below (not recommended);
//3       =1 - forces to use Webber'93 code (no renormalization etc.);
//3       =2 - forces to use TS00 code (no renormalization etc.);
//3       =3 - forces to use a const cross section fitted to the data.
//3       =10- forces to use Webber'93 code (renormalized if data exists);
//3       =11- forces to use cross section fit if exists (otherwise equiv. 10);
//3       =12- forces to use a numerical table if exists (otherwise equiv. 11);
//3       =20- forces to use TS00 code (renormalized if data exists).
//3       =21- forces to use cross section fit if exists (otherwise equiv. 20).
//3       =22- forces to use a numerical table if exists (otherwise equiv. 21);
//3 The best values recommended are kopt = 12, 22 (12 is preferrable).
//3 uses file_no[1] and file_no[3] as indicators of the data array and fit params.
//3    output:
//3 info  =0 if no data exists, the Webber'93 or TS00 evaluation is used;
//3       =1 if data exist and the renormalized Webber's formulae were used;
//3       =2 if data exist and the renormalized TS00 formulae were used;
//3       =3 if data exist and a const cross section is fitted to the data;
//3       =10 Webber'93 code used in all cases (renormalized if data exists);
//3       =20 TS00 code used in all cases (renormalized if data exists);
//3       =-11,-21 if used the cross section fit;
//3    comments:
//3 The CS_DATA array containts the following data:
//3 Zi.Ai, Zf.Af, emev/n, CSmb, CSerr.
//3 if CSerr>0, it is the absolute error of the measurement;
//3 if CSerr<0, it is the relative error and the abs. err. will be calculated 
//3 automatically in the subroutine.
//3 Any character in the 1st column can be used to comment the line.
//3 isobaric cross sections (zf = 0):
//3 *-marks the errors corrected by imos
//3 #-marks where the data summed (e.g. B11+C11) are taken at close 
//3   (but not equal) energies. 
//3 kopt =0 alghorithm (not recommended): when no data exists for the required
//3 channel, the Webber's evaluation is used. When Webber's code give 0, the result 
//3 of TS00 code is used. When the data exist, it calculates the renormalization
//3 coefficient for the two approximations (Webber's and TS00) using the
//3 least-square method. Then Xi2 test is used to choose the best approximation.
//3 When both (Webber and TS98) approximations give 0, but the experimental
//3 value is non zero, uses a const cross section fitted to the data.
//3    references:
//3 [Ab94] Abdullin S.K. et al. 1994, Nucl. Phys. A 569, 753
//3 [Gl93] Glagolev V.V. et al. 1993, Z. Phys. C 60, 421
//3 [Ko99] Korejwo A. et al. 1999, Proc. 26th ICRC (Salt Lake City), OG 3.2.22 
//3 [LM69] Lebowitz E., Miller J.M. 1969, Phys. Rev. 177, 1548
//3 [Ni72] Nicholls J.E. et al. 1972, Nucl. Phys. A 181, 329
//3 [Ol83] Olson D.L., et al. 1983, Phys.Rev.C 28, 1602 (the same as [ 17.])
//3 [Ra79] Radin J.R., Gradsztajn E., Smith A.R. 1979, Phys. Rev. C 20, 787
//3 [RV84] Read S.M., Viola V.E., Jr. 1984, Atom. Data Nucl. Data Tables 31, 359
//3        if no err. is shown, relative err. 0.1 (cs>10mb), 0.2 (cs<10mb), 
//3        0.3 (cs<1mb) are assumed; sometimes I made it differently /imos.
//3 [TS98] Silberberg R., Tsao C.H., Barghouty A.F. 1998, ApJ 501, 911 
//3 [We90] Webber W.R. 1990, in AIP Conf. Proc. 203, ed. W.V.Jones et al.
//3          (NY: AIP), 294
//3 [We96] Webber W.R. 1996, private comm.(from Ramaty R. et al.1997,ApJ 488,730)
//3    !  relative err. 0.1 (cs>10mb), 0.2 (cs<10mb) are assumed.
//3 [We98] Webber W.R., et al. 1998, ApJ 508, 949
//3    !!RELATIVE ERRORS adopted: B =0.04; C =0.08; D =0.16; E =0.26
//3 [We98prc] Webber W.R. et al. 1998, PRC 58, 3539 (Tables 7,12-15)
//3    ! For Table 7 relative err. 0.1 (cs>10mb), 0.2 (cs<10mb) are assumed;
//3    !!REL. ERR. adopted: A =0.03; B =0.04; C =0.07; D =0.10; E =0.18; F =0.26
//3 -- All the following data are taken from the Transport Collab. on 21.04.99:
//3 [ 17.] EXCLUDED - Lindstorm P.J., et al. 1975, Report LBL-3650 (see [Ol83])
//3 [105.] Webber W.R., Kish J.C., Schrier D.A. 1990, Phys.Rev. C 41, 547
//3 [108.] Chen C.-X., et al. 1997, ApJ 479, 504
//3 [109.] Knott C.N., et al. 1997, Phys.Rev. C 56, 398
//3 [110.] Chen C.-X., et al. 1997, Phys.Rev. C 56, 1536
//3 [111.] ???
//3 [MM03] Moskalenko I.V., Mashnik S.G. 2003, Proc. 28th ICRC (Tsukuba), p.1969
//
//4 double eval_cs(double emev,int za1,int za2,int* info)
//4 the routine checks if there exists a channel with a cross section (mb) given
//4 by a table and if so interpolates it linearly.
//4    input:
//4 emev  - energy of the primary nucleus in MeV/nucleon;
//4 za1= (100*z1+a1) - primary charge and atomic number;
//4 za2= (100*z2+a2) - final charge and atomic number;
//4    output:
//4 info =-1 - no channel found in the table, returns 0;
//4      = 1 - channel found, but emev < the lowest grid pt, returns 0;
//4      = 2 - channel found, but emev > the highest grid pt, returns cs at max energy given;
//4      = 3 - emev falls exactly on the grid, returns cs;
//4      = 4 - returns interpolated value.
//=="====!===="====!===="====!===="====!===="====!===="====!===="====!===="====!

//**.****|****.****|****.****|****.****|****.****|****.****|****.****|****.****|
// GLOBAL definitions

#include <iostream>//AWS20031223
#include <cstdlib> //IMOS20020112 AWS20050624
#include <string>
#include <cctype>               //AWS20050624
#include <fstream>
#include <sstream>
#include <cmath>

#include "fort_interface.h"
#include "constants.h"

using namespace std;

#include <ErrorLogger.h>

#define max(a,b) (((a) > (b)) ? (a) : (b))
#define fnuc(z,a) (100 * (z) + (a))
#define inuc(b) (int)(100 * (b) + 0.1)

#define N_DATA_FILES 4                 // total number of data-files to read

char* data_filename[N_DATA_FILES] = {"isotope_cs.dat","nucdata.dat","p_cs_fits.dat","eval_iso_cs.dat"};
int n_data[3][N_DATA_FILES],         // their dimensions n1,n2,n3  
  file_no[N_DATA_FILES] = {0,1,2,3}; // 0=isotope_cs.dat, 1=nucdata.dat etc.
float* data_file[N_DATA_FILES];        // pointers to the data arrays

double eval_cs(double,int,int,int*);

//**.****|****.****|****.****|****.****|****.****|****.****|****.****|****.****|

void read_nucdata(const string& path) {

  float *tmp,*tmp1;
  int i,j,k,size;
  const int BufferSize=200;
  char readBuffer[BufferSize];
  ifstream data;
  
  for(j=0; j<N_DATA_FILES; j++)
    {

      const string fullFilename = path + "/" + data_filename[j];

      //data.open(data_filename[j]);                    // open file if exists
      data.open(fullFilename.c_str());
      if(data.fail())
	{

	  ostringstream errBuf;
	  errBuf << "Error opening file " << fullFilename;
	  ERROR(errBuf.str());
	  exit(1);
	  //cerr<<"read_nucdata: Error opening file "<<data_filename[j]<<endl;
	  //exit(1);
	}
      
      while(!isspace(data.get()) && !data.eof())      // skip comments:
	data.getline(readBuffer,BufferSize,'\n');    // any symbol in 1st col. 
      
      for(i=0; i<3; data >> n_data[i++][j]);          // read array's sizes
      data.getline(readBuffer,BufferSize,'\n');       // skip the rest of line
      
      for(size=1, i=0; i<3; size*=n_data[i++][j]);    // allocate space
      //      data_file[j] = (float*) calloc( (size_t) size, (size_t) sizeof(float));
      data_file[j] = new float[size];
      
      for(k = 0; k < size && !data.eof();)            // read data loop
	{
	  while(!isspace(data.get()) && !data.eof())   // skip comments:
            data.getline(readBuffer,BufferSize,'\n'); // any symbol in 1st col. 
	  for(i=0; i < n_data[0][j]; i++) data >> *(data_file[j]+k++);
	  data.getline(readBuffer,BufferSize,'\n');    // skip the rest of line
	}
      data.close();
    }
  return;
}

//**.****|****.****|****.****|****.****|****.****|****.****|****.****|****.****|
//IMOS20060420 cleanup arrays to prevent memory leaks
void cleanup_nucdata()
{
  for(int j=0; j<N_DATA_FILES; j++) 
    delete[] data_file[j];  //Gulli20070810  added [] to delete
  return;
}

//**.****|****.****|****.****|****.****|****.****|****.****|****.****|****.****|

double nucdata(int ksp,int iz,int ia,int K_electron,int izf,int iaf,int* izl,int* ial,double* To)
{
   int i,j,k,l,m,n,iy,iz0,ia0,iz4,ia4,iz5,ia5,iw[121],
       nksp=ksp*n_data[0][file_no[1]]*n_data[1][file_no[1]];
   float w[2][121], *decay = data_file[file_no[1]]+nksp;
   double b, xxx;

// STABLE & LONG-LIVED ISOTOPES (numbers in the table are the proton numbers)
// The long-lived isotopes listed in "longliv" table are included as stable;
   int stable[64][3] = {      // second index changes faster
       1,  0,  1,     1,  0,  1,     1,  0,  2,     2,  0,  2, //  A = 1- 4
       0,  0,  0,     3,  0,  3,     3,  0,  4,     0,  0,  0, //  A = 5- 8
       4,  0,  4,     4,  0,  5,     5,  0,  5,     6,  0,  6, //  A = 9-12
       6,  0,  6,     6,  0,  7,     7,  0,  7,     8,  0,  8, //  A =13-16
       8,  0,  8,     8,  0,  8,     9,  0,  9,    10,  0, 10, //  A =17-20
      10,  0, 10,    10,  0, 11,    11,  0, 11,    12,  0, 12, //  A =21-24
      12,  0, 12,    12,  0, 13,    13,  0, 13,    14,  0, 14, //  A =25-28
      14,  0, 14,    14,  0, 14,    15,  0, 15,    14,  0, 16, //  A =29-32
      16,  0, 16,    16,  0, 16,    17,  0, 17,    16, 17, 18, //  A =33-36
      17,  0, 18,    18,  0, 18,    18,  0, 19,    18, 19, 20, //  A =37-40
      19,  0, 20,    18,  0, 20,    20,  0, 20,    20,  0, 22, //  A =41-44
      21,  0, 21,    20,  0, 22,    22,  0, 22,    20,  0, 22, //  A =45-48
      22,  0, 23,    22, 23, 24,    23,  0, 24,    24,  0, 24, //  A =49-52
      24,  0, 25,    24, 25, 26,    25,  0, 26,    26,  0, 28, //  A =53-56
      26,  0, 27,    26,  0, 28,    27,  0, 28,    26, 27, 28, //  A =57-60
      28,  0, 28,    28,  0, 28,    28,  0, 29,    28,  0, 28  //  A =61-64
   };

// LONG-LIVED ISOTOPES (>~1y): Zi.Ai T_1/2(y) Zf.Af - [][][0] half-life shown for naked nucleus
   int nll = 25;                                 // - [][][1] half-life shown for H2-like atoms
   float longliv[25][2][3] = {   // third index changes faster
      1.03,  12.33,    2.03,     //  3H (b-) 3He   100% [ToI]
      1.03,  12.33,    2.03,     // no EC

      4.07,  0.,       4.07,     // stable
      4.07,  0.1459,   3.07,     //  7Be(EC) 7Li   100% [ToI]

      4.10,  1.60e6,   5.10,     // 10Be(b-)10B    100% [ToI]
      4.10,  1.60e6,   5.10,     // no EC

      6.14,  5.73e3,   7.14,     // 14C (b-)14N    100% [ToI]
      6.14,  5.73e3,   7.14,     // no EC

     11.22,  4.80e3,  10.22,     // 22Na(b+)22Ne        [M98]
     11.22,  2.60e0,  10.22,     // 22Na(EC?)22Ne       [ToI] T1/2(Lab)=2.60e0 y

     13.26,  9.10e5,  12.26,     // 26Al(b+)26Mg        [M98]
     13.26,  4.075e6, 12.26,     // 26Al(EC)26Mg        [M98] T1/2(Lab)=7.4e5 y [ToI]

     14.32,  172.,    16.32,     // 32Si(2b-)32S   100% [ToI] Si-P -S 
     14.32,  172.,    16.32,     // no EC

     17.36,  3.07e5,  18.36,     // 36Cl(b-)36Ar        [ToI]
     17.36,  1.58e7,  16.36,     // 36Cl(EC)36S         [ToI] T1/2(Lab)=3.01e5 y

     18.37,  0.,      18.37,     // stable
     18.37,  0.1,     17.37,     // 37Ar(EC)37Cl   100% [ToI] T1/2(Lab)=35.04 d

     18.39,  2.69e2,  19.39,     // 39Ar(b-)39K    100% [ToI]
     18.39,  2.69e2,  19.39,     // no EC

     19.40,  1.43e9,  20.40,     // 40K (b-)40Ca   89.3%[ToI] T1/2(Lab)=1.277e9 y incl 10.7% ECb+
     19.40,  1.43e9,  20.40,     // no EC

     20.41,  0.,      20.41,     // stable
     20.41,  1.03e5,  19.41,     // 41Ca(EC)41K    100% [ToI]

     18.42,  32.9,    20.42,     // 42Ar(2b-)42Ca  100% [ToI] Ar-K -Ca
     18.42,  32.9,    20.42,     // no EC

     22.44,  0.,      22.44,     // stable
     22.44,  49.,     20.44,     // 44Ti(ECb+)44Ca 100% [ToI] Ti(EC)Sc(b+)Ca

     23.49,  0.,      23.49,     // stable
     23.49,  0.903,   22.49,     // 49V (EC)49Ti   100% [ToI] 

     24.51,  0.,      24.51,     // stable
     24.51,  0.076,   23.51,     // 51Cr(EC)51V   <100% [ToI] 

     25.53,  0.,      25.53,     // stable
     25.53,  3.74e6,  24.53,     // 53Mn(EC)53Cr   100% [ToI]

     25.54,  6.30e5,  26.54,     // 54Mn(b-)54Fe        [W98]
     25.54,  0.855,   24.54,     // 54Mn(EC)54Cr        [ToI] T1/2(Lab)=312.3 d

     26.55,  0.,      26.55,     // stable
     26.55,  2.73e0,  25.55,     // 55Fe(EC)55Mn   100% [ToI]

     28.56,  4.00e4,  26.56,     // 56Ni(2b+)56Fe <100% [F99] Ni-Co-Fe
     28.56,  0.1,     26.56,     // 56Ni(ECb+)56Fe      [ToI] T1/2(Lab)=~30 d Ni(EC)Co(b+)Fe

     27.57,  0.,      27.57,     // stable
     27.57,  0.744,   26.57,     // 57Co(EC)57Fe   100% [ToI]

     28.59,  0.,      28.59,     // stable
     28.59,  7.60e4,  27.59,     // 59Ni(EC)59Co  <100% [ToI] [B76]

     27.60,  5.27e0,  28.60,     // 60Co(b-)60Ni   100% [ToI]
     27.60,  5.27e0,  28.60,     // no EC

     26.60,  1.50e6,  27.60,     // 60Fe(b-)60Co   100% [ToI]
     26.60,  1.50e6,  27.60,     // no EC

     28.63,  1.00e2,  29.63,     // 63Ni(b-)63Cu   100% [ToI]
     28.63,  1.00e2,  29.63      // no EC
   };
// K-capture nuclei - factor of 2 because only 1 electron
   for(i=0;i<nll;i++) if(longliv[i][0][1]!=longliv[i][1][1]) longliv[i][1][1] *=2.; 
   
// BOUNDARY NUCLEI 
// on the left side from the left boundary, the proton-emission is assumed;
// on the right side from the right boundary, the neutron-emission is assumed.
// ZB.AB; left boundary[][0]: Nn=0(1)28; right boundary[][1]: Np=1(1)28 
    int nb = 29;
    float boundary[29][2] = {  // second index changes faster
       1.01,  1.04,    3.04,  2.08,
       3.05,  3.11,    6.09,  4.14,
       6.10,  5.15,    8.13,  6.16,
       8.14,  7.21,   10.17,  8.22,
      12.20,  9.24,   12.21, 10.26,
      14.24, 11.30,   14.25, 12.30,
      15.27, 13.34,   16.29, 14.35,
      16.30, 15.38,   18.33, 16.40,
      20.36, 17.43,   20.37, 18.46,
      20.38, 19.48,   22.41, 20.51,
      22.42, 21.52,   24.45, 22.53,
      24.46, 23.55,   24.47, 24.59,
      25.49, 25.63,   26.51, 26.65,
      27.53, 27.69,   27.54, 28.69,
      28.56, 00.99
   };

   b = *To = *izl = *ial = 0;
   if(iz <= 0 || ia <= 0) return(0.);    // check against negative numbers,
   if(iz*ia > 1 && iz >= ia) return(0.); // non-existed nuclei,
   if(2864 < fnuc(iz,ia)) return(0.);    // Ni64 is the heaviest nucleus
   if(64 < ia) return(0.);               // A=64 is the maximal atomic number

// CHECK FOR NUCLEI OUTSIDE THE BOUNDARIES (p/n decay)
   iz0 = iz;
   ia0 = ia;
   if(ia>inuc(boundary[iz-1][1]-iz)) ia0=inuc(boundary[iz-1][1]-iz); // n -decay
   if(29>ia-iz) if(ia>inuc(modf(boundary[ia-iz][0], &xxx)))          // p -decay
   { 
      iz0=(int)boundary[ia-iz][0];
      ia0=inuc(boundary[ia-iz][0]-iz0);
   }

   for(i=0; i<121; iw[i++]=-1)  for(j=0; j<2; w[j++][i]=0.);

// SEARCH FOR A SPECIAL CASE (non beta decay)
   for(i=0; i<n_data[1][file_no[1]]; i++)
      if(fnuc(iz0,ia0) == inuc(*(decay +i*n_data[0][file_no[1]] +0)))
      {
          iw[0]  = i;            // if found, save the line number
          w[1][0]= 1.00;         // assign 1 to the branching ratio
          break;
      }

// STANDARD CASE (beta decay & long-lived isotopes)
   if(iw[0] < 0)
   {
      iz5 = iz0;
      ia5 = ia0;
// *** BETA DECAY ***
      if(iz0 > stable[ia0-1][2]) iz5 = stable[ia0-1][2];   // b+ decay
      if(iz0 < stable[ia0-1][0]) iz5 = stable[ia0-1][0];   // b- decay
// *** LONG-LIVED ISOTOPES (>~1 y) ***
      for(i=0; i<nll; i++)
         if(fnuc(iz5,ia5) == inuc(longliv[i][K_electron][0]))
         {
            *izl = iz5;
            *ial = ia5;
            *To = longliv[i][K_electron][1]*year2sec;
            if(!*To) *izl=*ial=0;
            iz5 = (int) longliv[i][K_electron][2];
            ia5 = inuc(longliv[i][K_electron][2]-iz5);
            break;
         }
      if(fnuc(izf,iaf)==fnuc(iz5,ia5) || fnuc(izf,iaf)==fnuc(*izl,*ial)) b = 1.;
      if(fnuc(iz0,ia0) == fnuc(*izl,*ial)) *izl = *ial = 0;
      return(b);
   }

// DEVELOPING A NETWORK OF DECAYS
   for(l=-1, m=0, ia4=1, i=0; i<4; ia4 =(int) pow(3.,++i))
   {
      for(l+=ia4, iy=0, n=0; n<ia4; n++, m++)
      {                                                      // check if there is
         if(iw[m] < 0) continue;                             // a required channel
         for(w[0][m]=0., k=2; k<8; k+=2)
         {
            w[0][l+3*n+k/2]                                  // store sec.nuclei
               =*(decay +iw[m]*n_data[0][file_no[1]] +k-1);
            w[1][l+3*n+k/2]                                  // store branchings
               =*(decay +iw[m]*n_data[0][file_no[1]] +k)*w[1][m];
            for(j=0; j<iw[m]; j++)                           // check if sec.nucleus
               if(*(decay +iw[m]*n_data[0][file_no[1]] +k-1) // also develops a
                  == *(decay +j*n_data[0][file_no[1]] +0))   // network of decays
               {
                  iw[l+3*n+k/2] = j;                         // store such a nucleus
                  iy = l+3*n+k/2;
               }  ///printf("%d %d %d %d %d %d\n",l,n,m,k,j,iy);
         }
      }
      if(iy == 0) break;
   }

// CHECK FOR STABILITY OF THE FINAL NUCLEI
   for(k=0; k<=l+3*n; k++)
   {
      *To = *izl = *ial = 0;
      if(w[0][k] == 0.) continue;
      iz4 = (int) w[0][k];
      ia4 = inuc(w[0][k]-iz4);
      iz5 = iz4;
      ia5 = ia4;
// *** BETA DECAY ***
      if(iz4 > stable[ia4-1][2]) iz5 = stable[ia4-1][2];   // b+ decay
      if(iz4 < stable[ia4-1][0]) iz5 = stable[ia4-1][0];   // b- decay
// *** LONG-LIVED ISOTOPES (>~1 y) ***
      for(i=0; i<nll; i++)
      {
         if(fnuc(iz5,ia5) != inuc(longliv[i][K_electron][0])) continue;
         *izl = iz5;
         *ial = ia5;
         *To = longliv[i][K_electron][1]*year2sec;
         if(!*To) *izl=*ial=0;
         iz5 = (int) longliv[i][K_electron][2];
         ia5 = inuc(longliv[i][K_electron][2]-iz5);
         break;
      }
      if(fnuc(izf,iaf) == fnuc(*izl,*ial) || fnuc(izf,iaf) == fnuc(iz5,ia5))
         return(w[1][k]);
   }
   return(b);
}

//**.****|****.****|****.****|****.****|****.****|****.****|****.****|****.****|

double isotope_cs(double emev,int iz,int ia,int izf,int iaf,int kopt,int* info)
{
   int a1,a2,i,j,size, itable=0, info1;
   float e1,y,err2,xi1,xi2, CSmb=0., f1=0., f2=0., T[11], a[3]={1.,1.,0.}, b[6];
   float *cs_data = data_file[file_no[0]], *p_cs = data_file[file_no[2]], *tp=T;
   double ej;
   
   e1 = emev;
   *info = kopt;

// CHECK if user wants to use specific program (the value of "kopt")
   if(kopt == 1) CSmb = wsigma_cc(iz,ia,izf,iaf,emev);       // Webber's code IMOS20020502
   if(kopt == 2) CSmb = yieldx_cc(iz,ia,izf,iaf,e1);         // TS code       IMOS20020502
   CSmb = max(0.,CSmb);
   if(kopt == 1 || kopt == 2) return(CSmb);

   a1 = fnuc(iz, ia);
   a2 = fnuc(izf,iaf);

// EVALUATED CROSS SECTIONS

   if(kopt == 12 || kopt == 22)
   {      
      CSmb = eval_cs(emev,a1,a2,&info1);
      if (info1 > 0) return(max(0.,CSmb));
      kopt--;                                          // if evaluation doesn't exist,  
   }                                                   // try other options

// if user wants, use THE CROSS SECTION FITS

   if(kopt == 11 || kopt == 21)
   {      
// special cases: Be, B - recursion calls
      if(izf != 0)
      {
// A = 10
         if(10 == iaf)  // B10 = B10 + C10 = a10 - Be10
         {
            b[0] = isotope_cs(emev,iz,ia,0,iaf,21,&j);
            if(j == -21)
            {                                          // B10
                if(510 == a2) 
                {
                   b[0]-=isotope_cs(emev,iz,ia,4,iaf,21,&j);
                   return(max(0.,b[0]));
                }
                if(5 < izf) return(0.);                // C10, =0
            }
         }
// A = 11
         if(11 == iaf)  // B11 = a11 = Be11 + B11 + C11
         {
            b[0] = isotope_cs(emev,iz,ia,0,iaf,21,&j);
            if(j == -21) 
            {
               if(511 == a2) return(max(0.,b[0]));     // B11
               return(0.);                             // =0 for the rest
            }
         }
      }
// straight search in the table
      for(i=0; i<n_data[1][file_no[2]]-1; i++, p_cs+=n_data[0][file_no[2]])
         if(a1 == inuc(*p_cs) && a2 == inuc(*(p_cs+1)))
         {
            for(p_cs+=2, j=0; j<6; b[j++]=*p_cs++);    // take the parameters
            if(b[0] >= 0.)                             // if positive use fit
            {
               *info=-kopt;
               if(emev < b[5]) return(0);              // fitting function
               b[0]*=(1.+sin(b[1]*pow(log10(emev),1.*b[2]))*exp(-b[3]*(emev-b[4])));
               return(max(0.,b[0]));
            }
            kopt = (int)(-b[0]+0.1);                   // negative b[0] gives kopt
         }
      if(izf == 0) return(0.);
   }

// CHECK if user wants to use specific program (the value of "kopt")
   if(kopt == 1) CSmb = wsigma_cc(iz,ia,izf,iaf,emev);       // Webber's code IMOS20020502
   if(kopt == 2) CSmb = yieldx_cc(iz,ia,izf,iaf,e1);         // TS code       IMOS20020502
   CSmb = max(0.,CSmb);
   if(kopt == 1 || kopt == 2) return(CSmb);

// STARTING THE ALGHORITHM

   for(i=0; i<11; T[i++] = 0.);

// CHECK the array: cs_data (is there a channel we are looking for ?)

   for(size=1, i=0; i<3; size*=n_data[i++][file_no[0]]);
   for(tp = T, i=0; i<size; i+=n_data[0][file_no[0]], tp = T, f1=0., f2=0.)
   {
      if(a1 != inuc(*(cs_data+i)))   continue;
      if(a2 != inuc(*(cs_data+i+1))) continue;

// if there is such a channel then the LEAST-SQUARE FIT

      itable++;
      if(*(cs_data+i+4) < 0.) *(cs_data+i+4) *= -*(cs_data+i+3);  // calc.abs.err.
      err2 = pow(*(cs_data+i+4),2);                               // err^2
      y = *(cs_data+i+3);                                         // cs measured
      ej = *(cs_data+i+2);                                        // @ energy

      if(kopt/10 != 2) f1=wsigma_cc(iz,ia,izf,iaf,ej);             // Webber IMOS20020502
      if(kopt/10 != 1) f2=yieldx_cc(iz,ia,izf,iaf,*(cs_data+i+2)); // TS     IMOS20020502

// calculations of the separate terms:
      *tp++ += f1*y /err2;       // Webber
      *tp++ += f1*f1/err2;
      *tp++ += f2*y /err2;       // TS
      *tp++ += f2*f2/err2;
      *tp++ += y    /err2;       // const cs
      *tp++ += 1.   /err2;
      
// calculation of terms for the Xi2 estimates
      *tp++ += y*y    /err2;
      *tp++ += 2.*f1*y/err2;
      *tp++ += f1*f1  /err2;
      *tp++ += 2.*f2*y/err2;
      *tp   += f2*f2  /err2;

// calculation of renormalization coefficients 
      for(j=0; j<3; j++) a[j]= (T[2*j+1] != 0.) ? T[2*j]/T[2*j+1]: a[j];
   }
   if(kopt == 3 && a[2] != 0.) return(a[2]);                  // const cr.sect.
   if(kopt/10 == 1) CSmb = wsigma_cc(iz,ia,izf,iaf,emev);     // Webber code IMOS20020502
   if(kopt/10 == 2) CSmb = yieldx_cc(iz,ia,izf,iaf,e1);       // TS code     IMOS20020502
   if(kopt/10 == 1 || kopt/10 == 2) return(max(0.,CSmb*a[kopt/10-1]));

// CHOOSE THE BEST APPROXIMATION (kopt = 0)
   if(itable < 2)                                         // no data or 1 pt.
   {  
      *info = itable;
      CSmb = a[0]*wsigma_cc(iz,ia,izf,iaf,emev);          // use Webber code     IMOS20020502
      if(CSmb <= 0.) 
      {                                                   // if W-code give 0,
         CSmb = a[1]*yieldx_cc(iz,ia,izf,iaf,e1);         // take the TS approx. IMOS20020502
         if(CSmb != 0. && itable == 1) *info = 2;
      }
   } else                                                 // data exists
   {
      xi1= T[6] -a[0]*T[7] +a[0]*a[0]*T[8];               // Xi2 evaluation 1
      xi2= T[6] -a[1]*T[9] +a[1]*a[1]*T[10];              // Xi2 evaluation 2
      if(xi1 < xi2)
      {
         *info = 1;
         CSmb = a[0]*wsigma_cc(iz,ia,izf,iaf,emev);       // renorm. Webber approx. IMOS20020502
      } else
      {
         *info = 2;
         CSmb = a[1]*yieldx_cc(iz,ia,izf,iaf,e1);         // renorm. TS approx.     IMOS20020502
      }
   }
   return(max(0.,CSmb));
}

//**.****|****.****|****.****|****.****|****.****|****.****|****.****|****.****|

double eval_cs(double emev,int za1,int za2,int* info)
{
   int i,size;
   float *eval = data_file[file_no[3]];
   double x[2]={-1.e10,1.e10},y[2]={0.,0.};

// CHECK the array: eval (is there a channel we are looking for ?)

   for(size=1, i=0; i<3; size*=n_data[i++][file_no[3]]);
   for(*info=0, i=0; i<size; i+=n_data[0][file_no[3]])
   {
      if(za1 != inuc(*(eval+i)))   continue;
      if(za2 != inuc(*(eval+i+1))) continue;

      if(x[0] < *(eval+i+2) && *(eval+i+2) <= emev)   // find lower energy pt
      { 
         x[0] = *(eval+i+2); 
         y[0] = *(eval+i+3);
      }
      if(emev <= *(eval+i+2) && *(eval+i+2) < x[1])   // find higher energy pt 
      { 
         x[1] = *(eval+i+2); 
         y[1] = *(eval+i+3);
      }
   }

   if(x[0]*x[1] < -1.e19) { *info = -1; return(0.); } // no evaluation found, return 0

   if(x[0] <   0.) { *info = 1; return(0.); }         // no lower grid pt, return 0
   if(x[1] > 9.e9) { *info = 2; return(y[0]); }       // no higher grid pt, extrapolate

   if(x[1]-x[0] == 0.) { *info = 3; return(y[1]); }   // emev falls exactly on the grid

   for(*info = 4, i=0; i<2; i++) x[i] = log10(x[i]);
   return(y[0]+(log10(emev)-x[0])*(y[1]-y[0])/(x[1]-x[0]));// interpolate
}

//**.****|****.****|****.****|****.****|****.****|****.****|****.****|****.****|

// A SAMPLE TEST PROGRAM
/*
#include<stdio.h>
#define min(a,b) (((a) > (b)) ? (b) : (a))
extern "C" void   set_sigma__(int*);                           // Webber's init IMOS20020502
extern "C" void   yieldx_(int*,int*,int*,int*,float*,float*);  // TS code       IMOS20020502
extern "C" double wsigma_(int*,int*,int*,int*,double*);        // Webber's code IMOS20020502
main()
{
   FILE *f1,*f2;
   int i,j,k,n,m,info,iz,ia,jz,ja,sz[9],sa[9],izl,ial, ksp=0, cdr=51, 
      kopt[9]={0,1,2,10,11,12,20,21,22};
   float *tmp;
   double ej,b,To,cs;
   char filename[20], ch;                  // IMOS20020502

   cout<<"\nTESTING THE NUCLEAR PACKAGE / version of 3/23/2000, imos \n\n";

// test read_nucdata
   
   cout<<"read_nucdata: I am starting\n";
   read_nucdata();                                // read data from files
   cout<<"read_nucdata: I am finishing\n";

   f1=fopen("test_cs.1","w");                     // print data to a file
   for(fprintf(f1,"\n"), n=0; n<N_DATA_FILES; n++, fprintf(f1,"\n")) 
      for(tmp = data_file[n], k=0; k<n_data[2][n]; k++, fprintf(f1,"\n")) 
         for(i=0; i<n_data[1][n]; i++, fprintf(f1,"\n")) 
            for(j=0; j<n_data[0][n]; j++) fprintf(f1,"%11.3f",*tmp++);
   fclose(f1);
   cout<<"read_nucdata: test finished, check file test_cs.1\n\n";
         
// test nucdata

   cout<<"test nuclear reaction network ?(y/n)\n";
   cin>>ch;
   if(ch != 'y' && ch != 'Y') goto test3;

   cout<<"nucdata: start testing\n";
   f1=fopen("test_cs.2","w");
   f2=fopen("test_cs.2K","w");
   for(iz=1; iz<29; iz++)                               // check all the reactions
   {
      for(i=0;i<20; fprintf(f1,"-"), i++); 
      for(fprintf(f1,"%4d  ",iz), i=0;i<20; fprintf(f1,"-"), i++); fprintf(f1,"\n");
      for(i=0;i<20; fprintf(f2,"-"), i++); 
      for(fprintf(f2,"%4d  ",iz), i=0;i<20; fprintf(f2,"-"), i++); fprintf(f2,"\n");
      for(ia=2*iz-2; ia<2.5*iz+4.2; ia++)
         for(jz=1; jz<29; jz++)
            for(ja=2*jz-4; ja<2.5*jz+4.2; ja++)
	    {
               b = nucdata(ksp,iz,ia,0,jz,ja,&izl,&ial,&To);
               if(b == 0.) continue;
               fprintf(f1,"%5d.%2d%5d.%2d%5d.%2d%14.3E%8.3f\n",iz,ia,jz,ja,izl,ial,To,b);
               printf("%5d.%2d%5d.%2d%5d.%2d%14.3E%8.3f\n",iz,ia,jz,ja,izl,ial,To,b);

               b = nucdata(ksp,iz,ia,1,jz,ja,&izl,&ial,&To);
               if(b == 0.) continue;
               fprintf(f2,"%5d.%2d%5d.%2d%5d.%2d%14.3E%8.3f\n",iz,ia,jz,ja,izl,ial,To,b);
               printf("%5d.%2d%5d.%2d%5d.%2d%14.3E%8.3f\n",iz,ia,jz,ja,izl,ial,To,b);
            }
   }
   cout<<"nucdata: test finished, check file test_cs.2\n\n";

// test isotope_cs

test3:
   cout<<"test cross sections ?(y/n)\n";
   cin>>ch;
   if(ch != 'y' && ch != 'Y') goto test4;
   
   cout<<"cross sections: start testing \n";
   set_sigma__(&cdr);           // initialization of Webber's code IMOS20020502
   
   f1=fopen("test_cs.3","w");
   for(ej=900., iz=1; iz<29; iz++)
      for(ia=2*iz-2; ia<2.5*iz+4.2; ia++)
      {
         if(!nucdata(0,iz,ia,0,iz,ia,&izl,&ial,&To)) continue;
         for(jz=2; jz<iz+1; jz++)
            for(ja=max(1,2*jz-2); ja<min(2*jz+4,ia); ja++, fprintf(f1,"\n"))
            {
               fprintf(f1,"%7.2f%7.2f%12.3E",iz+ia/100.,jz+ja/100.,ej); 
               printf("%7.2f%7.2f%12.3e",iz+ia/100.,jz+ja/100.,ej); 
               for(i=0; i<3; i++) 
                  printf("%12.3e",isotope_cs(ej,iz,ia,jz,ja,kopt[i],&j)); 
               cout<<endl;
               for(i=0; i<9; i++) 
                  fprintf(f1,"%12.3E",isotope_cs(ej,iz,ia,jz,ja,kopt[i],&j));
            }
      }
   fclose(f1);
   cout<<"cross sections: test finished, check file test_cs.3\n\n";

// testing cross section fits

test4:
   cout<<"test cross section fits ?(y/n)\n";
   cin>>ch;
   if(ch != 'y' && ch != 'Y') exit(0);
   set_sigma__(&cdr);          // initialization of Webber's code IMOS20020502
test4_1:   
   cout<<"cross section fits: start testing \n\n"
       <<"give IZ IA  IZF IAF IZF IAF... 0 0:\n";
   cin>>iz>>ia>>jz>>ja;
   for(n=0;n<3;n++) 
   {
      cin>>sz[n]>>sa[n];
      if(!(sz[n]*sa[n])) break;
   }
   filename[0]='\0';
   sprintf(filename,"x%1d%1d%1d%1d_%1d%1d%1d%1d",iz/10,iz%10,ia/10,ia%10,jz/10,jz%10,ja/10,ja%10);
   
   f1=fopen(filename,"w");
   for(ej=1, k=0; k<100; k++, ej*=1.1, printf("%12.3e\n",ej))
   {
      fprintf(f1,"\n%12.3e",ej);
      for(i=0; i<9; i++) 
      { 
         cs=isotope_cs(ej,iz,ia,jz,ja,kopt[i],&j);
	 for(m=0; m<n; m++) cs+=isotope_cs(ej,iz,ia,sz[m],sa[m],kopt[i],&j);
         fprintf(f1,"%12.3e",cs);
      }
   }
   fprintf(f1,"\n");
   fclose(f1);
   cout<<"output file: "<<filename<<endl;   
   cout<<"do you want one more channel ?(y/n)\n";
   cin>>ch;
   if(ch == 'y' || ch == 'Y') goto test4_1;
}
//**.****|****.****|****.****|****.****|****.****|****.****|****.****|****.****|

// Webber's isotopic production cross section  IMOS20020502
double wsigma_cc(int IZ, int IA, int IZF, int IAF, double E)
{
   return( wsigma_(&IZ,&IA,&IZF,&IAF,&E) );
}

//**.****|****.****|****.****|****.****|****.****|****.****|****.****|****.****|

// Silberberg & Tsao isotopic production cross section  IMOS20020502
double yieldx_cc(int IZ, int IA, int IZF, int IAF, float E)
{
   float CSmb;
   yieldx_(&IZ,&IA,&IZF,&IAF,&E,&CSmb);
   return( 1.*CSmb );
}
*/


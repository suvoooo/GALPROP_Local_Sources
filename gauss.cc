// a simple Gaussian random number generator using accept/reject method
// samples within 3 sigma of mean

using namespace std;//AWS20050624
#include<iostream>  //AWS20050624
#include<cmath>     //AWS20050624
#include<cstdlib> // IMOS20010816 AWS20050624
static int init=0;

double gauss(double mean, double sigma){
double x,y;
double xmin,xmax;
double v;
unsigned seed;


 if(sigma==0.0){x=mean;return x;} // otherwise loops forever for this case

xmin=mean-sigma*3.0;
xmax=mean+sigma*3.0;

if(init==0){
init =1;
seed=1234;
srand(seed);
}

int accept=0;

while(accept==0){
x=xmin+double(rand())/RAND_MAX*(xmax-xmin);
y=double(rand())/RAND_MAX;
v=exp(-(x-mean)*(x-mean)/(2.0*sigma*sigma));
if(y < v) accept=1;

//cout<<x<<" "<<y<<" "<<v<<" "<<accept<<endl;
}
return x;
}
/////////////////////////////////////////////
int test_gauss(){

double mean, sigma;
mean=20.0;
sigma=5.0;
double g;

int nrept=20;
int n=1000;

cout<<"test_gauss:"<<endl;
cout<<"number of samples="<<n<<endl;
cout<<"number of runs   ="<<nrept<<endl;

cout<<"  true mean="<< mean<<"        true sigma="<< sigma<<endl;

for (int irept=0;irept<nrept;irept++){

float sum=0.;
float sumsq=0.;

for (int i=0;i<n;i++){
g=gauss(mean,sigma);
sum+=g;
sumsq+=g*g;
//cout<<g<<endl;
}
float gmean=sum/n;
float gsigma=sqrt((sumsq-2*gmean*sum + gmean*gmean*n)/n);

cout<<"sample mean="<<gmean<<     " sample sigma="<<gsigma<<endl;
}
cout<<endl;

return 0;}

//int main(){test_gauss();return 0;} // use for standalone test

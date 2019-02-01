#include <string.h>
#include <complex.h>
#include <fstream>

typedef double vec3[3];
using namespace std;

#ifndef calc2DIR_H
#define calc2DIR_H
#define IR2DOK 0
#define PI M_PI
#define HBAR 5.308837367 // in cm-1*ps
#define PSTR        "||||||||||||||||||||||||||||||||||||||||||||||||||"
#define PWID        50
class IR2D
{
    public:
        // default constructor and destructor
        IR2D( string _inpf_ );
        ~IR2D();

        // variables -- set default parameters
        string _efile_="e.dat";         // name for energy file
        string _dfile_="d.dat";         // name for dipole file
        string _ofile_="spec";          // name for output file
        double dt = 0.010 ;             // time step in ps
        int    nchrom=1;                // number of chromophores
        int    n1ex;                    // number of one-exciton states
        int    n2ex;                    // number of one-exciton states
        double t1t3_max=2.5;            // t1 and t3 max in ps
        double t2=0.;                   // t2 in ps
        double lifetime_T1=0.2;         // lifetime T1 in ps
        double lifetime_T2;             // lifetime T2 in ps
        int    nsamples=1 ;             // number of samples
        double sample_every=10.;        // how often to take a new sample in ps
        int    trjlen=0;                // number of frames in trajectory files
        int    fftlen=4096;             // number of data points for 1 dimension of fft
        int    t1t3_npoints;            // number of data points for t1 and t3 dimensions
        double window0 = 1400;          // lower limit of spectral window in cm-1
        double window1 = 1700;          // upper limit if spectral window in cm-1
        double shift;                   // reference frequency in cm-1

        // arrays to hold Hamiltonian and dipole
        double *H1;                     // one-exciton hamiltonian at current t
        complex<double> *eiH1_t0t1;     // one-exciton hamiltonian integrated from time t0 to t1
        complex<double> *eiH1_t1t2;     // one-exciton hamiltonian integrated from time t1 to t2
        complex<double> *eiH1_t2t3;     // one-exciton hamiltonian integrated from time t2 to t3
        complex<double> *eiH1_t1_last;  // one-exciton hamiltonian at previous timestep
        complex<double> *eiH1_t2_last;  // one-exciton hamiltonian at previous timestep
        complex<double> *eiH1_t3_last;  // one-exciton hamiltonian at previous timestep
        double *mu_eg_x;                // ground-to-one-exciton dipole vector
        double *mu_eg_y;                // ground-to-one-exciton dipole vector
        double *mu_eg_z;                // ground-to-one-exciton dipole vector
        double *mu_eg_t0_x;             // ground-to-one-exciton dipole vector at t0 - x component
        double *mu_eg_t1_x;             // ground-to-one-exciton dipole vector at t1 - x component
        double *mu_eg_t2_x;             // ground-to-one-exciton dipole vector at t2 - x component
        double *mu_eg_t3_x;             // ground-to-one-exciton dipole vector at t3 - x component
        double *mu_eg_t0_y;             // ground-to-one-exciton dipole vector at t0 - y component
        double *mu_eg_t1_y;             // ground-to-one-exciton dipole vector at t1 - y component
        double *mu_eg_t2_y;             // ground-to-one-exciton dipole vector at t2 - y component
        double *mu_eg_t3_y;             // ground-to-one-exciton dipole vector at t3 - y component
        double *mu_eg_t0_z;             // ground-to-one-exciton dipole vector at t0 - z component
        double *mu_eg_t1_z;             // ground-to-one-exciton dipole vector at t1 - z component
        double *mu_eg_t2_z;             // ground-to-one-exciton dipole vector at t2 - z component
        double *mu_eg_t3_z;             // ground-to-one-exciton dipole vector at t3 - z component

        // complex constants
        const complex<double> img          = {0.,1.};   
        const complex<double> complex_one  = {1.,0.};   
        const complex<double> complex_zero = {0.,0.};
      
        // response functions
        complex<double> *R1D;           // Linear response function
        complex<double> *R2D_R1;        // third order response function rephasing
        complex<double> *R2D_R2;        // third order response function non-rephasing

        // file handles
        ifstream efile;
        ifstream dfile;

        // member functions
        void fileOpenErr( string _fn_ );
        void fileReadErr( string _fn_ );
        int  readParam( string _inpf_ );
        int  readEframe( int frame );
        int  readDframe( int frame );
        int  setMUatT( string which );
        int  propigateH1( int t0, int t1, string which );
        int  doeiH( complex<double> *eiH1, double *H1, int n1ex );
        int  writeR1D();
        int  writeR2D();
        int  write1Dfft();
        int  write2DRabs();
        int  write2Dout( complex<double> *data, string fn, string which, int n);
        double dot3( vec3 a, vec3 b );
        template<class T> void tellParam( string param, T value );
        complex<double> getR1D( );
        complex<double> getR2D( int t1, int t2, string which );
        
};
#endif

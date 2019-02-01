#include <string.h>
#include <fstream>
#include <sstream>
#include <iostream>
#include <stdio.h>
#include <complex.h>
#include <fftw3.h>
#include <math.h>
#include <algorithm>
#include "calc2DIR.h"
#define MKL_Complex16 complex<double>
#include <mkl.h>

using namespace std;

IR2D::IR2D( string _inpf_ )
// Default Constructor
{
    string line;

    // read input file
    if ( readParam( _inpf_ ) != IR2DOK ) exit( EXIT_FAILURE );

    // allocate variable arrays
    R1D             = new complex<double>[t1t3_npoints]();
    R2D_R1          = new complex<double>[t1t3_npoints*t1t3_npoints]();
    R2D_R2          = new complex<double>[t1t3_npoints*t1t3_npoints]();
    H1              = new double [n1ex*n1ex]();
    eiH1_t0t1       = new complex<double>[n1ex*n1ex]();
    eiH1_t1t2       = new complex<double>[n1ex*n1ex]();
    eiH1_t2t3       = new complex<double>[n1ex*n1ex]();
    eiH1_t1_last    = new complex<double>[n1ex*n1ex]();
    eiH1_t2_last    = new complex<double>[n1ex*n1ex]();
    eiH1_t3_last    = new complex<double>[n1ex*n1ex]();
    mu_eg_x         = new double[n1ex]();
    mu_eg_y         = new double[n1ex]();
    mu_eg_z         = new double[n1ex]();
    mu_eg_t0_x      = new double[n1ex]();
    mu_eg_t1_x      = new double[n1ex]();
    mu_eg_t2_x      = new double[n1ex]();
    mu_eg_t3_x      = new double[n1ex]();
    mu_eg_t0_y      = new double[n1ex]();
    mu_eg_t1_y      = new double[n1ex]();
    mu_eg_t2_y      = new double[n1ex]();
    mu_eg_t3_y      = new double[n1ex]();
    mu_eg_t0_z      = new double[n1ex]();
    mu_eg_t1_z      = new double[n1ex]();
    mu_eg_t2_z      = new double[n1ex]();
    mu_eg_t3_z      = new double[n1ex]();

    // open the energy and dipole files
    efile.open( _efile_, ios::binary );    
    if ( ! efile.is_open() ) { fileOpenErr( _efile_ ); exit(EXIT_FAILURE);}
    dfile.open( _dfile_, ios::binary );    
    if ( ! dfile.is_open() ) { fileOpenErr( _dfile_ ); exit(EXIT_FAILURE);}
}

IR2D::~IR2D()
// Default Destructor
{
    // destroy arrays
    delete [] R1D;
    delete [] R2D_R1;
    delete [] R2D_R2;
    delete [] H1;
    delete [] eiH1_t0t1;
    delete [] eiH1_t1t2;
    delete [] eiH1_t2t3;
    delete [] eiH1_t1_last;
    delete [] eiH1_t2_last;
    delete [] eiH1_t3_last;
    delete [] mu_eg_x;
    delete [] mu_eg_y;
    delete [] mu_eg_z;
    delete [] mu_eg_t0_x;
    delete [] mu_eg_t1_x;
    delete [] mu_eg_t2_x;
    delete [] mu_eg_t3_x;
    delete [] mu_eg_t0_y;
    delete [] mu_eg_t1_y;
    delete [] mu_eg_t2_y;
    delete [] mu_eg_t3_y;
    delete [] mu_eg_t0_z;
    delete [] mu_eg_t1_z;
    delete [] mu_eg_t2_z;
    delete [] mu_eg_t3_z;

    // close files
    efile.close();
    dfile.close();
}

int IR2D::readParam( string _inpf_ )
// read the input file to get parameters
{
    string line, para, value;

    // open file for reading
    ifstream inpf( _inpf_.c_str() );
    if ( ! inpf.is_open() ) { fileOpenErr( _inpf_ ); return 1; }
    cout << ">>> Reading simulation parameters from " << _inpf_ << endl;

    // parse input file
    while ( getline( inpf, line ) )
    {
        // create string stream and parse parameter name and value
        istringstream line2( line );
        line2 >> para >> value;

        // skip comments
        if ( para[0] == '#' ) continue;

        // transform parameter name to lower case
        transform( para.begin(), para.end(), para.begin(), ::tolower);

        // save parameters as class variable
        if      ( para.compare("energy_file") == 0 )  _efile_      = value;
        else if ( para.compare("dipole_file") == 0 )  _dfile_      = value;
        else if ( para.compare("trjlen")      == 0 )  trjlen       = stoi(value);
        else if ( para.compare("output_file") == 0 )  _ofile_      = value;
        else if ( para.compare("time_step")   == 0 )  dt           = stof(value);
        else if ( para.compare("nchrom")      == 0 )  nchrom       = stoi(value);
        else if ( para.compare("t1t3_max")    == 0 )  t1t3_max     = stof(value);
        else if ( para.compare("t2")          == 0 )  t2           = stof(value);
        else if ( para.compare("lifetimet1")  == 0 )  lifetime_T1  = stof(value);
        else if ( para.compare("nsamples")    == 0 )  nsamples     = stoi(value);
        else if ( para.compare("sample_every") == 0 ) sample_every = stof(value);
        else if ( para.compare("fftlen")      == 0 )  fftlen       = stoi(value);
        else if ( para.compare("window0")     == 0 )  window0      = stof(value);
        else if ( para.compare("window1")     == 0 )  window1      = stof(value);
        else cerr << "\tWARNING:: parameter " << para << " is not recognized." << endl;
    }
    inpf.close();

    // some output to confirm parameters
    tellParam<string>( "energy_file", _efile_ );
    tellParam<string>( "dipole_file", _dfile_ );
    tellParam<int>( "trjlen", trjlen );
    tellParam<string>( "output_file", _ofile_ );
    tellParam<double>( "time_step", dt );
    tellParam<int>( "nchrom", nchrom );
    tellParam<double>( "t1t3_max", t1t3_max );
    tellParam<double>( "t2", t2 );
    tellParam<double>( "lifetimeT1", lifetime_T1 );
    tellParam<int>( "nsamples", nsamples );
    tellParam<int>( "sample_every", sample_every );
    tellParam<int>( "fftlen", fftlen );
    tellParam<double>( "window0", window0 );
    tellParam<double>( "window1", window1 );
    cout << ">>> Done reading simulation parameters from " << _inpf_ << endl;

    if ( trjlen < static_cast<int>(((nsamples-1)*sample_every + (2*t1t3_max + t2))/dt) ){
        cout << "WARNING:: The given trajectory length is not long enough.\n" << 
                      "\t  Must be " <<  ((nsamples-1)*sample_every + (2*t1t3_max + t2))/dt << 
                      " frames long.\n\t  Check input file. Aborting." << endl;
        exit(EXIT_FAILURE);
    }

    // set number of points in response functions based on t1t3_max and dt
    t1t3_npoints = static_cast<int>(t1t3_max/dt);

    // set number of one- and two-exciton states
    n1ex = nchrom;
    n2ex = nchrom*(nchrom+1)/2;

    // shift energies by mean of spectral window to avoid high frequency oscillations
    shift = 0.5*(window1 + window0);
    cout << ">>> Shifting frequencies by: " << shift << " cm." << endl;

    // set T2 lifetime as 2T1 (see Hamm and Zanni p29 for discussion)
    // also see Liang and Jansen JCTC 2012 eq 14 and 16
    lifetime_T2 = 2*lifetime_T1;
    cout << ">>> Setting T2 to 2 T1: " << lifetime_T2 << " ps." << endl;

    return IR2DOK;
}

template<class T>
void IR2D::tellParam( string param, T value )
// write parameter and value to output
{
    cout << "\tSetting " << param << " to " << value << "." << endl; 
}

void IR2D::fileOpenErr( string _fn_ )
// give an error message when I cant open a file
{
    cerr << "ERROR:: Could not open " << _fn_ << "." << endl;
}

void IR2D::fileReadErr( string _fn_ )
// give an error message when I cant open a file
{
    cerr << "ERROR:: Reading " << _fn_ << " failed. Probably reached EOF...Aborting." << endl;
}

int IR2D::readEframe( int frame )
// Read the energy file
{
    int     frameTmp, nelm, i, j;
    float   *Htmp;
    int64_t file_offset;

    nelm = n1ex*(n1ex+1)/2;
    Htmp = new float[nelm]();
    file_offset = frame*(sizeof(int)+sizeof(float)*(nelm));
    efile.seekg( file_offset );

    // the energy for each frame is stored on a seperate line
    // with an integer frame number at the beginning of the line
    efile.read( (char*)&frameTmp, sizeof(int) );
    if ( not efile.good() ){ fileReadErr( _efile_ ); return 1;}
    efile.read( (char*)Htmp, nelm*sizeof(float) );
    if ( not efile.good() ){ fileReadErr( _efile_ ); return 1;}

    // save the one-exciton Hamiltonian
    for( i = 0; i < n1ex; i ++ ){
    for( j = i; j < n1ex; j++ ){
        H1[i*n1ex + j] = Htmp[i*n1ex + (j-i)];
        H1[j*n1ex + i] = H1[i*n1ex + j];
    }
        H1[i*n1ex + i] -= shift; // subtract the frequency shift
    }

    delete [] Htmp;
    return IR2DOK;
}

int IR2D::readDframe( int frame )
// Read the dipole file
{
    int     frameTmp, i, ex;
    float   *dipoleTmp;
    int64_t file_offset;

    dipoleTmp = new float[ 3*n1ex ];
    file_offset = frame*(sizeof(int)+sizeof(float)*3*n1ex);
    dfile.seekg( file_offset );

    dfile.read( (char*)&frameTmp, sizeof(int) );
    if ( not dfile.good() ){ fileReadErr( _dfile_ ); return 1;}
    dfile.read( (char*)dipoleTmp, sizeof(float)*3*n1ex );
    if ( not dfile.good() ){ fileReadErr( _dfile_ ); return 1;}

    // put these into the dipole vector variable
    // note that in the bin file all x's come first, then y's, etc
    for ( ex = 0; ex < n1ex; ex ++ ){
        mu_eg_x[ex] = dipoleTmp[ex];
        mu_eg_y[ex] = dipoleTmp[ex+n1ex];
        mu_eg_z[ex] = dipoleTmp[ex+2*n1ex];
    }

    delete[] dipoleTmp;
    return IR2DOK;
}

int IR2D::setMUatT( string which )
// set dipole vector for a given time
{
    int ex;

    if ( which.compare("t0") == 0 ){
        memcpy( mu_eg_t0_x, mu_eg_x, sizeof(double)*n1ex );
        memcpy( mu_eg_t0_y, mu_eg_y, sizeof(double)*n1ex );
        memcpy( mu_eg_t0_z, mu_eg_z, sizeof(double)*n1ex );
    }
    else if ( which.compare("t1") == 0 ){
        memcpy( mu_eg_t1_x, mu_eg_x, sizeof(double)*n1ex );
        memcpy( mu_eg_t1_y, mu_eg_y, sizeof(double)*n1ex );
        memcpy( mu_eg_t1_z, mu_eg_z, sizeof(double)*n1ex );
    }
    else if ( which.compare("t2") == 0 ){
        memcpy( mu_eg_t2_x, mu_eg_x, sizeof(double)*n1ex );
        memcpy( mu_eg_t2_y, mu_eg_y, sizeof(double)*n1ex );
        memcpy( mu_eg_t2_z, mu_eg_z, sizeof(double)*n1ex );
    }
    else if ( which.compare("t3") == 0 ){
        memcpy( mu_eg_t3_x, mu_eg_x, sizeof(double)*n1ex );
        memcpy( mu_eg_t3_y, mu_eg_y, sizeof(double)*n1ex );
        memcpy( mu_eg_t3_z, mu_eg_z, sizeof(double)*n1ex );
    }

    else{
        cout << "ERROR:: IR2D::setMUatT which= " << which << " is unknown." << endl;
        return 1;
    }

    return IR2DOK;
}

int IR2D::propigateH1( int t0, int t1, string which )
// integrate the 1-exciton Hamiltonian
{
    int i;
    complex<double> *eiH1, *tmp;
    
    eiH1 = new complex<double>[n1ex*n1ex];
    tmp  = new complex<double>[n1ex*n1ex];
    
    // deterimine e^-iH1dt/(2*Hbar), return to eiH1
    doeiH( eiH1, H1, n1ex );

    // do the integration
    if ( which.compare("t0-t1") == 0 ){
        // do the integration
        if ( t1 == t0 ){
            // at t1==t0, the integral is 0, so e^0 is Identity matrix
            for (i = 0; i < n1ex; i ++ ){
                eiH1_t0t1[i*n1ex +i] = complex_one;
            }
        }
        else{
            // integrate with the trapezoid rule using the last and current
            // exponentiated hamiltonians
            cblas_zgemm( CblasRowMajor, CblasNoTrans, CblasNoTrans, \
                         n1ex, n1ex, n1ex, &complex_one, eiH1_t0t1, n1ex, \
                         eiH1_t1_last, n1ex, &complex_zero, tmp, n1ex );
            cblas_zgemm( CblasRowMajor, CblasNoTrans, CblasNoTrans, \
                         n1ex, n1ex, n1ex, &complex_one, tmp, n1ex, \
                         eiH1, n1ex, &complex_zero, eiH1_t0t1, n1ex );
        }
        memcpy( eiH1_t1_last, eiH1, sizeof(complex<double>)*n1ex*n1ex );
    }
    else if ( which.compare("t1-t2") == 0 ){
        // do the integration
        if ( t1 == t0 ){
            // at t1==t0, the integral is 0, so e^0 is Identity matrix
            for (i = 0; i < n1ex; i ++ ){
                eiH1_t1t2[i*n1ex +i] = complex_one;
            }
        }
        else{
            // integrate with the trapezoid rule using the last and current
            // exponentiated hamiltonians
            cblas_zgemm( CblasRowMajor, CblasNoTrans, CblasNoTrans, \
                         n1ex, n1ex, n1ex, &complex_one, eiH1_t1t2, n1ex, \
                         eiH1_t2_last, n1ex, &complex_zero, tmp, n1ex );
            cblas_zgemm( CblasRowMajor, CblasNoTrans, CblasNoTrans, \
                         n1ex, n1ex, n1ex, &complex_one, tmp, n1ex, \
                         eiH1, n1ex, &complex_zero, eiH1_t1t2, n1ex );
        }
        memcpy( eiH1_t2_last, eiH1, sizeof(complex<double>)*n1ex*n1ex );
    }
    else if ( which.compare("t2-t3") == 0 ){
        // do the integration
        if ( t1 == t0 ){
            // at t1==t0, the integral is 0, so e^0 is Identity matrix
            for (i = 0; i < n1ex; i ++ ){
                eiH1_t2t3[i*n1ex +i] = complex_one;
            }
        }
        else{
            // integrate with the trapezoid rule using the last and current
            // exponentiated hamiltonians
            cblas_zgemm( CblasRowMajor, CblasNoTrans, CblasNoTrans, \
                         n1ex, n1ex, n1ex, &complex_one, eiH1_t2t3, n1ex, \
                         eiH1_t3_last, n1ex, &complex_zero, tmp, n1ex );
            cblas_zgemm( CblasRowMajor, CblasNoTrans, CblasNoTrans, \
                         n1ex, n1ex, n1ex, &complex_one, tmp, n1ex, \
                         eiH1, n1ex, &complex_zero, eiH1_t2t3, n1ex );
        }
        memcpy( eiH1_t3_last, eiH1, sizeof(complex<double>)*n1ex*n1ex );
    }
    else{
        cout << "ERROR:: IR2D::propigateH1 which= " << which << \
                " is unknown." << endl; return 1;
    }

    delete [] eiH1;
    delete [] tmp;
    return IR2DOK;
}

int IR2D::doeiH( complex<double> *eiH, double *H, int N )
// get the exponential of exp(i H dt/hbar)
{
    complex<double> *evec, *tmp, arg;
    double *W;
    int i = 0, j = 0, info;

    // allocate arrays
    evec = new complex<double>[N*N]();
    tmp  = new complex<double>[N*N]();
    W    = new double[N]();

    // get eigenvalues and eigenvector of H
    if ( info = LAPACKE_dsyevd( LAPACK_ROW_MAJOR, 'V', 'U', N, H, N, W ) != 0 ){
        cout << "WARNING:: LAPACKE_dsyeved returned " << info << endl;
        return 1;
    }

    // exponentiate eigenvalues
    for ( i = 0; i < N; i ++ ){
        arg = -img*W[i]*dt/(2.*HBAR); // divide by two to do trapezoid integration
        eiH[i*N+i] = exp(arg);
    }

    // promote eigenvectors to complex
    for ( i = 0; i < N; i ++ ){
    for ( j = 0; j < N; j ++ ){
        evec[i*N+j].real(H[i*N+j]);
        evec[i*N+j].imag(0.);
    }
    }
    
    // convert back to original basis using matrix multiplications
    cblas_zgemm( CblasRowMajor, CblasNoTrans, CblasNoTrans, \
                 N, N, N, &complex_one, evec, N, eiH, N, \
                 &complex_zero, tmp, N );
    cblas_zgemm( CblasRowMajor, CblasNoTrans, CblasTrans, \
                 N, N, N, &complex_one, tmp, N, evec, N, \
                 &complex_zero, eiH, N );

    // delete arrays
    delete [] evec;
    delete [] tmp;
    delete [] W;

    return IR2DOK;
}

complex<double> IR2D::getR1D()
// return the linear response function at a given t1
{
    int i;
    complex<double> *mu0, *mu1, *tmp;
    complex<double> R1D, res;

    mu0 = new complex<double>[n1ex];
    mu1 = new complex<double>[n1ex];
    tmp = new complex<double>[n1ex];

    R1D = complex_zero;

    // get the isotropically averaged spectrum -- x component
    for ( i = 0; i < n1ex; i ++ ) mu0[i] = mu_eg_t0_x[i];
    for ( i = 0; i < n1ex; i ++ ) mu1[i] = mu_eg_t1_x[i];
    cblas_zgemv( CblasRowMajor, CblasNoTrans, n1ex, n1ex,\
                &complex_one, eiH1_t0t1, n1ex, mu0, 1, \
                &complex_zero, tmp, 1 );
    cblas_zdotu_sub( n1ex, mu1, 1, tmp, 1, &res );
    R1D += res;

    // get the isotropically averaged spectrum -- y component
    for ( i = 0; i < n1ex; i ++ ) mu0[i] = mu_eg_t0_y[i];
    for ( i = 0; i < n1ex; i ++ ) mu1[i] = mu_eg_t1_y[i];
    cblas_zgemv( CblasRowMajor, CblasNoTrans, n1ex, n1ex,\
                &complex_one, eiH1_t0t1, n1ex, mu0, 1, \
                &complex_zero, tmp, 1 );
    cblas_zdotu_sub( n1ex, mu1, 1, tmp, 1, &res );
    R1D += res;

    // get the isotropically averaged spectrum -- x component
    for ( i = 0; i < n1ex; i ++ ) mu0[i] = mu_eg_t0_z[i];
    for ( i = 0; i < n1ex; i ++ ) mu1[i] = mu_eg_t1_z[i];
    cblas_zgemv( CblasRowMajor, CblasNoTrans, n1ex, n1ex,\
                &complex_one, eiH1_t0t1, n1ex, mu0, 1, \
                &complex_zero, tmp, 1 );
    cblas_zdotu_sub( n1ex, mu1, 1, tmp, 1, &res );
    R1D += res;

    delete [] mu0;
    delete [] mu1;
    delete [] tmp;

    return R1D;
}

complex<double> IR2D::getR2D( int t1, int t3, string which )
// return the third order response function at a given t1, t3
// See Eq 7.35 from Hamm and Zanni
{
    double mu;
    complex<double> R2D;
    vec3   polx={1.,0.,0.}, poly={0.,1.,0.}, polz={0.,0.,1.};
    double lifetime_T2_12 = 2.*lifetime_T1/3.;// See Hamm and Zanni eq 4.21
    lifetime_T2_12 = lifetime_T2;             // set equal here, see Jansen 2012 

    /*
    R2D = complex_zero;
    // get dipole part for isotropically averaged ZZZZ spectrum
    // NOTE: The Condon approximation is NOT made here
    // all four pulses have the same polarization 
    // see Hamm and Zanni eq 5.35
    mu = 0;
    mu+=dot3(dipole_t0,polx)*dot3(dipole_t1,polx)*
        dot3(dipole_t2,polx)*dot3(dipole_t3,polx);
    mu+=dot3(dipole_t0,poly)*dot3(dipole_t1,poly)*
        dot3(dipole_t2,poly)*dot3(dipole_t3,poly);
    mu+=dot3(dipole_t0,polz)*dot3(dipole_t1,polz)*
        dot3(dipole_t2,polz)*dot3(dipole_t3,polz);

    // get the response function
    if ( which.compare("R1") == 0 ){ // rephasing
        // dipole
        // first pulse oscillating in 01 coherence
        // second pulse population relaxation (note t2 is in ps already)
        // third pulse oscillating in 01 coherence
        R2D += img*mu* 
               conj(eint_t1)*exp(-t1*dt/lifetime_T2)
                            *exp(-t2   /lifetime_T1)
                   *eint_t3 *exp(-t3*dt/lifetime_T2);
    }
    else if ( which.compare("R4") == 0 ){ // non-rephasing
        // dipole
        // first pulse oscillating in 01 coherence
        // second pulse population relaxation (note t2 is in ps already)
        // third pulse oscillating in 01 coherence
        R2D += img*mu* 
               eint_t1*exp(-t1*dt/lifetime_T2)
                      *exp(-t2   /lifetime_T1)
              *eint_t3*exp(-t3*dt/lifetime_T2);
    }
    else if ( which.compare("R3") == 0 ){ // rephasing
        // dipole, the factor of 2 assumes the transition dipoles scale 
        // like a harmonic oscillator (see p 68 of Hamm and Zanni)
        // first pulse oscillating in 01 coherence
        // second pulse population relaxation (note t2 is in ps already)
        // third pulse oscillating in 12 coherence -- see eq 4.21 for relaxation
        // include anharmonicity term for the 12 transition
        R2D -= 2.*img*mu*
               conj(eint_t1)*exp(-t1*dt/lifetime_T2)
                            *exp(-t2   /lifetime_T1)
                   *eint_t3 *exp(-t3*dt/lifetime_T2_12)
                   *exp(img*(dt*t3)*anharm/HBAR);              
    }
    else if ( which.compare("R6") == 0 ){ // non-rephasing
        // dipole, the factor of 2 assumes the transition dipoles scale 
        // like a harmonic oscillator (see p 68 of Hamm and Zanni)
        // first pulse oscillating in 01 coherence
        // second pulse population relaxation (note t2 is in ps already)
        // third pulse oscillating in 12 coherence -- see eq 4.21 for relaxation
        // include anharmonicity term for the 12 transition
        R2D -= 2.*img*mu*
               eint_t1*exp(-t1*dt/lifetime_T2)
                      *exp(-t2   /lifetime_T1)
              *eint_t3*exp(-t3*dt/lifetime_T2_12)
              *exp(img*(dt*t3)*anharm/HBAR);
    }
    else {
        cout << "ERROR:: IR2D::getR2D which= " << which << " is unknown. Aborting." << endl;
        exit(EXIT_FAILURE);
    }
    */

    return R2D;
}

int IR2D::writeR1D()
// write R1D to file
{
    string   fn=_ofile_+"-R1D.dat";
    ofstream ofile;
    int t1;

    ofile.open( fn );
    if ( ! ofile.is_open() ) { fileOpenErr( fn ); return 1;}
    cout << ">>> Writing " << fn << "." << endl;

    ofile << "# time (ps) Real Imag" << endl;
    for ( t1 = 0; t1 < t1t3_npoints; t1 ++ ){
        ofile << t1 * dt << " " << R1D[t1].real() << " " << R1D[t1].imag() << endl;
    }
    ofile.close();

    return IR2DOK;
}

int IR2D::writeR2D()
// write R2D to file
{
    string fn;
    ofstream ofile;
    int t1, t3;
    complex<double> Rtmp;

    /*
    // write rephasing response function
    fn=_ofile_+"-RparI.dat";
    ofile.open( fn );
    if ( ! ofile.is_open() ) { fileOpenErr( fn ); return 1;}
    cout << ">>> Writing " << fn << "." << endl;

    ofile << "# Rephasing parallel ZZZZ polarized response function, t2 = " << t2 << endl;
    ofile << "# t1 (ps) t3 (ps) Real Imag" << endl;
    for ( t1 = 0; t1 < t1t3_npoints; t1 ++ ){
        for ( t3 = 0; t3 < t1t3_npoints; t3 ++ ){
            // see Hamm and Zanni eq 4.23 and note R1=R2, hence the factor of 2
            // The experiment can only see all rephasing or non-rephasing, not the individual
            // response functions, so add them here.
            Rtmp = 2.*R2D_R1[t1 * t1t3_npoints + t3] + R2D_R3[t1 * t1t3_npoints + t3];
            ofile << t1 * dt << " " << t3 * dt << " " << Rtmp.real() << " " << Rtmp.imag() << endl;
        }
    }
    ofile.close();

    // write non-rephasing response function
    fn=_ofile_+"-RparII.dat";
    ofile.open( fn );
    if ( ! ofile.is_open() ) { fileOpenErr( fn ); return 1;}
    cout << ">>> Writing " << fn << "." << endl;

    ofile << "# Non-rephasing parallel ZZZZ polarized response function, t2 = " << t2 << endl;
    ofile << "# t1 (ps) t3 (ps) Real Imag" << endl;
    for ( t1 = 0; t1 < t1t3_npoints; t1 ++ ){
        for ( t3 = 0; t3 < t1t3_npoints; t3 ++ ){
            // see Hamm and Zanni eq 4.23 and note R4=R5, hence the factor of 2
            // The experiment can only see all rephasing or non-rephasing, not the individual
            // response functions, so add them here.
            Rtmp = 2.*R2D_R4[t1 * t1t3_npoints + t3] + R2D_R6[t1 * t1t3_npoints + t3];
            ofile << t1 * dt << " " << t3 * dt << " " << Rtmp.real() << " " << Rtmp.imag() << endl;
        }
    }
    ofile.close();
    */

    return IR2DOK;
}

int IR2D::write1Dfft()
// write fft of R1D
{
    string   fn;
    ofstream ofile;
    complex<double> *fftIn, *fftOut, *res;
    fftw_plan plan;
    double   freq, scale;
    int      t1, i;

    // allocate arrays
    fftIn  = new complex<double>[fftlen]();
    fftOut = new complex<double>[fftlen]();
    res    = new complex<double>[fftlen]();

    // scale fftOut to keep total area of the spectrum conserved and normalize
    // by root(fftlen) since FFTW3 does not
    // the -1 makes it look like absorption spectrum
    scale = dt*fftlen/(2.*PI*HBAR*sqrt(fftlen));
     
    // do fft
    plan = fftw_plan_dft_1d( fftlen, reinterpret_cast<fftw_complex*>(fftIn), \
                                     reinterpret_cast<fftw_complex*>(fftOut),\
                                     FFTW_BACKWARD, FFTW_ESTIMATE );

    if ( 2*t1t3_npoints > fftlen ){
        cout << "ERROR:: fftlen = " << fftlen << " < " << "2*t1t3_max/dt = " 
             << 2*t1t3_npoints << endl;
        cout << "Specify longer fftlen in input file. Aborting." << endl;
        exit(EXIT_FAILURE);
    }
    
    // Absorptive part
    for ( i = 0; i < fftlen ; i ++ ) fftIn[i] = complex_zero;
    for ( t1 = 0; t1 < t1t3_npoints; t1 ++ ){
        fftIn[t1]          = R1D[t1];
        if ( t1 == 0 ) continue; // dont take hermitian at t=0
        fftIn[fftlen-t1] = conj(R1D[t1]);
    }
    fftw_execute(plan); 
    for ( t1 = 0; t1 < fftlen; t1 ++ ) res[t1].real(fftOut[t1].real());

    // Dispersive part
    for ( i = 0; i < fftlen ; i ++ ) fftIn[i] = complex_zero;
    for ( t1 = 0; t1 < t1t3_npoints; t1 ++ ){
        fftIn[t1]          = img*R1D[t1];
        if ( t1 == 0 ) continue; // dont take hermitian at t=0
        fftIn[fftlen-t1] = conj(img*R1D[t1]);
    }
    fftw_execute(plan); 
    fftw_destroy_plan(plan);
    for ( t1 = 0; t1 < fftlen; t1 ++ ) res[t1].imag(fftOut[t1].real());

    // write file
    fn = _ofile_+"-R1Dw.dat";
    ofile.open(fn);
    if ( ! ofile.is_open() ) { fileOpenErr( fn ); return 1;}
    cout << ">>> Writing " << fn << "." << endl;
    ofile << "# Frequency (cm-1) Absorptive Dispersive" << endl;

    // negative frequencies are stored at the end of the fftOut array
    for ( i = fftlen/2; i < fftlen; i ++ )
    {
        // get frequency in wavenumber, add back the shift
        freq   = 2.*PI*HBAR*(i-fftlen)/(dt*fftlen) + shift;
        if ( freq < window0 or freq > window1 ) continue;
        ofile << freq << " " << scale*res[i].real() << " " << scale*res[i].imag() << endl;
    }
    // positive frequencies are stored at the beginning of fftOut
    for ( i = 0; i < fftlen/2; i ++ ){
        // get frequency in wavenumber, add back the shift
        freq   = 2.*PI*HBAR*i/(dt*fftlen) + shift;
        if ( freq < window0 or freq > window1 ) continue;
        ofile << freq << " " << scale*res[i].real() << " " << scale*res[i].imag() << endl;
    }
    ofile.close();

    // delete arrays
    delete [] fftIn;
    delete [] fftOut;
    delete [] res;

    return IR2DOK;
}

int IR2D::write2DRabs()
// write the purely absorptive 2D IR spectrum
{
    string    fn;
    complex<double> *fftIn, *fftOut, *res;
    fftw_plan plan;
    int       t1, t3, i, j;

    /*
    // allocate arrays
    fftIn  = new complex<double>[fftlen*fftlen]();
    fftOut = new complex<double>[fftlen*fftlen]();
    res    = new complex<double>[fftlen*fftlen]();

    // do fft plan
    plan = fftw_plan_dft_2d( fftlen, fftlen, \
                             reinterpret_cast<fftw_complex*>(fftIn), \
                             reinterpret_cast<fftw_complex*>(fftOut),\
                             FFTW_BACKWARD, FFTW_ESTIMATE );
    
    // fourier transform rephasing response functions, see Hamm and Zanni eq 4.31
    // see Hamm and Zanni eq 4.23 and note R1=R2, hence the factor of 2
    for ( i = 0; i < fftlen*fftlen; i ++ ) fftIn[i] = complex_zero;
    for ( t1 = 0; t1 < t1t3_npoints; t1 ++ ){
        for ( t3 = 0; t3 < t1t3_npoints; t3 ++ ){
            fftIn[ t1*fftlen + t3 ] = 2.*img*R2D_R1[ t1*t1t3_npoints + t3 ]\
                                      +  img*R2D_R3[ t1*t1t3_npoints + t3 ];
            // divide t=0 point by 2 (see Hamm and Zanni sec 9.5.3)
            if ( t1 == 0 and t3 == 0 ) fftIn[ t1*fftlen + t3 ] /=2.;
        }
    }
    fftw_execute(plan);
    fn=_ofile_+"-RparIw.dat";
    if ( write2Dout( fftOut, fn, "rephasing", fftlen ) != IR2DOK ) return 1;

    // Save rephasing contribution to purly absorptive spectrum
    // here we map w1 to -w1
    // See Hamm and Zanni eq 4.36
    for ( i = 0; i < fftlen; i ++ ){
        for ( j = 0; j < fftlen; j ++ ){
            if ( i == 0 ) res[ i*fftlen + j ] = fftOut[ i*fftlen + j ]; // w1=0 goes to w1=0
            else          res[ i*fftlen + j ] = fftOut[ (fftlen - i)*fftlen + j ]; // w1 to -w1
        }
    }

    // fourier transform non-rephasing response functions, see Hamm and Zanni eq 4.31
    // see Hamm and Zanni eq 4.23 and note R4=R5, hence the factor of 2
    for ( i = 0; i < fftlen*fftlen; i ++ ) fftIn[i] = complex_zero;
    for ( t1 = 0; t1 < t1t3_npoints; t1 ++ ){
        for ( t3 = 0; t3 < t1t3_npoints; t3 ++ ){
            fftIn[ t1*fftlen + t3 ] = 2.*img*R2D_R4[ t1*t1t3_npoints + t3 ]\
                                      +  img*R2D_R6[ t1*t1t3_npoints + t3 ];
            // divide t=0 point by 2 (see Hamm and Zanni sec 9.5.3)
            if ( t1 == 0 and t3 == 0 ) fftIn[ t1*fftlen + t3 ] /=2.;
        }
    }
    fftw_execute(plan);
    fftw_destroy_plan(plan);
    fn=_ofile_+"-RparIIw.dat";
    if ( write2Dout( fftOut, fn, "non-rephasing", fftlen ) != IR2DOK ) return 1;
    
    // Save rephasing contribution to purly absorptive spectrum
    // here w1 goes to w1
    // See Hamm and Zanni eq 4.36
    for ( i = 0; i < fftlen; i ++ ){
        for ( j = 0; j < fftlen; j ++ ){
            res[ i*fftlen + j ] += fftOut[ i*fftlen + j ];
        }
    }

    // write purly absorptive spectrum 
    fn=_ofile_+"-RparAbs.dat";
    if ( write2Dout( res, fn, "rabs", fftlen ) != IR2DOK ) return 1;

    // delete arrays
    delete [] fftIn;
    delete [] fftOut;
    delete [] res;
    */

    return IR2DOK;
}


int IR2D::write2Dout( complex<double> *data, string fn, string which, int n )
// write 2D fourier transformed output
{
    int i, j;
    double shift_w1, shift_w3, window0_w1, window0_w3, window1_w1, window1_w3;
    double scale, w1, w3;
    ofstream ofile;

    // scaling
    scale = dt*n/(2.*PI*HBAR*sqrt(n));
    scale*= scale;  // since 2D scaling

    ofile.open( fn );
    if ( ! ofile.is_open() ) { fileOpenErr( fn ); return 1;}
    cout << ">>> Writing " << fn << "." << endl;

    // assign shifts and spectral window limits
    shift_w1   = shift;
    shift_w3   = shift;
    window0_w1 = window0;
    window1_w1 = window1;
    window0_w3 = window0;
    window1_w3 = window1;

    if ( which.compare("rephasing") == 0 ){
        ofile << "# Rephasing parallel ZZZZ polarized response function, t2 = " << t2 << endl;
        shift_w1   = -shift;
        window0_w1 = -window1;
        window1_w1 = -window0;
    }
    else if ( which.compare("non-rephasing") == 0 ){
        ofile << "# Non-rephasing parallel ZZZZ polarized response function, t2 = " << t2 << endl;
    }
    else if ( which.compare("rabs") == 0 ){
        ofile << "# Pure parallel ZZZZ polarized absorption, t2 = " << t2 << endl;
    }
    else{
        cout << "ERROR:: write2Dout which= " << which << " unknown." << endl;
        return 1;
    }

    ofile << "# w1 (cm-1) w3 (cm-1) Real Imag" << endl;
    // negative frequencies are stored at the end of data array
    for ( i = n/2; i < n; i ++ ){
        w1 = 2.*PI*HBAR*(i-n)/(dt*n) + shift_w1;
        if ( w1 < window0_w1 or w1 > window1_w1 ) continue;
        for ( j = n/2; j < n; j ++ ){
            w3 = 2.*PI*HBAR*(j-n)/(dt*n) + shift_w3;
            if ( w3 < window0_w3 or w3 > window1_w3 ) continue;
            ofile << w1 << " " << w3 << " " << scale*data[i*n+j].real() 
                        << " " << scale*data[i*n+j].imag() << endl;
        }
        for ( j = 0; j < n/2; j ++ ){
            w3 = 2.*PI*HBAR*j/(dt*n) + shift_w3;
            if ( w3 < window0_w3 or w3 > window1_w3 ) continue;
            ofile << w1 << " " << w3 << " " << scale*data[i*n+j].real() 
                        << " " << scale*data[i*n+j].imag() << endl;
        }
    }
    // positive frequencies are stored at the beginning of data array
    for ( i = 0; i < n/2; i ++ ){
        w1 = 2.*PI*HBAR*i/(dt*n) + shift_w1;
        if ( w1 < window0_w1 or w1 > window1_w1 ) continue;
        for ( j = n/2; j < n; j ++ ){
            w3 = 2.*PI*HBAR*(j-n)/(dt*n) + shift_w3;
            if ( w3 < window0_w3 or w3 > window1_w3 ) continue;
            ofile << w1 << " " << w3 << " " << scale*data[i*n+j].real() 
                  << " " << scale*data[i*n+j].imag() << endl;
        }
        for ( j = 0; j < n/2; j ++ ){
            w3 = 2.*PI*HBAR*j/(dt*n) + shift_w3;
            if ( w3 < window0_w3 or w3 > window1_w3 ) continue;
            ofile << w1 << " " << w3 << " " << scale*data[i*n+j].real() 
                  << " " << scale*data[i*n+j].imag() << endl;
        }
    }
    ofile.close();

    return IR2DOK;
}

void printProgress( int currentStep, int totalSteps )
// print a progress bar to keep updated on usage
{
    float percentage = (float) currentStep / (float) totalSteps;
    int lpad = (int) (percentage*PWID);
    int rpad = PWID - lpad;
    fprintf(stderr, "\r [%.*s%*s]%3d%%", lpad, PSTR, rpad, "",(int) (percentage*100));
}

int main( int argc, char* argv[] )
{
    int sample;
    int it0, it1, it2, it3;
    int it1_max, it2_max, it3_max;

    // check user input to program
    if ( argc != 2 ){
        cout << "ERROR:: Program expects the name of the input" << 
                "file as the only argument. Aborting." << endl;
        exit( EXIT_FAILURE );
    }

    // get input file name and initialize IR2D class
    IR2D spectrum( argv[1] ); 

    // set t2 waiting time in units of frames
    it2_max = static_cast<int>(spectrum.t2/spectrum.dt);

    // Loop over the trajectory
    for ( sample = 0; sample < spectrum.nsamples; sample ++ ){

        // get frame number, read and save dipole at t0
        it0 = sample*static_cast<int>(spectrum.sample_every/spectrum.dt);
        fprintf(stderr, "    Now processing sample %d/%d starting at %.2f ps\n", \
                sample+1, spectrum.nsamples, it0*spectrum.dt ); fflush(stderr);
        if ( spectrum.readDframe( it0 ) != IR2DOK ) exit(EXIT_FAILURE);
        if ( spectrum.setMUatT( "t0" )  != IR2DOK ) exit(EXIT_FAILURE);

        // loop over t1
        it1_max = it0 + spectrum.t1t3_npoints;
        for ( it1 = it0; it1 < it1_max; it1 ++ ){

            // read in energy and dipole
            if ( spectrum.readEframe(it1) != IR2DOK ) exit(EXIT_FAILURE);
            if ( spectrum.propigateH1( it0, it1, "t0-t1" ) != IR2DOK ) exit(EXIT_FAILURE);
            if ( spectrum.readDframe(it1) != IR2DOK ) exit(EXIT_FAILURE);
            if ( spectrum.setMUatT("t1")  != IR2DOK ) exit(EXIT_FAILURE);

            // get exponential integral and 1D response function
            spectrum.R1D[it1-it0] += spectrum.getR1D();

            // get frame number and dipole for time t2
            /*
            frame_t2 = frame_t1 + it2;
            if ( spectrum.readDframe(it2) != IR2DOK ) exit(EXIT_FAILURE);
            if ( spectrum.setMUatT("t2")  != IR2DOK ) exit(EXIT_FAILURE);
            */


            // loop over t3 at current t1
            /*
            for ( it3 = 0; it3 < spectrum.t1t3_npoints; it3 ++ ){
                // get frame number for current t3
                frame_t3 = frame_t2 + it3;

                // read in energy and dipole
                if ( spectrum.readEframe(frame_t3, "t3") != IR2DOK ) exit(EXIT_FAILURE);
                if ( spectrum.readDframe(frame_t3, "t3") != IR2DOK ) exit(EXIT_FAILURE); 

                // get exponential integral and 2D response function 
                if ( spectrum.get_eint(t3, "t3") != IR2DOK ) exit(EXIT_FAILURE);
                spectrum.R2D_R1[ it1 * spectrum.t1t3_npoints + it3 ] += spectrum.getR2D(it1, it3, "R1" );
                spectrum.R2D_R2[ it1 * spectrum.t1t3_npoints + it3 ] += spectrum.getR2D(it1, it3, "R2" );
            }
            */
            printProgress( it1+1, spectrum.t1t3_npoints );
        }
        cerr << endl;
    }

    // normalize response functions by number of samples and add phenomelogical decay
    for ( it1 = 0; it1 < spectrum.t1t3_npoints; it1 ++ ){
        spectrum.R1D[it1]*=exp(-it1*spectrum.dt/spectrum.lifetime_T2)/(1.*spectrum.nsamples);
        /*
        for ( it3 = 0; it3 < spectrum.t1t3_npoints; it3 ++ ){
            spectrum.R2D_R1[ it1 * spectrum.t1t3_npoints + it3 ] /= ( 1.*spectrum.nsamples );
            spectrum.R2D_R2[ it1 * spectrum.t1t3_npoints + it3 ] /= ( 1.*spectrum.nsamples );
        }
        */
    }

    // do fourier transforms and write them out
    if ( spectrum.writeR1D()    != IR2DOK ) exit(EXIT_FAILURE);
    //if ( spectrum.writeR2D()    != IR2DOK ) exit(EXIT_FAILURE);
    if ( spectrum.write1Dfft()  != IR2DOK ) exit(EXIT_FAILURE);
    //if ( spectrum.write2DRabs() != IR2DOK ) exit(EXIT_FAILURE);
    cout << ">>> Done!" << endl;
}

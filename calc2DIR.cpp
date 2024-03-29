#include <string.h>
#include <fstream>
#include <sstream>
#include <iostream>
#include <iomanip>
#include <stdio.h>
#include <complex.h>
#include <fftw3.h>
#include <math.h>
#include <algorithm>
#include "calc2DIR.h"
#define MKL_Complex16 complex<double>
#include <mkl.h>
#include <omp.h>
#include <pthread.h>

using namespace std;

IR2D::IR2D( string _inpf_ )
// Default Constructor
{
    string line;

    // read input file
    if ( readParam( _inpf_ ) != IR2DOK ) exit( EXIT_FAILURE );

    // allocate variable arrays
    // response functions
    R1D             = new complex<double>[t1t3_npoints]();
    R2D_R1          = new complex<double>[t1t3_npoints*t1t3_npoints]();
    R2D_R2          = new complex<double>[t1t3_npoints*t1t3_npoints]();

    // one exciton
    H1              = new double [n1ex*n1ex]();
    eiH1_t0t1       = new complex<double>[n1ex*n1ex]();
    eiH1_t1t2       = new complex<double>[n1ex*n1ex]();
    eiH1_t2t3       = new complex<double>[n1ex*n1ex]();
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
    eiH1_t0t1_mu_eg_t0_x = new complex<double>[n1ex*t1t3_npoints]();
    eiH1_t0t1_mu_eg_t0_y = new complex<double>[n1ex*t1t3_npoints]();
    eiH1_t0t1_mu_eg_t0_z = new complex<double>[n1ex*t1t3_npoints]();

    // two exciton
    H2              = new double [n2ex*n2ex]();
    eiH2_t1t2       = new complex<double>[n2ex*n2ex]();
    eiH2_t2t3       = new complex<double>[n2ex*n2ex]();
    mu_ce_x         = new double[n1ex*n2ex]();
    mu_ce_y         = new double[n1ex*n2ex]();
    mu_ce_z         = new double[n1ex*n2ex]();
    mu_ce_t2_x      = new double[n1ex*n2ex]();
    mu_ce_t3_x      = new double[n1ex*n2ex]();
    mu_ce_t2_y      = new double[n1ex*n2ex]();
    mu_ce_t3_y      = new double[n1ex*n2ex]();
    mu_ce_t2_z      = new double[n1ex*n2ex]();
    mu_ce_t3_z      = new double[n1ex*n2ex]();

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
    // response functions
    delete [] R1D;
    delete [] R2D_R1;
    delete [] R2D_R2;

    // one-exciton
    delete [] H1;
    delete [] eiH1_t0t1;
    delete [] eiH1_t1t2;
    delete [] eiH1_t2t3;
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
    delete [] eiH1_t0t1_mu_eg_t0_x;
    delete [] eiH1_t0t1_mu_eg_t0_y;
    delete [] eiH1_t0t1_mu_eg_t0_z;

    // two-exciton
    delete [] H2;
    delete [] eiH2_t1t2;
    delete [] eiH2_t2t3;
    delete [] mu_ce_x;
    delete [] mu_ce_y;
    delete [] mu_ce_z;
    delete [] mu_ce_t2_x;
    delete [] mu_ce_t3_x;
    delete [] mu_ce_t2_y;
    delete [] mu_ce_t3_y;
    delete [] mu_ce_t2_z;
    delete [] mu_ce_t3_z;

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
        else if ( para.compare("anharm")      == 0 )  anharm       = stof(value);
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
    tellParam<double>( "anharm", anharm);
    tellParam<int>( "nsamples", nsamples );
    tellParam<double>( "sample_every", sample_every );
    tellParam<int>( "fftlen", fftlen );
    tellParam<double>( "window0", window0 );
    tellParam<double>( "window1", window1 );
    cout << ">>> Done reading simulation parameters from " << _inpf_ << endl;

    if ( trjlen < static_cast<int>(round(((nsamples-1)*sample_every + (2*t1t3_max + t2))/dt)) ){
        cout << "WARNING:: The given trajectory length is not long enough.\n" << 
                      "\t  Must be " <<  ((nsamples-1)*sample_every + (2*t1t3_max + t2))/dt << 
                      " frames long.\n\t  Check input file. Aborting." << endl;
        exit(EXIT_FAILURE);
    }

    // set number of points in response functions based on t1t3_max and dt
    t1t3_npoints = static_cast<int>(round(t1t3_max/dt));

    // set number of one- and two-exciton states
    n1ex = nchrom;
    n2ex = nchrom*(nchrom+1)/2;

    cout << ">>> There are " << n1ex << " one-exciton states and " \
         << n2ex << " two-exciton states." << endl;

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
    int     frameTmp, nelm, i, j, inx;
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
        inx = i*n1ex - i*(i+1)/2; // the indexing is a little unintuitive since 
                                  // only the upper tridiagonal is given as input
    for( j = i; j < n1ex; j++ ){
        H1[i*n1ex + j] = Htmp[inx  + j];
        H1[j*n1ex + i] = H1[i*n1ex + j];
    }
    }

    // build the 2-exciton hamiltonian
    // do it here to avoid possible problems doing it later since H1 can change
    if ( buildH2() != IR2DOK ) return 2;

    // shift energies so eigenvalues are close to zero, these are added back later
    // after taking the fft
    for ( i = 0; i < n1ex; i ++ ) H1[i*n1ex + i] -= shift;
    for ( i = 0; i < n2ex; i ++ ) H2[i*n2ex + i] -= 2.*shift;

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

    // build the one-exciton to two-exciton transition dipole matrix
    if ( build_ce_mu() != IR2DOK ) return 2;

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
        memcpy( mu_ce_t2_x, mu_ce_x, sizeof(double)*n1ex*n2ex );
        memcpy( mu_ce_t2_y, mu_ce_y, sizeof(double)*n1ex*n2ex );
        memcpy( mu_ce_t2_z, mu_ce_z, sizeof(double)*n1ex*n2ex );
    }
    else if ( which.compare("t3") == 0 ){
        memcpy( mu_eg_t3_x, mu_eg_x, sizeof(double)*n1ex );
        memcpy( mu_eg_t3_y, mu_eg_y, sizeof(double)*n1ex );
        memcpy( mu_eg_t3_z, mu_eg_z, sizeof(double)*n1ex );
        memcpy( mu_ce_t3_x, mu_ce_x, sizeof(double)*n1ex*n2ex );
        memcpy( mu_ce_t3_y, mu_ce_y, sizeof(double)*n1ex*n2ex );
        memcpy( mu_ce_t3_z, mu_ce_z, sizeof(double)*n1ex*n2ex );
    }

    else{
        cout << "ERROR:: IR2D::setMUatT which= " << which << " is unknown." << endl;
        return 1;
    }

    return IR2DOK;
}

int IR2D::get2nx(int i,int j)
// return index for 2-exciton states
{
    // since |i,j>=|j,i> have to be careful to not overcount
    if ( i <=j ){
        return n1ex*i - i*(i+1)/2 + j;
    }
    else{
        return n1ex*j - j*(j+1)/2 + i;
    }
}

int IR2D::buildH2()
// build the two-exciton hamiltonian from the one-exciton hamiltonian
{
    int i, j, k, l;
    int nx, mx;

    // set to zero
    for ( i = 0; i < n2ex*n2ex; i ++ ) H2[i] = 0.;

    for ( i = 0; i < n1ex; i ++ ){
        // diagonal energies for double excited states <i,i|H|i,i>
        nx = get2nx( i, i ); // |i,i>
        H2[ nx*n2ex + nx ] = 2*H1[ i*n1ex + i ] - anharm; // <i,i|H|i,i>
    for ( j = i+1; j < n1ex; j ++ ){
        // diagonal energies for states with two single excitations <i,j|H|i,j>
        nx = get2nx( i, j ); // |i,j>
        H2[ nx*n2ex + nx ] = H1[ i*n1ex + i ] + H1[ j*n1ex + j ]; // <i,j|H|i,j>
    }
    }

    // off diagonal couplings between single and double excited states <i,i|H|i,j>
    // make harmonic approximation for the couplings --> gives sqrt(2) factor
    // see Hamm and Zanni chapter 6.1
    for ( i = 0; i < n1ex; i ++ ){
        nx = get2nx(i,i);
    for ( j = 0; j < n1ex; j ++ ){
        mx = get2nx(i,j);
        if ( nx == mx ) continue; // diagonal done above
        H2[ nx*n2ex + mx ] = sqrt(2.0)*H1[ i*n1ex + j ];    // <i,i|H|i,j>
        H2[ mx*n2ex + nx ] = H2[ nx*n2ex + mx ];            // c.c.
    }
    }
        
    // off diagonal couplings between singly excited states <i,j|H|i,k>
    for ( i = 0; i < n1ex; i ++ ){
    for ( j = i+1; j < n1ex; j ++ ){
        nx = get2nx(i,j);
        for ( k = 0; k < n1ex; k ++ ){
            if ( i == k ) continue; // double excited states done above
            mx = get2nx(i,k);
            if ( nx == mx ) continue; // diagonal done above
            H2[ nx*n2ex + mx ] = H1[ j*n1ex + k ];      // <i,j|H|i,k>
            H2[ mx*n2ex + nx ] = H2[ nx*n2ex + mx ];    // c.c.
        }
        for ( k = 0; k < n1ex; k ++ ){
            if ( j == k ) continue; // double excited states done above
            mx = get2nx(j,k);
            if ( nx == mx ) continue; // diagonal done above
            H2[ nx*n2ex + mx ] = H1[ i*n1ex + k ];      // <i,j|H|k,j>
            H2[ mx*n2ex + nx ] = H2[ nx*n2ex + mx ];    // c.c.
    }
    }
    }

    return IR2DOK;
}

int IR2D::build_ce_mu()
// build the two-exciton to one-exciton transition dipole moment
{
    int i, j, k, mx;

    // zero
    for ( i = 0; i < n1ex*n2ex; i ++ ){
        mu_ce_x[i] = 0.;
        mu_ce_y[i] = 0.;
        mu_ce_z[i] = 0.;
    }

    // use harmonic rules to build transition dipole moment matrix
    // this is where the sqrt(2) factor comes from
    for ( i = 0; i < n1ex; i ++ ){
        mx = get2nx(i,i);
        mu_ce_x[i*n2ex + mx] = sqrt(2.0)*mu_eg_x[i]; // <i|u|i,i>
        mu_ce_y[i*n2ex + mx] = sqrt(2.0)*mu_eg_y[i]; // <i|u|i,i>
        mu_ce_z[i*n2ex + mx] = sqrt(2.0)*mu_eg_z[i]; // <i|u|i,i>
    for ( j = 0; j < n1ex; j ++ ){
        if ( j == i ) continue; // did single to double excitations above
        mx = get2nx(i,j);
        mu_ce_x[i*n2ex + mx] = mu_eg_x[j]; // <i|u|i,j>
        mu_ce_y[i*n2ex + mx] = mu_eg_y[j]; // <i|u|i,j>
        mu_ce_z[i*n2ex + mx] = mu_eg_z[j]; // <i|u|i,j>
    }
    }

    return IR2DOK;
}

int IR2D::reset_eiH1( string which )
// propigate the one exciton hamiltonian to the next time frame
{
    int i;
    
    // reset eiH1 to complex_one
    if ( which.compare("t0-t1") == 0 ){
        for ( i = 0; i < n1ex*n1ex; i ++ ) eiH1_t0t1[i]   = complex_zero;
        for ( i = 0; i < n1ex; i ++ ) eiH1_t0t1[i*n1ex+i] = complex_one;
    }
    else if ( which.compare("t1-t2") == 0 ){
        for ( i = 0; i < n1ex*n1ex; i ++ ) eiH1_t1t2[i]   = complex_zero;
        for ( i = 0; i < n1ex; i ++ ) eiH1_t1t2[i*n1ex+i] = complex_one;
    }
    else if ( which.compare("t2-t3") == 0 ){
        for ( i = 0; i < n1ex*n1ex; i ++ ) eiH1_t2t3[i]   = complex_zero;
        for ( i = 0; i < n1ex; i ++ ) eiH1_t2t3[i*n1ex+i] = complex_one;
    }
    else{
        cout << "ERROR:: IR2D::reset_eiH1 which= " << which << \
                " is unknown." << endl; return 1;
    }

    return IR2DOK;
}

int IR2D::reset_eiH2( string which )
// propigate the two exciton hamiltonian to the next time frame
{
    int i;
    
    // reset eiH1 to complex_one
    if ( which.compare("t1-t2") == 0 ){
        for ( i = 0; i < n2ex*n2ex; i ++ ) eiH2_t1t2[i]   = complex_zero;
        for ( i = 0; i < n2ex; i ++ ) eiH2_t1t2[i*n2ex+i] = complex_one;
    }
    else if ( which.compare("t2-t3") == 0 ){
        for ( i = 0; i < n2ex*n2ex; i ++ ) eiH2_t2t3[i]   = complex_zero;
        for ( i = 0; i < n2ex; i ++ ) eiH2_t2t3[i*n2ex+i] = complex_one;
    }
    else{
        cout << "ERROR:: IR2D::reset_eiH2 which= " << which << \
                " is unknown." << endl; return 1;
    }

    return IR2DOK;
}

int IR2D::prop_eiH1( int it0, string which )
// propigate the one exciton hamiltonian to the next time frame
{
    int i;
    complex<double> *eiH1, *work;
    
    eiH1 = new complex<double>[n1ex*n1ex];
    work = new complex<double>[n1ex*n1ex];
    
    // read current energy frame and then
    // deterimine e^-iH1dt/(2*Hbar), return to eiH1
    readEframe(it0);
    doeiH(eiH1,H1,n1ex);

    // do the integration
    if ( which.compare("t0-t1") == 0 ){
        // add current term to the propigator -- note that in the main
        // routine t0 iterates backwards from t0_max = t1, so the current
        // propigator is multiplied on the right
        // work=eiH1_t0t1*eiH1
        cblas_zgemm( CblasRowMajor, CblasNoTrans, CblasNoTrans, \
                     n1ex, n1ex, n1ex, &complex_one, eiH1_t0t1, n1ex, \
                     eiH1, n1ex, &complex_zero, work, n1ex );
        memcpy( eiH1_t0t1, work, sizeof(complex<double>)*n1ex*n1ex );
    }
    else if ( which.compare("t1-t2") == 0 ){
        // work=eiH1*eiH1_t1t2 -- note that in the main routine t2 iterates
        // forward from t1_max to t2_max, so the current
        // propigator is multiplied on the left
        // work=eiH1*eiH1_t1t2
        cblas_zgemm( CblasRowMajor, CblasNoTrans, CblasNoTrans, \
                     n1ex, n1ex, n1ex, &complex_one, eiH1, n1ex, \
                     eiH1_t1t2, n1ex, &complex_zero, work, n1ex );
        memcpy( eiH1_t1t2, work, sizeof(complex<double>)*n1ex*n1ex );
    }
    else if ( which.compare("t2-t3") == 0 ){
        // work=eiH1*eiH1_t2t3 -- note that in the main routine t3 iterates
        // forward from t2_max to t3_max, so the current
        // propigator is multiplied on the left
        // work=eiH1*eiH1_t2t3
        cblas_zgemm( CblasRowMajor, CblasNoTrans, CblasNoTrans, \
                     n1ex, n1ex, n1ex, &complex_one, eiH1, n1ex, \
                     eiH1_t2t3, n1ex, &complex_zero, work, n1ex );
        memcpy( eiH1_t2t3, work, sizeof(complex<double>)*n1ex*n1ex );
    }
    else{
        cout << "ERROR:: IR2D::prop_eiH1 which= " << which << \
                " is unknown." << endl; return 1;
    }

    delete [] eiH1;
    delete [] work;
    return IR2DOK;
}

int IR2D::save_eiH1_t0t1_mu0(int it0 )
// save eiH1_t0t1*mu(t0) for all t0
// note, it0 is relative time; it0 = it1-it0
{
    int i;
    complex<double> *work, *mu;
    
    work = new complex<double>[n1ex]();
    mu   = new complex<double>[n1ex]();
    
    // get mu_eg_0 at it0 
    // x component
    for ( i = 0; i < n1ex; i ++ ) mu[i] = complex_one*mu_eg_t0_x[i];
    // work=eiH1_t0t1*mu
    cblas_zgemv( CblasRowMajor, CblasNoTrans, n1ex, n1ex,\
                 &complex_one, eiH1_t0t1, n1ex, mu, 1, \
                 &complex_zero, work, 1 );
    memcpy( &eiH1_t0t1_mu_eg_t0_x[it0*n1ex], work, sizeof(complex<double>)*n1ex );

    // y component
    for ( i = 0; i < n1ex; i ++ ) mu[i] = complex_one*mu_eg_t0_y[i];
    // work=eiH1_t0t1*mu
    cblas_zgemv( CblasRowMajor, CblasNoTrans, n1ex, n1ex,\
                 &complex_one, eiH1_t0t1, n1ex, mu, 1, \
                 &complex_zero, work, 1 );
    memcpy( &eiH1_t0t1_mu_eg_t0_y[it0*n1ex], work, sizeof(complex<double>)*n1ex );

    // z component
    for ( i = 0; i < n1ex; i ++ ) mu[i] = complex_one*mu_eg_t0_z[i];
    // work=eiH1_t0t1*mu
    cblas_zgemv( CblasRowMajor, CblasNoTrans, n1ex, n1ex,\
                 &complex_one, eiH1_t0t1, n1ex, mu, 1, \
                 &complex_zero, work, 1 );
    memcpy( &eiH1_t0t1_mu_eg_t0_z[it0*n1ex], work, sizeof(complex<double>)*n1ex );

    delete [] work;
    delete [] mu;

    return IR2DOK;
}

int IR2D::prop_eiH2( int it0, string which )
// propigate the one exciton hamiltonian to the next time frame
{
    int i;
    complex<double> *eiH2, *work;
    
    eiH2 = new complex<double>[n2ex*n2ex];
    work = new complex<double>[n2ex*n2ex];
    
    // read current energy frame and then
    // deterimine e^-iH2dt/(2*Hbar), return to eiH1
    readEframe(it0);
    doeiH(eiH2,H2,n2ex);

    // do the integration
    if ( which.compare("t1-t2") == 0 ){
        // work=eiH2*eiH2_t1t2 -- note that in the main routine t2 iterates
        // forward from t1_max to t2_max, so the current propigator 
        // is multiplied on the left
        // work=eiH2*eiH2_t1t2
        cblas_zgemm( CblasRowMajor, CblasNoTrans, CblasNoTrans, \
                     n2ex, n2ex, n2ex, &complex_one, eiH2, n2ex, \
                     eiH2_t1t2, n2ex, &complex_zero, work, n2ex );
        memcpy( eiH2_t1t2, work, sizeof(complex<double>)*n2ex*n2ex );
    }
    else if ( which.compare("t2-t3") == 0 ){
        // work=eiH2*eiH2_t2t3 -- note that in the main routine t3 iterates
        // forward from t2_max to t3_max, so the current propigator 
        // is multiplied on the left
        // work=eiH2*eiH2_t2t3
        cblas_zgemm( CblasRowMajor, CblasNoTrans, CblasNoTrans, \
                     n2ex, n2ex, n2ex, &complex_one, eiH2, n2ex, \
                     eiH2_t2t3, n2ex, &complex_zero, work, n2ex );
        memcpy( eiH2_t2t3, work, sizeof(complex<double>)*n2ex*n2ex );
    }
    else{
        cout << "ERROR:: IR2D::prop_eiH2 which= " << which << \
                " is unknown." << endl; return 1;
    }

    delete [] eiH2;
    delete [] work;
    return IR2DOK;
}

int IR2D::doeiH( complex<double> *eiH, double *H, int N )
// get the exponential of exp(i H dt/hbar)
{
    complex<double> *evec, *work, arg;
    double *W;
    int i = 0, j = 0, info;

    // allocate arrays
    evec = new complex<double>[N*N]();
    work = new complex<double>[N*N]();
    W    = new double[N]();

    // get eigenvalues and eigenvector of H
    if ( info = LAPACKE_dsyevd( LAPACK_ROW_MAJOR, 'V', 'U', N, H, N, W ) != 0 ){
        cout << "WARNING:: LAPACKE_dsyeved returned " << info << endl;
        return 1;
    }

    // exponentiate eigenvalues
    for ( i = 0; i < N; i ++ ){
        arg = -img*W[i]*dt/HBAR;
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
    // work = evec*eiH
    cblas_zgemm( CblasRowMajor, CblasNoTrans, CblasNoTrans, \
                 N, N, N, &complex_one, evec, N, eiH, N, \
                 &complex_zero, work , N );
    // eiH = work*evec'
    cblas_zgemm( CblasRowMajor, CblasNoTrans, CblasTrans, \
                 N, N, N, &complex_one, work, N, evec, N, \
                 &complex_zero, eiH, N );

    // delete arrays
    delete [] evec;
    delete [] work;
    delete [] W;

    return IR2DOK;
}

complex<double> IR2D::getR1D(int it0)
// return the linear response function
{
    int i, k;
    complex<double> *mu0, *mu1;
    complex<double> R1D, res;

    mu0  = new complex<double>[n1ex];
    mu1  = new complex<double>[n1ex];

    R1D = complex_zero;

    // get the isotropic average
    // x component
    for ( i = 0; i < n1ex; i ++ ) mu0[i] = eiH1_t0t1_mu_eg_t0_x[it0*n1ex+i];
    for ( i = 0; i < n1ex; i ++ ) mu1[i] = mu_eg_t1_x[i];
    // res=mu1*mu0
    cblas_zdotu_sub( n1ex, mu1, 1, mu0, 1, &res );
    R1D += res;

    // y component
    for ( i = 0; i < n1ex; i ++ ) mu0[i] = eiH1_t0t1_mu_eg_t0_y[it0*n1ex+i];
    for ( i = 0; i < n1ex; i ++ ) mu1[i] = mu_eg_t1_y[i];
    // res=mu1*mu0
    cblas_zdotu_sub( n1ex, mu1, 1, mu0, 1, &res );
    R1D += res;

    // z component
    for ( i = 0; i < n1ex; i ++ ) mu0[i] = eiH1_t0t1_mu_eg_t0_z[it0*n1ex+i];
    for ( i = 0; i < n1ex; i ++ ) mu1[i] = mu_eg_t1_z[i];
    // res=mu1*mu0
    cblas_zdotu_sub( n1ex, mu1, 1, mu0, 1, &res );
    R1D += res;

    delete [] mu0;
    delete [] mu1;

    return R1D;
}

int IR2D::getMU4pol( complex<double> *mu0_eg, complex<double> *mu1_eg,\
                     complex<double> *mu2_eg, complex<double> *mu3_eg,\
                     complex<double> *mu2_ce, complex<double> *mu3_ce,\
                     int it0, int k )
// return the relevant set of dipole matrices for a given polarization combination
{
    int i;

    if ( k == 0 ){ // XXXX
        for ( i = 0; i < n1ex; i ++ )      mu0_eg[i] = eiH1_t0t1_mu_eg_t0_x[it0*n1ex+i];
        for ( i = 0; i < n1ex; i ++ )      mu1_eg[i] = mu_eg_t1_x[i];
        for ( i = 0; i < n1ex; i ++ )      mu2_eg[i] = mu_eg_t2_x[i];
        for ( i = 0; i < n1ex; i ++ )      mu3_eg[i] = mu_eg_t3_x[i];
        for ( i = 0; i < n1ex*n2ex; i ++ ) mu2_ce[i] = mu_ce_t2_x[i];
        for ( i = 0; i < n1ex*n2ex; i ++ ) mu3_ce[i] = mu_ce_t3_x[i];
    }
    else if ( k == 1 ){ // YYYY
        for ( i = 0; i < n1ex; i ++ )      mu0_eg[i] = eiH1_t0t1_mu_eg_t0_y[it0*n1ex+i];
        for ( i = 0; i < n1ex; i ++ )      mu1_eg[i] = mu_eg_t1_y[i];
        for ( i = 0; i < n1ex; i ++ )      mu2_eg[i] = mu_eg_t2_y[i];
        for ( i = 0; i < n1ex; i ++ )      mu3_eg[i] = mu_eg_t3_y[i];
        for ( i = 0; i < n1ex*n2ex; i ++ ) mu2_ce[i] = mu_ce_t2_y[i];
        for ( i = 0; i < n1ex*n2ex; i ++ ) mu3_ce[i] = mu_ce_t3_y[i];
    }
    else if ( k == 2 ){ // ZZZZ
        for ( i = 0; i < n1ex; i ++ )      mu0_eg[i] = eiH1_t0t1_mu_eg_t0_z[it0*n1ex+i];
        for ( i = 0; i < n1ex; i ++ )      mu1_eg[i] = mu_eg_t1_z[i];
        for ( i = 0; i < n1ex; i ++ )      mu2_eg[i] = mu_eg_t2_z[i];
        for ( i = 0; i < n1ex; i ++ )      mu3_eg[i] = mu_eg_t3_z[i];
        for ( i = 0; i < n1ex*n2ex; i ++ ) mu2_ce[i] = mu_ce_t2_z[i];
        for ( i = 0; i < n1ex*n2ex; i ++ ) mu3_ce[i] = mu_ce_t3_z[i];
    }
    else if ( k == 3 ){ // XXYY
        for ( i = 0; i < n1ex; i ++ )      mu0_eg[i] = eiH1_t0t1_mu_eg_t0_x[it0*n1ex+i];
        for ( i = 0; i < n1ex; i ++ )      mu1_eg[i] = mu_eg_t1_x[i];
        for ( i = 0; i < n1ex; i ++ )      mu2_eg[i] = mu_eg_t2_y[i];
        for ( i = 0; i < n1ex; i ++ )      mu3_eg[i] = mu_eg_t3_y[i];
        for ( i = 0; i < n1ex*n2ex; i ++ ) mu2_ce[i] = mu_ce_t2_y[i];
        for ( i = 0; i < n1ex*n2ex; i ++ ) mu3_ce[i] = mu_ce_t3_y[i];
    }
    else if ( k == 4 ){ // XXZZ
        for ( i = 0; i < n1ex; i ++ )      mu0_eg[i] = eiH1_t0t1_mu_eg_t0_x[it0*n1ex+i];
        for ( i = 0; i < n1ex; i ++ )      mu1_eg[i] = mu_eg_t1_x[i];
        for ( i = 0; i < n1ex; i ++ )      mu2_eg[i] = mu_eg_t2_z[i];
        for ( i = 0; i < n1ex; i ++ )      mu3_eg[i] = mu_eg_t3_z[i];
        for ( i = 0; i < n1ex*n2ex; i ++ ) mu2_ce[i] = mu_ce_t2_z[i];
        for ( i = 0; i < n1ex*n2ex; i ++ ) mu3_ce[i] = mu_ce_t3_z[i];
    }
    else if ( k == 5 ){ // YYXX
        for ( i = 0; i < n1ex; i ++ )      mu0_eg[i] = eiH1_t0t1_mu_eg_t0_y[it0*n1ex+i];
        for ( i = 0; i < n1ex; i ++ )      mu1_eg[i] = mu_eg_t1_y[i];
        for ( i = 0; i < n1ex; i ++ )      mu2_eg[i] = mu_eg_t2_x[i];
        for ( i = 0; i < n1ex; i ++ )      mu3_eg[i] = mu_eg_t3_x[i];
        for ( i = 0; i < n1ex*n2ex; i ++ ) mu2_ce[i] = mu_ce_t2_x[i];
        for ( i = 0; i < n1ex*n2ex; i ++ ) mu3_ce[i] = mu_ce_t3_x[i];
    }
    else if ( k == 6 ){ // YYZZ
        for ( i = 0; i < n1ex; i ++ )      mu0_eg[i] = eiH1_t0t1_mu_eg_t0_y[it0*n1ex+i];
        for ( i = 0; i < n1ex; i ++ )      mu1_eg[i] = mu_eg_t1_y[i];
        for ( i = 0; i < n1ex; i ++ )      mu2_eg[i] = mu_eg_t2_z[i];
        for ( i = 0; i < n1ex; i ++ )      mu3_eg[i] = mu_eg_t3_z[i];
        for ( i = 0; i < n1ex*n2ex; i ++ ) mu2_ce[i] = mu_ce_t2_z[i];
        for ( i = 0; i < n1ex*n2ex; i ++ ) mu3_ce[i] = mu_ce_t3_z[i];
    }
    else if ( k == 7 ){ // ZZXX
        for ( i = 0; i < n1ex; i ++ )      mu0_eg[i] = eiH1_t0t1_mu_eg_t0_z[it0*n1ex+i];
        for ( i = 0; i < n1ex; i ++ )      mu1_eg[i] = mu_eg_t1_z[i];
        for ( i = 0; i < n1ex; i ++ )      mu2_eg[i] = mu_eg_t2_x[i];
        for ( i = 0; i < n1ex; i ++ )      mu3_eg[i] = mu_eg_t3_x[i];
        for ( i = 0; i < n1ex*n2ex; i ++ ) mu2_ce[i] = mu_ce_t2_x[i];
        for ( i = 0; i < n1ex*n2ex; i ++ ) mu3_ce[i] = mu_ce_t3_x[i];
    }
    else if ( k == 8 ){ // ZZYY
        for ( i = 0; i < n1ex; i ++ )      mu0_eg[i] = eiH1_t0t1_mu_eg_t0_z[it0*n1ex+i];
        for ( i = 0; i < n1ex; i ++ )      mu1_eg[i] = mu_eg_t1_z[i];
        for ( i = 0; i < n1ex; i ++ )      mu2_eg[i] = mu_eg_t2_y[i];
        for ( i = 0; i < n1ex; i ++ )      mu3_eg[i] = mu_eg_t3_y[i];
        for ( i = 0; i < n1ex*n2ex; i ++ ) mu2_ce[i] = mu_ce_t2_y[i];
        for ( i = 0; i < n1ex*n2ex; i ++ ) mu3_ce[i] = mu_ce_t3_y[i];
    }
    else if ( k == 9 ){ // XYXY
        for ( i = 0; i < n1ex; i ++ )      mu0_eg[i] = eiH1_t0t1_mu_eg_t0_x[it0*n1ex+i];
        for ( i = 0; i < n1ex; i ++ )      mu1_eg[i] = mu_eg_t1_y[i];
        for ( i = 0; i < n1ex; i ++ )      mu2_eg[i] = mu_eg_t2_x[i];
        for ( i = 0; i < n1ex; i ++ )      mu3_eg[i] = mu_eg_t3_y[i];
        for ( i = 0; i < n1ex*n2ex; i ++ ) mu2_ce[i] = mu_ce_t2_x[i];
        for ( i = 0; i < n1ex*n2ex; i ++ ) mu3_ce[i] = mu_ce_t3_y[i];
    }
    else if ( k == 10){ // XZXZ
        for ( i = 0; i < n1ex; i ++ )      mu0_eg[i] = eiH1_t0t1_mu_eg_t0_x[it0*n1ex+i];
        for ( i = 0; i < n1ex; i ++ )      mu1_eg[i] = mu_eg_t1_z[i];
        for ( i = 0; i < n1ex; i ++ )      mu2_eg[i] = mu_eg_t2_x[i];
        for ( i = 0; i < n1ex; i ++ )      mu3_eg[i] = mu_eg_t3_z[i];
        for ( i = 0; i < n1ex*n2ex; i ++ ) mu2_ce[i] = mu_ce_t2_x[i];
        for ( i = 0; i < n1ex*n2ex; i ++ ) mu3_ce[i] = mu_ce_t3_z[i];
    }
    else if ( k == 11){ // YXYX
        for ( i = 0; i < n1ex; i ++ )      mu0_eg[i] = eiH1_t0t1_mu_eg_t0_y[it0*n1ex+i];
        for ( i = 0; i < n1ex; i ++ )      mu1_eg[i] = mu_eg_t1_x[i];
        for ( i = 0; i < n1ex; i ++ )      mu2_eg[i] = mu_eg_t2_y[i];
        for ( i = 0; i < n1ex; i ++ )      mu3_eg[i] = mu_eg_t3_x[i];
        for ( i = 0; i < n1ex*n2ex; i ++ ) mu2_ce[i] = mu_ce_t2_y[i];
        for ( i = 0; i < n1ex*n2ex; i ++ ) mu3_ce[i] = mu_ce_t3_x[i];
    }
    else if ( k == 12){ // YZYZ
        for ( i = 0; i < n1ex; i ++ )      mu0_eg[i] = eiH1_t0t1_mu_eg_t0_y[it0*n1ex+i];
        for ( i = 0; i < n1ex; i ++ )      mu1_eg[i] = mu_eg_t1_z[i];
        for ( i = 0; i < n1ex; i ++ )      mu2_eg[i] = mu_eg_t2_y[i];
        for ( i = 0; i < n1ex; i ++ )      mu3_eg[i] = mu_eg_t3_z[i];
        for ( i = 0; i < n1ex*n2ex; i ++ ) mu2_ce[i] = mu_ce_t2_y[i];
        for ( i = 0; i < n1ex*n2ex; i ++ ) mu3_ce[i] = mu_ce_t3_z[i];
    }
    else if ( k == 13){ // ZXZX
        for ( i = 0; i < n1ex; i ++ )      mu0_eg[i] = eiH1_t0t1_mu_eg_t0_z[it0*n1ex+i];
        for ( i = 0; i < n1ex; i ++ )      mu1_eg[i] = mu_eg_t1_x[i];
        for ( i = 0; i < n1ex; i ++ )      mu2_eg[i] = mu_eg_t2_z[i];
        for ( i = 0; i < n1ex; i ++ )      mu3_eg[i] = mu_eg_t3_x[i];
        for ( i = 0; i < n1ex*n2ex; i ++ ) mu2_ce[i] = mu_ce_t2_z[i];
        for ( i = 0; i < n1ex*n2ex; i ++ ) mu3_ce[i] = mu_ce_t3_x[i];
    }
    else if ( k == 14){ // ZYZY
        for ( i = 0; i < n1ex; i ++ )      mu0_eg[i] = eiH1_t0t1_mu_eg_t0_z[it0*n1ex+i];
        for ( i = 0; i < n1ex; i ++ )      mu1_eg[i] = mu_eg_t1_y[i];
        for ( i = 0; i < n1ex; i ++ )      mu2_eg[i] = mu_eg_t2_z[i];
        for ( i = 0; i < n1ex; i ++ )      mu3_eg[i] = mu_eg_t3_y[i];
        for ( i = 0; i < n1ex*n2ex; i ++ ) mu2_ce[i] = mu_ce_t2_z[i];
        for ( i = 0; i < n1ex*n2ex; i ++ ) mu3_ce[i] = mu_ce_t3_y[i];
    }
    else if ( k == 15){ // XYYX
        for ( i = 0; i < n1ex; i ++ )      mu0_eg[i] = eiH1_t0t1_mu_eg_t0_x[it0*n1ex+i];
        for ( i = 0; i < n1ex; i ++ )      mu1_eg[i] = mu_eg_t1_y[i];
        for ( i = 0; i < n1ex; i ++ )      mu2_eg[i] = mu_eg_t2_y[i];
        for ( i = 0; i < n1ex; i ++ )      mu3_eg[i] = mu_eg_t3_x[i];
        for ( i = 0; i < n1ex*n2ex; i ++ ) mu2_ce[i] = mu_ce_t2_y[i];
        for ( i = 0; i < n1ex*n2ex; i ++ ) mu3_ce[i] = mu_ce_t3_x[i];
    }
    else if ( k == 16){ // XZZX
        for ( i = 0; i < n1ex; i ++ )      mu0_eg[i] = eiH1_t0t1_mu_eg_t0_x[it0*n1ex+i];
        for ( i = 0; i < n1ex; i ++ )      mu1_eg[i] = mu_eg_t1_z[i];
        for ( i = 0; i < n1ex; i ++ )      mu2_eg[i] = mu_eg_t2_z[i];
        for ( i = 0; i < n1ex; i ++ )      mu3_eg[i] = mu_eg_t3_x[i];
        for ( i = 0; i < n1ex*n2ex; i ++ ) mu2_ce[i] = mu_ce_t2_z[i];
        for ( i = 0; i < n1ex*n2ex; i ++ ) mu3_ce[i] = mu_ce_t3_x[i];
    }
    else if ( k == 17){ // YXXY
        for ( i = 0; i < n1ex; i ++ )      mu0_eg[i] = eiH1_t0t1_mu_eg_t0_y[it0*n1ex+i];
        for ( i = 0; i < n1ex; i ++ )      mu1_eg[i] = mu_eg_t1_x[i];
        for ( i = 0; i < n1ex; i ++ )      mu2_eg[i] = mu_eg_t2_x[i];
        for ( i = 0; i < n1ex; i ++ )      mu3_eg[i] = mu_eg_t3_y[i];
        for ( i = 0; i < n1ex*n2ex; i ++ ) mu2_ce[i] = mu_ce_t2_x[i];
        for ( i = 0; i < n1ex*n2ex; i ++ ) mu3_ce[i] = mu_ce_t3_y[i];
    }
    else if ( k == 18){ // YZZY
        for ( i = 0; i < n1ex; i ++ )      mu0_eg[i] = eiH1_t0t1_mu_eg_t0_y[it0*n1ex+i];
        for ( i = 0; i < n1ex; i ++ )      mu1_eg[i] = mu_eg_t1_z[i];
        for ( i = 0; i < n1ex; i ++ )      mu2_eg[i] = mu_eg_t2_z[i];
        for ( i = 0; i < n1ex; i ++ )      mu3_eg[i] = mu_eg_t3_y[i];
        for ( i = 0; i < n1ex*n2ex; i ++ ) mu2_ce[i] = mu_ce_t2_z[i];
        for ( i = 0; i < n1ex*n2ex; i ++ ) mu3_ce[i] = mu_ce_t3_y[i];
    }
    else if ( k == 19){ // ZXXZ
        for ( i = 0; i < n1ex; i ++ )      mu0_eg[i] = eiH1_t0t1_mu_eg_t0_z[it0*n1ex+i];
        for ( i = 0; i < n1ex; i ++ )      mu1_eg[i] = mu_eg_t1_x[i];
        for ( i = 0; i < n1ex; i ++ )      mu2_eg[i] = mu_eg_t2_x[i];
        for ( i = 0; i < n1ex; i ++ )      mu3_eg[i] = mu_eg_t3_z[i];
        for ( i = 0; i < n1ex*n2ex; i ++ ) mu2_ce[i] = mu_ce_t2_x[i];
        for ( i = 0; i < n1ex*n2ex; i ++ ) mu3_ce[i] = mu_ce_t3_z[i];
    }
    else if ( k == 20){ // ZYYZ
        for ( i = 0; i < n1ex; i ++ )      mu0_eg[i] = eiH1_t0t1_mu_eg_t0_z[it0*n1ex+i];
        for ( i = 0; i < n1ex; i ++ )      mu1_eg[i] = mu_eg_t1_y[i];
        for ( i = 0; i < n1ex; i ++ )      mu2_eg[i] = mu_eg_t2_y[i];
        for ( i = 0; i < n1ex; i ++ )      mu3_eg[i] = mu_eg_t3_z[i];
        for ( i = 0; i < n1ex*n2ex; i ++ ) mu2_ce[i] = mu_ce_t2_y[i];
        for ( i = 0; i < n1ex*n2ex; i ++ ) mu3_ce[i] = mu_ce_t3_z[i];
    }
    else{
        cout << "ERROR:: getMU4pol k= " << k << " unknown." << endl;
        return 1;
    }

    return IR2DOK;

}

complex<double> IR2D::getR2D_R1(int it0)
// return the third order rephasing response function
{
    complex<double> *mu0_eg, *mu1_eg, *mu2_eg, *mu3_eg, *work1a, *work1b, *work1c;
    complex<double> *mu2_ce, *mu3_ce, *work2a, *work2b;
    complex<double> R2D, work0a, work0b;
    int i, k;
    double weight;

    mu0_eg = new complex<double>[n1ex];
    mu1_eg = new complex<double>[n1ex];
    mu2_eg = new complex<double>[n1ex];
    mu3_eg = new complex<double>[n1ex];
    mu2_ce = new complex<double>[n1ex*n2ex];
    mu3_ce = new complex<double>[n1ex*n2ex];
    work1a = new complex<double>[n1ex];
    work1b = new complex<double>[n1ex];
    work1c = new complex<double>[n1ex];
    work2a = new complex<double>[n2ex];
    work2b = new complex<double>[n2ex];

    R2D = complex_zero;

    // get the isotropically averaged spectrum -- see notes sec 5
    for ( k = 0; k < 21; k ++ ){
        if ( getMU4pol( mu0_eg, mu1_eg, mu2_eg, mu3_eg, mu2_ce, mu3_ce, it0, k ) != IR2DOK ){
            exit(EXIT_FAILURE);
        }
        // assign weights for parallel isotropic spectrum -- see notes sec 5
        if ( k <=2 ) weight = 1./5.;
        else weight = 1./15;

        // ground state bleach -- eq (27)
        // ############################################################

        // start matrix algebra from right
        // work1a=eiH1_t2t3*mu2_eg
        cblas_zgemv( CblasRowMajor, CblasNoTrans, n1ex, n1ex,\
                     &complex_one, eiH1_t2t3, n1ex, mu2_eg, 1, \
                     &complex_zero, work1a, 1 );
        // work0a=mu3_eg*work1a
        cblas_zdotu_sub( n1ex, mu3_eg, 1, work1a, 1, &work0a );

        // now do matrix algebra from left
        // work0b=conj(mu0_eg)*mu1_eg
        for ( int i = 0; i < n1ex; i ++ ) work1a[i] = conj(mu0_eg[i]);
        cblas_zdotu_sub( n1ex, work1a, 1, mu1_eg, 1, &work0b );

        // now meet in the middle to get the result
        R2D   += weight*work0b*work0a;
        // ############################################################


        // stimulated emission -- eq (28)
        // ############################################################

        // start matrix algebra from right
        // work1a=eiH1_t1t2*mu1_eg
        cblas_zgemv( CblasRowMajor, CblasNoTrans, n1ex, n1ex,\
                     &complex_one, eiH1_t1t2, n1ex, mu1_eg, 1, \
                     &complex_zero, work1a, 1 );
        // work1b=eiH1_t2t3*work1a
        cblas_zgemv( CblasRowMajor, CblasNoTrans, n1ex, n1ex,\
                     &complex_one, eiH1_t2t3, n1ex, work1a, 1, \
                     &complex_zero, work1b, 1 );
        // work0a=mu3_eg*work1b
        cblas_zdotu_sub( n1ex, mu3_eg, 1, work1b, 1, &work0a );

        // now do matrix algebra from left
        // work1a=eiH1_t1t2*mu0_eg
        cblas_zgemv( CblasRowMajor, CblasNoTrans, n1ex, n1ex,\
                     &complex_one, eiH1_t1t2, n1ex, mu0_eg, 1, \
                     &complex_zero, work1a, 1 );
        // work0b=conj(work1a)*mu2_eg
        for ( int i = 0; i < n1ex; i ++ ) work1b[i] = conj(work1a[i]);
        cblas_zdotu_sub( n1ex, work1b, 1, mu2_eg, 1, &work0b );
    
        // now meet in the middle to get the result
        R2D   += weight*work0b*work0a;
        // ############################################################


        // excited state absorption -- eq (29)
        // ############################################################

        // start matrix algebra from right
        // work1a=eiH1_t1t2*mu1_eg
        cblas_zgemv( CblasRowMajor, CblasNoTrans, n1ex, n1ex,\
                     &complex_one, eiH1_t1t2, n1ex, mu1_eg, 1, \
                     &complex_zero, work1a, 1 );
        // work2a=mu2_ce'*work1a
        cblas_zgemv( CblasRowMajor, CblasTrans, n1ex, n2ex,\
                     &complex_one, mu2_ce, n2ex, work1a, 1, \
                     &complex_zero, work2a, 1 );
        // work2b=eiH2_t2t3*work2a
        cblas_zgemv( CblasRowMajor, CblasNoTrans, n2ex, n2ex,\
                     &complex_one, eiH2_t2t3, n2ex, work2a, 1, \
                     &complex_zero, work2b, 1 );
        // work1c=mu3_ce*work2b
        cblas_zgemv( CblasRowMajor, CblasNoTrans, n1ex, n2ex,\
                     &complex_one, mu3_ce, n2ex, work2b, 1, \
                     &complex_zero, work1c, 1 );

        // now do matrix algebra from the left
        // work1a=eiH1_t1t2*mu0_eg
        cblas_zgemv( CblasRowMajor, CblasNoTrans, n1ex, n1ex,\
                     &complex_one, eiH1_t1t2, n1ex, mu0_eg, 1, \
                     &complex_zero, work1a, 1 );
        // work1b=eiH1_t2t3*work1a
        cblas_zgemv( CblasRowMajor, CblasNoTrans, n1ex, n1ex,\
                     &complex_one, eiH1_t2t3, n1ex, work1a, 1, \
                     &complex_zero, work1b, 1 );
        
        // now meet in the middle and combine the result
        // work0a = conj(work1b)*work1c
        for ( i = 0; i < n1ex; i ++ ) work1a[i] = conj(work1b[i]);
        cblas_zdotu_sub( n1ex, work1a, 1, work1c, 1, &work0a );
        R2D    -= weight*work0a;
        // ############################################################
    }

    delete [] mu0_eg;
    delete [] mu1_eg;
    delete [] mu2_eg;
    delete [] mu3_eg;
    delete [] mu2_ce;
    delete [] mu3_ce;
    delete [] work1a;
    delete [] work1b;
    delete [] work1c;
    delete [] work2a;
    delete [] work2b;

    return R2D;
}

complex<double> IR2D::getR2D_R2(int it0)
// return the third order rephasing response function
{
    complex<double> *mu0_eg, *mu1_eg, *mu2_eg, *mu3_eg, *work1a, *work1b, *work1c;
    complex<double> *mu2_ce, *mu3_ce, *work2a, *work2b;
    complex<double> R2D, work0a, work0b;
    int i, k;
    double weight;

    mu0_eg = new complex<double>[n1ex];
    mu1_eg = new complex<double>[n1ex];
    mu2_eg = new complex<double>[n1ex];
    mu3_eg = new complex<double>[n1ex];
    mu2_ce = new complex<double>[n1ex*n2ex];
    mu3_ce = new complex<double>[n1ex*n2ex];
    work1a = new complex<double>[n1ex];
    work1b = new complex<double>[n1ex];
    work1c = new complex<double>[n1ex];
    work2a = new complex<double>[n2ex];
    work2b = new complex<double>[n2ex];

    R2D = complex_zero;

    // get the isotropically averaged spectrum -- see notes sec 5
    for ( k = 0; k < 21; k ++ ){
        if ( getMU4pol( mu0_eg, mu1_eg, mu2_eg, mu3_eg, mu2_ce, mu3_ce, it0, k ) != IR2DOK ){
            exit(EXIT_FAILURE);
        }
        // assign weights for parallel isotropic spectrum -- see notes sec 5
        if ( k <=2 ) weight = 1./5.;
        else weight = 1./15;

        // ground state bleach -- eq (30)
        // ############################################################

        // start matrix algebra from left
        // work1a=eiH1_t2t3*mu2_eg
        cblas_zgemv( CblasRowMajor, CblasNoTrans, n1ex, n1ex,\
                     &complex_one, eiH1_t2t3, n1ex, mu2_eg, 1, \
                     &complex_zero, work1a, 1 );
        // work0a=mu3_eg*work1a
        cblas_zdotu_sub( n1ex, mu3_eg, 1, work1a, 1, &work0a );

        // now do matrix algebra from right
        // work0b=mu0_eg*mu1_eg
        cblas_zdotu_sub( n1ex, mu0_eg, 1, mu1_eg, 1, &work0b );

        // now meet in the middle to get the result
        R2D   += weight*work0b*work0a;
        // ############################################################


        // stimulated emission -- eq (31)
        // ############################################################
 
        // now do matrix algebra from right
        // work1a=eiH1_t1t2*mu0_eg
        cblas_zgemv( CblasRowMajor, CblasNoTrans, n1ex, n1ex,\
                     &complex_one, eiH1_t1t2, n1ex, mu0_eg, 1, \
                     &complex_zero, work1a, 1 );
        // work1b=eiH1_t2t3*work1a
        cblas_zgemv( CblasRowMajor, CblasNoTrans, n1ex, n1ex,\
                     &complex_one, eiH1_t2t3, n1ex, work1a, 1, \
                     &complex_zero, work1b, 1 );
        // work0b=work1b*mu3_eg
        cblas_zdotu_sub( n1ex, work1b, 1, mu3_eg, 1, &work0b );


        // now do matrix algebra from left
        // work1a=eiH1_t1t2*mu1_eg
        cblas_zgemv( CblasRowMajor, CblasNoTrans, n1ex, n1ex,\
                     &complex_one, eiH1_t1t2, n1ex, mu1_eg, 1, \
                     &complex_zero, work1a, 1 );
        // work0a=conj(work1a)*mu2_eg
        for ( int i = 0; i < n1ex; i ++ ) work1b[i] = conj(work1a[i]);
        cblas_zdotu_sub( n1ex, work1b, 1, mu2_eg, 1, &work0a );
        R2D   += weight*work0b*work0a;
        // ############################################################


        // excited state absorption -- eq (32)
        // ############################################################

        // start matrix algebra from right
        // work1a=eiH1_t1t2*mu0_eg
        cblas_zgemv( CblasRowMajor, CblasNoTrans, n1ex, n1ex,\
                     &complex_one, eiH1_t1t2, n1ex, mu0_eg, 1, \
                     &complex_zero, work1a, 1 );
        // work2a=mu2_ce*work1a
        cblas_zgemv( CblasRowMajor, CblasTrans, n1ex, n2ex,\
                     &complex_one, mu2_ce, n2ex, work1a, 1, \
                     &complex_zero, work2a, 1 );
        // work2b=eiH2_t2t3*work2a
        cblas_zgemv( CblasRowMajor, CblasNoTrans, n2ex, n2ex,\
                     &complex_one, eiH2_t2t3, n2ex, work2a, 1, \
                     &complex_zero, work2b, 1 );
        // work1c=mu3_ce*work2b
        cblas_zgemv( CblasRowMajor, CblasNoTrans, n1ex, n2ex,\
                     &complex_one, mu3_ce, n2ex, work2b, 1, \
                     &complex_zero, work1c, 1 );

        // now do matrix algebra from the left
        // work1a=eiH1_t1t2*mu1_eg
        cblas_zgemv( CblasRowMajor, CblasNoTrans, n1ex, n1ex,\
                     &complex_one, eiH1_t1t2, n1ex, mu1_eg, 1, \
                     &complex_zero, work1a, 1 );
        // work1b=eiH1_t2t3*work1a
        cblas_zgemv( CblasRowMajor, CblasNoTrans, n1ex, n1ex,\
                     &complex_one, eiH1_t2t3, n1ex, work1a, 1, \
                     &complex_zero, work1b, 1 );

        // now meet in the middle and combine the results
        // work0a = conj(work1b)*work1c
        for ( i = 0; i < n1ex; i ++ ) work1a[i] = conj(work1b[i]);
        cblas_zdotu_sub( n1ex, work1a, 1, work1c, 1, &work0a );
        R2D    -= weight*work0a;
        // ############################################################
    }

    delete [] mu0_eg;
    delete [] mu1_eg;
    delete [] mu2_eg;
    delete [] mu3_eg;
    delete [] mu2_ce;
    delete [] mu3_ce;
    delete [] work1a;
    delete [] work1b;
    delete [] work1c;
    delete [] work2a;
    delete [] work2b;

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

    // write rephasing response function
    fn=_ofile_+"-RparI.dat";
    ofile.open( fn );
    if ( ! ofile.is_open() ) { fileOpenErr( fn ); return 1;}
    cout << ">>> Writing " << fn << "." << endl;

    ofile << "# Rephasing parallel ZZZZ polarized response function, t2 = " << t2 << endl;
    ofile << "# t1 (ps) t3 (ps) Real Imag" << endl;
    for ( t1 = 0; t1 < t1t3_npoints; t1 ++ ){
        for ( t3 = 0; t3 < t1t3_npoints; t3 ++ ){
            ofile << t1 * dt << " " << t3 * dt << " " \
                  << R2D_R1[t1*t1t3_npoints + t3].real() << " " \
                  << R2D_R1[t1*t1t3_npoints + t3].imag() << endl;
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
            ofile << t1 * dt << " " << t3 * dt << " " \
                << R2D_R2[t1*t1t3_npoints + t3].real() << " " \
                << R2D_R2[t1*t1t3_npoints + t3].imag() << endl;
        }
    }
    ofile.close();

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
    int      it1, i;

    // allocate arrays
    fftIn  = new complex<double>[fftlen]();
    fftOut = new complex<double>[fftlen]();
    res    = new complex<double>[fftlen]();

    // scale fftOut to keep total area of the spectrum conserved and normalize
    // by root(fftlen) since FFTW3 does not
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
    for ( it1 = 0; it1 < t1t3_npoints; it1 ++ ){
        fftIn[it1]        = R1D[it1];
        if ( it1 == 0 ) continue; // dont take hermitian at t=0
        fftIn[fftlen-it1] = conj(R1D[it1]);
    }
    fftw_execute(plan); 
    for ( it1 = 0; it1 < fftlen; it1 ++ ) res[it1].real(fftOut[it1].real());

    // Dispersive part
    for ( i = 0; i < fftlen ; i ++ ) fftIn[i] = complex_zero;
    for ( it1 = 0; it1 < t1t3_npoints; it1 ++ ){
        fftIn[it1]        = img*R1D[it1];
        if ( it1 == 0 ) continue; // dont take hermitian at t=0
        fftIn[fftlen-it1] = conj(img*R1D[it1]);
    }
    fftw_execute(plan); 
    fftw_destroy_plan(plan);
    for ( it1 = 0; it1 < fftlen; it1 ++ ) res[it1].imag(fftOut[it1].real());

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
        ofile << freq << " " << scale*res[i].real() << " " \
              << scale*res[i].imag() << endl;
    }
    // positive frequencies are stored at the beginning of fftOut
    for ( i = 0; i < fftlen/2; i ++ ){
        // get frequency in wavenumber, add back the shift
        freq   = 2.*PI*HBAR*i/(dt*fftlen) + shift;
        if ( freq < window0 or freq > window1 ) continue;
        ofile << freq << " " << scale*res[i].real() << " " \
              << scale*res[i].imag() << endl;
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
    int       it1, it3, i, j;


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
    for ( i = 0; i < fftlen*fftlen; i ++ ) fftIn[i] = complex_zero;
    for ( it1 = 0; it1 < t1t3_npoints; it1 ++ ){
        for ( it3 = 0; it3 < t1t3_npoints; it3 ++ ){
            fftIn[ it1*fftlen + it3 ] = R2D_R1[ it1*t1t3_npoints + it3 ];
            // divide t=0 point by 2 (see Hamm and Zanni sec 9.5.3)
            if ( it1 == 0 and it3 == 0 ) fftIn[ it1*fftlen + it3 ] /=2.;
        }
    }
    fftw_execute(plan);
    fn=_ofile_+"-RparIw.dat";
    if ( write2Dout( fftOut, fn, "rephasing", fftlen ) != IR2DOK ) return 1;

    // Save rephasing contribution to purly absorptive spectrum
    // here we map w1 to -w1, see Hamm and Zanni eq 4.36
    for ( i = 0; i < fftlen; i ++ ){
        for ( j = 0; j < fftlen; j ++ ){
            if ( i == 0 ) res[ i*fftlen + j ] = fftOut[ i*fftlen + j ]; // w1=0 goes to w1=0
            else          res[ i*fftlen + j ] = fftOut[ (fftlen - i)*fftlen + j ]; // w1 to -w1
        }
    }

    // fourier transform non-rephasing response functions, see Hamm and Zanni eq 4.31
    for ( i = 0; i < fftlen*fftlen; i ++ ) fftIn[i] = complex_zero;
    for ( it1 = 0; it1 < t1t3_npoints; it1 ++ ){
        for ( it3 = 0; it3 < t1t3_npoints; it3 ++ ){
            fftIn[ it1*fftlen + it3 ] = R2D_R2[ it1*t1t3_npoints + it3 ];
            // divide t=0 point by 2 (see Hamm and Zanni sec 9.5.3)
            if ( it1 == 0 and it3 == 0 ) fftIn[ it1*fftlen + it3 ] /=2.;
        }
    }
    fftw_execute(plan);
    fftw_destroy_plan(plan);
    fn=_ofile_+"-RparIIw.dat";
    if ( write2Dout( fftOut, fn, "non-rephasing", fftlen ) != IR2DOK ) return 1;
    
    // Save rephasing contribution to purly absorptive spectrum
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
    scale*= -1.*scale;  // since 2D scaling, multiply by -1 so ground state bleach 
                        // and stimulated emission are negative and excited state 
                        // absorption is positive

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
    int it1_max, it2_max, it3_max, it0_min;
    int ndx;

    // check user input to program
    if ( argc != 2 ){
        cout << "ERROR:: Program expects the name of the input" << 
                "file as the only argument. Aborting." << endl;
        exit( EXIT_FAILURE );
    }

    // get input file name and initialize IR2D class
    IR2D spectrum( argv[1] ); 

    // TODO: If you really want to increase performance, you can use sparse matrices and GPU

    // Loop over the samples
    for ( sample = 0; sample < spectrum.nsamples; sample ++ ){

        // set bounds on time points
        it0_min = sample*static_cast<int>(round(spectrum.sample_every/spectrum.dt));
        it1_max = it0_min + spectrum.t1t3_npoints - 1;
        it2_max = it1_max + static_cast<int>(round(spectrum.t2/spectrum.dt));
        it3_max = it2_max + spectrum.t1t3_npoints - 1;
        //cout << it0_min << "-" << it1_max << " " << it2_max << "-" << it3_max << endl;

        fprintf(stderr, "    Now processing sample %d/%d with t1 at %.2f ps\n", \
                sample+1, spectrum.nsamples, it1_max*spectrum.dt ); fflush(stderr);

        // do it0 to it1
        // pay attention to order -- propigate H1 to next t1 after propigating mu0
        spectrum.reset_eiH1("t0-t1");
        for ( it0 = it1_max; it0 >= it0_min; it0-- ){
            if ( spectrum.readDframe(it0) != IR2DOK ) exit(EXIT_FAILURE);
            if ( spectrum.setMUatT("t0")  != IR2DOK ) exit(EXIT_FAILURE);
            if ( spectrum.save_eiH1_t0t1_mu0(it1_max-it0) != IR2DOK ) exit(EXIT_FAILURE);
            if ( it0 == it0_min ) continue; // don't propigate the last frame
            if ( spectrum.prop_eiH1(it0,"t0-t1") != IR2DOK ) exit(EXIT_FAILURE);
        }

        // read in dipoles at t1
        if ( spectrum.readDframe(it1_max) != IR2DOK ) exit(EXIT_FAILURE);
        if ( spectrum.setMUatT("t1")      != IR2DOK ) exit(EXIT_FAILURE);

        // get linear response function
        for ( it0 = 0; it0 < spectrum.t1t3_npoints; it0 ++ ){
            spectrum.R1D[it0] += spectrum.getR1D(it0);
        }
        
        // get eiH_t1t2
        spectrum.reset_eiH1("t1-t2");
        spectrum.reset_eiH2("t1-t2");
        for ( it2 = it1_max; it2 < it2_max; it2 ++ ){
            if ( spectrum.prop_eiH1(it2,"t1-t2" ) != IR2DOK ) exit(EXIT_FAILURE);
            if ( spectrum.prop_eiH2(it2,"t1-t2" ) != IR2DOK ) exit(EXIT_FAILURE);
        }

        // read in dipoles at t2
        if ( spectrum.readDframe(it2_max) != IR2DOK ) exit(EXIT_FAILURE);
        if ( spectrum.setMUatT("t2")      != IR2DOK ) exit(EXIT_FAILURE);

        // now loop over t3 and get the response function
        spectrum.reset_eiH1("t2-t3");
        spectrum.reset_eiH2("t2-t3");
        for ( it3 = it2_max; it3 <= it3_max; it3 ++ ){
            if ( spectrum.readDframe(it3) != IR2DOK ) exit(EXIT_FAILURE);
            if ( spectrum.setMUatT("t3")  != IR2DOK ) exit(EXIT_FAILURE);
            // get response function for all it0
            #pragma omp parallel for private(ndx)
            for ( it0 = 0; it0 < spectrum.t1t3_npoints; it0 ++ ){
                ndx = it0*spectrum.t1t3_npoints + it3-it2_max;
                spectrum.R2D_R1[ndx] += spectrum.getR2D_R1(it0);
                spectrum.R2D_R2[ndx] += spectrum.getR2D_R2(it0);
            }
            printProgress( it3-it2_max+1, spectrum.t1t3_npoints );
            if ( it3 == it3_max ) continue; // dont propigate the last frame
            if ( spectrum.prop_eiH1(it3,"t2-t3" ) != IR2DOK ) exit(EXIT_FAILURE);
            if ( spectrum.prop_eiH2(it3,"t2-t3" ) != IR2DOK ) exit(EXIT_FAILURE);
        }
        cerr << endl;
    }

    // account for dephasing phenomonelogically and normalize
    for ( it1 = 0; it1 < spectrum.t1t3_npoints; it1 ++ ){
        // dephasing and normalization -- linear
        spectrum.R1D[it1]*=exp(-it1*spectrum.dt/spectrum.lifetime_T2)/(1.*spectrum.nsamples);
        for ( it3 = 0; it3 < spectrum.t1t3_npoints; it3 ++ ){
            ndx = it1 * spectrum.t1t3_npoints + it3;
            // dephasing and normalization -- rephasing
            // note that the relaxation from the excited state absorption pathways
            // should have a different relaxation time, but we ignore that here...
            spectrum.R2D_R1[ ndx ]*= exp(-(it1+2.*it2+it3)*spectrum.dt/\
                                     spectrum.lifetime_T2) / ( 1.*spectrum.nsamples );
            // dephasing and normalization -- non-rephasing
            spectrum.R2D_R2[ ndx ]*= exp(-(it1+2.*it2+it3)*spectrum.dt/\
                                     spectrum.lifetime_T2) / ( 1.*spectrum.nsamples );
        }
    }

    // do fourier transforms and write them out
    if ( spectrum.writeR1D()    != IR2DOK ) exit(EXIT_FAILURE);
    if ( spectrum.writeR2D()    != IR2DOK ) exit(EXIT_FAILURE);
    if ( spectrum.write1Dfft()  != IR2DOK ) exit(EXIT_FAILURE);
    if ( spectrum.write2DRabs() != IR2DOK ) exit(EXIT_FAILURE);
    cout << ">>> Done!" << endl;
}

// ===========================================================================
//  lattice_stub.hpp — NOT a lattice implementation.
//
//  The C++ sources in this repository are meant to be #included by the
//  production driver, which supplies a full parallel lattice library (gauge
//  fields, site arithmetic, MPI, OpenMP).  That library is external and is not
//  distributed here.
//
//  This header declares just enough of that interface for the sources to
//  COMPILE, so that the parts which do not touch gauge fields at all -- the
//  trajectory symmetry orbits and the trie construction -- can be built and
//  exercised on their own.  That is what src/hpe_cost.cpp does when it counts
//  matrix multiplications.
//
//  Every operation below is a placeholder.  Nothing here computes anything:
//  linking against it and evaluating a Wilson loop yields garbage by design.
//  Use it only for cost counting and compile checking.
// ===========================================================================
#ifndef HPE_LATTICE_STUB_HPP
#define HPE_LATTICE_STUB_HPP

#include <complex.h>
#include <vector>
#include <array>
#include <iostream>
#include <cstdlib>

using std::vector;
using std::endl;

template<int NR, int NC>
class Matrix
{
public:
  double _Complex* address;
  double _Complex  own[NR*NC];
  Matrix() : address(own) {}
  explicit Matrix( double _Complex* p ) : address(p) {}
  Matrix( const Matrix& o ) : address(own) {
    for( int i=0 ; i<NR*NC ; i++ ) own[i] = o.address[i];
  }
  Matrix& operator=( const Matrix& o ){
    for( int i=0 ; i<NR*NC ; i++ ) address[i] = o.address[i];
    return *this;
  }
  Matrix& operator=( double _Complex s ){
    for( int i=0 ; i<NR*NC ; i++ ) address[i] = 0.;
    for( int i=0 ; i<NR && i<NC ; i++ ) address[i*NC+i] = s;
    return *this;
  }
  Matrix operator*( const Matrix& ) const { return Matrix(); }
  Matrix operator+( const Matrix& ) const { return Matrix(); }
  Matrix& operator+=( const Matrix& ){ return *this; }
};

template<int NR, int NC> Matrix<NR,NC> Hermite  ( const Matrix<NR,NC>& ){ return Matrix<NR,NC>(); }
template<int NR, int NC> double _Complex Trace  ( const Matrix<NR,NC>& ){ return 0.; }
template<int NR, int NC> double _Complex TraceProd( const Matrix<NR,NC>& , const Matrix<NR,NC>& ){ return 0.; }

struct StubLattice { long size() const { return 1; } };

struct Site
{
  StubLattice* lat = nullptr;
  Site() {}
  explicit Site( const StubLattice& l ) : lat( const_cast<StubLattice*>(&l) ) {}
  Site operator+( int ) const { return *this; }
  Site operator-( int ) const { return *this; }
};

struct Gauge
{
  StubLattice lat;
  const StubLattice& lattice() const { return lat; }
  Matrix<3,3> operator()( const Site& , int ) const { return Matrix<3,3>(); }
};

struct StubMPI {
  void add( double& ) {}
  void add( double _Complex& ) {}
  void abort() { std::exit(1); }
};
static StubMPI mpi;

#define mpicout std::cout
#define ForAllSitesOMP(x) for( int _stub_once = 0 ; _stub_once < 1 ; _stub_once++ )

#endif

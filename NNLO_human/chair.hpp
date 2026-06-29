template< int mu , int nu , int rho >
inline
double Chair_local( const Gauge& V , const Site& x )
{
  double _Complex ptr1[9] , ptr2[9] , ptr3[9];
  Const_Matrix<3,3,false,false,1> tmp1( ptr1 );
  Const_Matrix<3,3,false,false,1> tmp2( ptr2 );
  Const_Matrix<3,3,false,false,1> tmp3( ptr3 );

  Site y(V.lattice()) , z(V.lattice());
  Site xm = x+mu;

  y = x+nu;
  z = x-nu;

  tmp1 = V(x,nu) * V(y,mu);
  tmp2 = tmp1 * Hermite( V(xm,nu) );
  tmp1 = Hermite( V(z,nu) ) * V(z,mu);
  tmp2 += tmp1 * V(z+mu,nu);
  
  y = x+rho;
  z = x-rho;

  tmp1 = V(x,rho) * V(y,mu);
  tmp3 = tmp1 * Hermite( V(xm,rho) );
  tmp1 = Hermite( V(z,rho) ) * V(z,mu);
  tmp3 += tmp1 * V(z+mu,rho);
  
  tmp1 = tmp2 * Hermite( tmp3 );

  return __real__( Trace(tmp1) ) / 3. / 4.;
}

double Chair_avS( const Gauge& V )
{
  double ans = 0.;
#pragma omp parallel reduction( + : ans )
  {
    Site x(V.lattice());
    ForAllSitesOMP(x){
      ans += Chair_local<1,2,3>( V , x );
      ans += Chair_local<2,3,1>( V , x );
      ans += Chair_local<3,1,2>( V , x );
    }
  } // end of omp parallel
  mpi.add(ans);
  return ans / V.lattice().size() / 3.;
}

double Chair_avT1( const Gauge& V )
{
  double ans = 0.;
#pragma omp parallel reduction( + : ans )
  {
    Site x(V.lattice());
    double _Complex ptr1[9] , ptr2[9] , ptr3[9];
    ForAllSitesOMP(x){
      ans += Chair_local<1,0,2>( V , x );
      ans += Chair_local<1,0,3>( V , x );
      ans += Chair_local<2,0,1>( V , x );
      ans += Chair_local<2,0,3>( V , x );
      ans += Chair_local<3,0,1>( V , x );
      ans += Chair_local<3,0,2>( V , x );
    }
  } // end of omp parallel
  mpi.add(ans);
  return ans / V.lattice().size() / 6.;
}
double Chair_avT2( const Gauge& V )
{
  double ans = 0.;
#pragma omp parallel reduction( + : ans )
  {
    Site x(V.lattice());
    ForAllSitesOMP(x){
      ans += Chair_local<0,1,2>( V , x );
      ans += Chair_local<0,2,3>( V , x );
      ans += Chair_local<0,3,1>( V , x );
    }
  } // end of omp parallel
  mpi.add(ans);
  return ans / V.lattice().size() / 3.;
}



template< int mu , int nu , int rho >
inline
double Crown_local( const Gauge& V , const Site& x )
{
  double _Complex ptr1[9] , ptr2[9];
  Const_Matrix<3,3,false,false,1> tmp1( ptr1 );
  Const_Matrix<3,3,false,false,1> tmp2( ptr2 );

  Site y(V.lattice()) , z(V.lattice());

  y = x+mu;
  z = x+rho;
  tmp1 = Hermite( V(x,rho) ) * V(x,mu);
  tmp1 = tmp1 * V(y,nu);
  tmp1 = tmp1 * V(y+nu,rho);
  tmp1 = tmp1 * Hermite( V(z+nu,mu) );
  tmp2 = tmp1 * Hermite( V(z,nu) );

  y = x+nu;
  z = x+rho;
  tmp1 = Hermite( V(x,rho) ) * V(x,nu);
  tmp1 = tmp1 * V(y,mu);
  tmp1 = tmp1 * V(y+mu,rho);
  tmp1 = tmp1 * Hermite( V(z+mu,nu) );
  tmp2 += tmp1 * Hermite( V(z,mu) );

  y = x+mu-rho;
  z = x+nu-rho;
  tmp1 = Hermite( V(x,nu) )  * V(x,mu);
  tmp1 = tmp1 * Hermite( V(y,rho) );
  tmp1 = tmp1 * V(y,nu);
  tmp1 = tmp1 * Hermite( V(z,mu) );
  tmp2 += tmp1 * V(z,rho);

  y = x+mu;
  z = x+nu;
  tmp1 = Hermite( V(x,nu) ) * V(x,mu);
  tmp1 = tmp1 * V(y,rho);
  tmp1 = tmp1 * V(y+rho,nu);
  tmp1 = tmp1 * Hermite( V(z+rho,mu) );
  tmp2 += tmp1 * Hermite( V(z,rho) );

  return __real__( Trace(tmp2) ) / 3. / 4.;
}
double Crown_avS( const Gauge& V )
{
  double ans = 0.;
#pragma omp parallel reduction( + : ans )
  {
    Site x(V.lattice());
    ForAllSitesOMP(x){
      ans += Crown_local<1,2,3>( V , x );
    }
  } // end of omp parallel
  mpi.add(ans);
  return ans / V.lattice().size();
}

double Crown_avT( const Gauge& V )
{
  double ans = 0.;
#pragma omp parallel reduction( + : ans )
  {
    Site x(V.lattice());
    ForAllSitesOMP(x){
      ans += Crown_local<0,1,2>( V , x );
      ans += Crown_local<0,1,3>( V , x );
      ans += Crown_local<0,2,3>( V , x );
    }
  } // end of omp parallel
  mpi.add(ans);
  return ans / V.lattice().size() / 3.;
}

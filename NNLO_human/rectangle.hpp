
template< int mu , int nu >
inline
double Rectangle_local( const Gauge& V , const Site& x , 
			double _Complex *plaq_ptr )
{
  Const_Matrix<3,3,false,false,1> rect( plaq_ptr );
  Site y = x - nu;
  rect = Hermite( V( x,nu ) ) * Hermite( V( y,nu ) );
  Site z = y + mu;
  rect = rect * V( y,mu );
  y = x + mu;
  rect = rect * V( z,nu );
  z = x + nu;
  rect = rect * V( y,nu );
  rect = rect * Hermite( V( z,mu ) );
  return __real__( Trace(rect) ) / 3.;
}

double Rectangle_avS( const Gauge& V )
{
  double ans = 0.;
#pragma omp parallel reduction( + : ans )
  {
    Site x(V.lattice());
    ForAllSitesOMP(x){
      double _Complex plaq_ptr[9];
      ans += Rectangle_local<1,2>( V , x , plaq_ptr );
      ans += Rectangle_local<2,1>( V , x , plaq_ptr );
      ans += Rectangle_local<2,3>( V , x , plaq_ptr );
      ans += Rectangle_local<3,2>( V , x , plaq_ptr );
      ans += Rectangle_local<3,1>( V , x , plaq_ptr );
      ans += Rectangle_local<1,3>( V , x , plaq_ptr );
    }
  } // end of omp parallel
  mpi.add(ans);
  return ans / V.lattice().size() / 6.;
}

double Rectangle_avT1( const Gauge& V )
{
  double ans = 0.;
#pragma omp parallel reduction( + : ans )
  {
    Site x(V.lattice());
    ForAllSitesOMP(x){
      double _Complex plaq_ptr[9];
      ans += Rectangle_local<0,1>( V , x , plaq_ptr );
      ans += Rectangle_local<0,2>( V , x , plaq_ptr );
      ans += Rectangle_local<0,3>( V , x , plaq_ptr );
    }
  } // end of omp parallel
  mpi.add(ans);
  return ans / V.lattice().size() / 3.;
}

double Rectangle_avT2( const Gauge& V )
{
  double ans = 0.;
#pragma omp parallel reduction( + : ans )
  {
    Site x(V.lattice());
    ForAllSitesOMP(x){
      double _Complex plaq_ptr[9];
      ans += Rectangle_local<1,0>( V , x , plaq_ptr );
      ans += Rectangle_local<2,0>( V , x , plaq_ptr );
      ans += Rectangle_local<3,0>( V , x , plaq_ptr );
    }
  } // end of omp parallel
  mpi.add(ans);
  return ans / V.lattice().size() / 3.;
}


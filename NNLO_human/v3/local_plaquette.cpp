class LocalPlaquette {

public:
  mutable double _Complex m[9*4*12];

  LocalPlaquette( const Gauge& V , const SingleStaple& S ,
		  const Site& x ){
    Calc( V , S , x );
  }
  //        nu
  //        ^
  //        |
  // +->-+  |  +-<-+
  // |   |  |  |   |
  // +-<-+  |  +->-+
  // -------+--------> mu
  // +-<-+  |  +->-+
  // |   |  |  |   |
  // +->-+  |  +-<-+
  //        |
  int index( int mu , int nu , int sign_mu , int sign_nu ) const {
    static constexpr int idx[4][4] = {
      {-1,0,1,2} ,
      {3,-1,4,5} ,
      {6,7,-1,8} ,
      {9,10,11,-1} };
    const int local = sign_mu+1 + (sign_nu+1)/2;
    return 4*idx[mu][nu]+local;
  }
    
  Matrix<3,3> operator()( int mu , int nu , int sign_mu , int sign_nu ){
    return Matrix<3,3>( m + 9*index(mu,nu,sign_mu,sign_nu) );
  }
  Matrix<3,3> operator()( int mu , int nu , int sign_mu , int sign_nu ) const {
    return Matrix<3,3>( m + 9*index(mu,nu,sign_mu,sign_nu) );
  }

  void Calc( const Gauge& V , const SingleStaple& S ,
	     const Site& x ){
    double _Complex ptr[9];
    Matrix<3,3> tmp(ptr);
    for( int mu=0 ; mu<4 ; mu++ ){
      Site xm = x - mu;
      for( int nu=mu+1 ; nu<4 ; nu++ ){
	tmp = V(x,mu) * Hermite(S(x,mu,nu,1));
	Matrix<3,3>(m+9*index(mu,nu,1,1)) = tmp;
	Matrix<3,3>(m+9*index(nu,mu,1,1)) = Hermite(tmp);
	tmp = V(x,mu) * Hermite(S(x,mu,nu,-1));
	Matrix<3,3>(m+9*index(mu,nu,1,-1)) = tmp;
	Matrix<3,3>(m+9*index(nu,mu,-1,1)) = Hermite(tmp);
	tmp = Hermite(V(xm,mu)) * S(xm,mu,nu,1);
	Matrix<3,3>(m+9*index(mu,nu,-1,1)) = tmp;
	Matrix<3,3>(m+9*index(nu,mu,1,-1)) = Hermite(tmp);
	tmp = Hermite(V(xm,mu)) * S(xm,mu,nu,-1);
	Matrix<3,3>(m+9*index(mu,nu,-1,-1)) = tmp;
	Matrix<3,3>(m+9*index(nu,mu,-1,-1)) = Hermite(tmp);
      }
    }
  }
};

void Calc_LocalPlaquette( const Gauge& V , const SingleStaple& S ){
  double ans = 0.;
#pragma omp parallel reduction( + : ans )
  {
    double _Complex ptr[27];
    Site x(V.lattice());
    ForAllSitesOMP(x){
      LocalPlaquette P( V , S , x );
    }
  }
}

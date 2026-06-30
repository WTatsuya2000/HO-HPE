class LocalElements{

public:
  LocalElements( const Gauge& V , const SingleStaple& S ,
		  const Site& x ){
    Calc_Plaquette( V , S , x );
    Calc_Sigma( V , S , x );
  }
  mutable COMPLEX plaq[9*12*4];
  mutable COMPLEX sigma[9*12*2*4*2*2];

  //Matrix<3,3> P( int mu , int nu , int smu , int snu ){
  //  return Matrix<3,3>( plaq + index(mu,nu,smu,snu) );
  //}
  Matrix<3,3> P( int mu , int nu , int smu , int snu ) const {
    return Matrix<3,3>( plaq + index(mu,nu,smu,snu) );
  }
  Matrix<3,3> Sig( int mu , int nu , int rho , int smu , int snu , int srho) const {
    return Matrix<3,3>( sigma + index(mu,nu,rho,smu,snu,srho) );
  }

  static constexpr int idx2[4][4] = {
    {-1,0,1,2} ,
    {3,-1,4,5} ,
    {6,7,-1,8} ,
    {9,10,11,-1} };
  static constexpr int idx3[4][4][4] = {
    { {-1,-1,-1,-1},
      {-1,-1, 0, 1},
      {-1, 2,-1, 3},
      {-1, 4, 5,-1} },
    { {-1,-1, 6, 7},
      {-1,-1,-1,-1},
      { 8,-1,-1, 9},
      {10,-1,11,-1} },
    { {-1,12,-1,13},
      {14,-1,-1,15},
      {-1,-1,-1,-1},
      {16,17,-1,-1} },
    { {-1,18,19,-1},
      {20,-1,21,-1},
      {22,23,-1,-1},
      {-1,-1,-1,-1} }
  };
  int idx_sign( int s1 , int s2 ) const {
    return 2 * (s1==1) + (s2==1);
  }
  int index( int mu , int nu , int smu , int snu ) const {
    return 9 * ( 4 * idx2[mu][nu] + this->idx_sign(smu,snu) );
  }
  int idx_sign( int s1 , int s2 , int s3 ) const {
    return 4 * (s1==1) + 2 * (s2==1) + (s3==1);
  }
  int index( int mu , int nu , int rho , int smu , int snu , int srho ) const {
    return 9 * ( 8 * idx3[mu][nu][rho] + this->idx_sign(smu,snu,srho) );
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
  void Calc_Plaquette( const Gauge& V , const SingleStaple& S ,
	     const Site& x ){
    double _Complex ptr[9];
    Matrix<3,3> tmp(ptr);
    for( int mu=0 ; mu<4 ; mu++ ){
      Site xm = x - mu;
      for( int nu=mu+1 ; nu<4 ; nu++ ){
	tmp = V(x,mu) * Hermite(S(x,mu,nu,1));
	Matrix<3,3>(plaq+index(mu,nu,1,1)) = tmp;
	Matrix<3,3>(plaq+index(nu,mu,1,1)) = Hermite(tmp);
	tmp = V(x,mu) * Hermite(S(x,mu,nu,-1));
	Matrix<3,3>(plaq+index(mu,nu,1,-1)) = tmp;
	Matrix<3,3>(plaq+index(nu,mu,-1,1)) = Hermite(tmp);
	tmp = Hermite(V(xm,mu)) * S(xm,mu,nu,1);
	Matrix<3,3>(plaq+index(mu,nu,-1,1)) = tmp;
	Matrix<3,3>(plaq+index(nu,mu,1,-1)) = Hermite(tmp);
	tmp = Hermite(V(xm,mu)) * S(xm,mu,nu,-1);
	Matrix<3,3>(plaq+index(mu,nu,-1,-1)) = tmp;
	Matrix<3,3>(plaq+index(nu,mu,-1,-1)) = Hermite(tmp);
      }
    }
  }
  // Sigma: (rho,mu,nu,-rho)
  // (smu,snu)=(+,+)
  //  ^nu
  //  |     |
  //  |     |
  //  | -->-+
  //--+--------> mu
  // (smu,snu)=(-,+)
  //  ^nu
  //  | --<-+
  //  |     |
  //  |     |
  //--+--------> mu
  // (smu,snu)=(-,-)
  //  ^nu
  //  | +-<--
  //  | |
  //  | |
  //--+--------> mu
  // (smu,snu)=(-,+)
  //  ^nu
  //  | |
  //  | |
  //  | +->--
  //--+--------> mu
  void Calc_Sigma( const Gauge& V , const SingleStaple& S ,
	     const Site& x ){
    double _Complex ptr[9];
    Matrix<3,3> tmp(ptr);
    for( int mu=0 ; mu<4 ; mu++ ){
      auto xm = x + mu;
      for( int nu=mu+1 ; nu<4 ; nu++ ){
	auto xn = x + nu;
	for( int rho=0 ; rho<4 ; rho++ ){
	  if( rho==mu or rho==nu ) continue;
	  for( int srho : {-1,1} ){
	    tmp = S(x,mu,rho,srho) * S(xm,nu,rho,srho);
	    Matrix<3,3>(sigma+index(mu,nu,rho, 1 ,1,srho)) = tmp;
	    Matrix<3,3>(sigma+index(nu,mu,rho,-1,-1,srho)) = Hermite(tmp);
	    tmp = S(xm,nu,rho,srho) * Hermite(S(xn,mu,rho,srho));
	    Matrix<3,3>(sigma+index(mu,nu,rho,-1, 1,srho)) = tmp;
	    Matrix<3,3>(sigma+index(nu,mu,rho,-1, 1,srho)) = Hermite(tmp);
	    tmp = Hermite(S(xn,mu,rho,srho)) * Hermite(S(x,nu,rho,srho));
	    Matrix<3,3>(sigma+index(mu,nu,rho,-1,-1,srho)) = tmp;
	    Matrix<3,3>(sigma+index(nu,mu,rho, 1, 1,srho)) = Hermite(tmp);
	    tmp = Hermite(S(x,nu,rho,srho)) * S(x,mu,rho,srho);
	    Matrix<3,3>(sigma+index(mu,nu,rho, 1,-1,srho)) = tmp;
	    Matrix<3,3>(sigma+index(nu,mu,rho, 1,-1,srho)) = Hermite(tmp);
	  }
	}
      }
    }
  }
};

void Calc_LocalElements( const Gauge& V , const SingleStaple& S ){
  double ans=0.;
#pragma omp parallel reduction( + : ans )
  {
    double _Complex ptr[27];
    Site x(V.lattice());
    ForAllSitesOMP(x){
      LocalElements LE( V , S , x );
    }
  }
}

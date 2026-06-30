//const int Dirac = { 0 , -128 , 32 , -32 , 64 , 32 , -64 , -128 , 32 ,32 , 0 , -64 , -32 , -32 , -32 ,-32 , -64 ,-64 , -32 , 16 , -32 , -16 , -16 , -32 , -32 ,-16 };

struct HPE_comb
{
  static constexpr int C2[6][2] = {
    {0,1} , {0,2} , {0,3} , {1,2} , {1,3} , {2,3} 
  };
  static constexpr int P2[12][2] = {
    {0,1} , {1,0} , {0,2} , {2,0} , {0,3} , {3,0} , {1,2} , {2,1} , {1,3} , {3,1} , {2,3} , {3,2} 
  };
  static constexpr int C3[4][3] = { 
    {0,1,2} , {1,2,3} , {2,3,0} , {3,0,1} 
  };
  static constexpr int C1C2[12][3] = { 
    {0,1,2} , {0,2,3} , {0,3,1} , 
    {1,2,3} , {1,3,0} , {1,0,2} , 
    {2,3,0} , {2,0,1} , {2,1,3} , 
    {3,0,1} , {3,1,2} , {3,2,0} 
  };
  static constexpr int C2C1[12][3] = { 
    {0,1,2} , {0,1,3} , {0,2,1} , {0,2,3} , {0,3,1} , {0,3,2} ,
    {1,2,0} , {1,2,3} , {1,3,0} , {1,3,2} , {2,3,0} , {2,3,1}
  };
  static constexpr int P3[24][3] = { 
    {0,1,2} , {0,2,1} , {1,2,0} , {1,0,2} , {2,0,1} , {2,1,0} , 
    {1,2,3} , {1,3,2} , {2,3,1} , {2,1,3} , {3,1,2} , {3,2,1} ,
    {2,3,0} , {2,0,3} , {3,0,2} , {3,2,0} , {0,2,3} , {0,3,2} ,
    {3,0,1} , {3,1,0} , {0,1,3} , {0,3,1} , {1,3,0} , {1,0,3} 
  };
  static constexpr int C22[3][4] = {
    {0,1,2,3} , {0,2,1,3} , {0,3,1,2} 
  };
  static constexpr int C2C1C1[12][4] = {
    {0,1,2,3} , {0,1,3,2} , {0,2,1,3} , {0,2,3,1} , {0,3,1,2} , {0,3,2,1} ,
    {1,2,0,3} , {1,2,3,0} , {1,3,0,2} , {1,3,2,0} , {2,3,0,1} , {2,3,1,0}
  };
  static constexpr int P4[24][4] = {
    {0,1,2,3} , {0,1,3,2} , {0,2,1,3} , {0,2,3,1} , {0,3,1,2} , {0,3,2,1} ,
    {1,0,2,3} , {1,0,3,2} , {1,2,0,3} , {1,2,3,0} , {1,3,0,2} , {1,3,2,0} ,
    {2,0,1,3} , {2,0,3,1} , {2,1,0,3} , {2,1,3,0} , {2,3,0,1} , {2,3,1,0} ,
    {3,0,1,2} , {3,0,2,1} , {3,1,0,2} , {3,1,2,0} , {3,2,0,1} , {3,2,1,0}
  };
  static constexpr int P4R[12][4] = {
    {0,1,2,3} , {1,2,3,0} , {2,3,0,1} , {3,0,1,2} , 
    {0,1,3,2} , {1,3,2,0} , {3,2,0,1} , {2,0,1,3} ,
    {0,2,1,3} , {2,1,3,0} , {1,3,0,2} , {3,0,2,1}
  };
};

double W1_local_elem( const Gauge& V , const SingleStaple& S , const Site& x , const int mu , const int nu , double _Complex *Uptr )
{
  Matrix<3,3> tmp( Uptr );
  Site xm = x+mu;
  Site xn = x+nu;
  tmp = Hermite( S(x,nu,mu,-1) ) * V(x,mu);
  tmp = tmp * S(xm,nu,mu,1);
  tmp = tmp * Hermite( V(xn,mu) );
  return __real__( Trace(tmp) );
}
double W1_local( const Gauge& V , const SingleStaple& S , const Site& x , double _Complex *Uptr )
{
  double ans = 0.;
  for( const auto& m : HPE_comb::P2 )
    ans += W1_local_elem( V , S , x , m[0] , m[1] , Uptr );
  return -128. * ans;
}

double W2_local_elem( const Gauge& V , const SingleStaple& S , const Site& x , const int mu , const int nu , double _Complex *Uptr )
{
  Matrix<3,3> tmp( Uptr );
  tmp = S(x,nu,mu,1) * Hermite( V(x,nu) );
  tmp = tmp * S(x,nu,mu,-1);
  tmp = tmp * Hermite( V(x,nu) );
  return __real__( Trace(tmp) );
}
double W2_local( const Gauge& V , const SingleStaple& S , const Site& x , double _Complex *Uptr )
{
  double ans = 0.;
  for( const auto& m : HPE_comb::P2 )
    ans += W2_local_elem( V , S , x , m[0] , m[1] , Uptr );
  return 32. * ans;
}


double W3_local_elem( const Gauge& V , const SingleStaple& S , const Site& x , const int mu , const int nu , double _Complex *Uptr )
{
  Matrix<3,3> tmp( Uptr );
  Site xn = x-nu;
  tmp = S(x,nu,mu,1) * Hermite( V(x,nu) );
  tmp = tmp * Hermite( S(xn,nu,mu,-1) );
  tmp = tmp * V(xn,nu);
  double ans = __real__( Trace(tmp) );
  tmp = S(x,nu,mu,-1) * Hermite( V(x,nu) );
  tmp = tmp * Hermite( S(xn,nu,mu,1) );
  tmp = tmp * V(xn,nu);
  return ans + __real__( Trace(tmp) );
}
double W3_local( const Gauge& V , const SingleStaple& S , const Site& x , double _Complex *Uptr )
{
  double ans = 0.;
  for( const auto& m : HPE_comb::C2 )
    ans += W3_local_elem( V , S , x , m[0] , m[1] , Uptr );
  return -32. * ans;
}

double W4_local_elem( const Gauge& V , const SingleStaple& S , const Site& x , const int mu , const int nu , double _Complex *Uptr )
{
  Matrix<3,3> tmp( Uptr );
  Site xn = x-nu;
  tmp = S(x,nu,mu,1) * Hermite( V(x,nu) );
  tmp = tmp * Hermite( V(xn,nu) );
  tmp = tmp * S(xn,nu,mu,-1);
  double ans = __real__( Trace(tmp) );
  tmp = S(x,nu,mu,-1) * Hermite( V(x,nu) );
  tmp = tmp * Hermite( V(xn,nu) );
  tmp = tmp * S(xn,nu,mu,1);
  return ans + __real__( Trace(tmp) );
}
double W4_local( const Gauge& V , const SingleStaple& S , const Site& x , double _Complex *Uptr )
{
  double ans = 0.;
  for( const auto& m : HPE_comb::C2 )
    ans += W4_local_elem( V , S , x , m[0] , m[1] , Uptr );
  return 64. * ans;
}

double W5_local_elem( const Gauge& V , const SingleStaple& S , const Site& x , const int mu , const int nu , double _Complex *Uptr )
{
  Matrix<3,3> tmp( Uptr );
  tmp = S(x,nu,mu,1) * Hermite( V(x,nu) );
  tmp = tmp * tmp;
  return __real__( Trace(tmp) );
}
double W5_local( const Gauge& V , const SingleStaple& S , const Site& x , double _Complex *Uptr )
{
  double ans = 0.;
  for( const auto& m : HPE_comb::C2 )
    ans += W5_local_elem( V , S , x , m[0] , m[1] , Uptr );
  return 32. * ans / 2.;
}

double W6_local_elem( const Gauge& V , const SingleStaple& S , const Site& x , const int mu , const int nu , double _Complex *Uptr )
{
  Matrix<3,3> tmp( Uptr );
  Site xm = x-mu;
  double ans = 0.;
  for( int s : {-1,1} ){
    tmp = V(x,mu) * Hermite( S(x,mu,nu,-s) );
    tmp = tmp * Hermite( S(xm,mu,nu,-s) );
    tmp = tmp * S(xm,mu,nu,s);
    ans += __real__( Trace(tmp) );
    tmp = S(x,mu,nu,s) * Hermite( S(x,mu,nu,-s) );
    tmp = tmp * Hermite( S(xm,mu,nu,-s) );
    tmp = tmp * V(xm,mu);
    ans += __real__( Trace(tmp) );
  }
  return ans;
}
double W6_local( const Gauge& V , const SingleStaple& S , const Site& x , double _Complex *Uptr )
{
  double ans = 0.;
  for( const auto& m : HPE_comb::C2 )
    ans += W6_local_elem( V , S , x , m[0] , m[1] , Uptr );
  return -64. * ans;
}

double W7_local_elem( const Gauge& V , const SingleStaple& S , const Site& x , const int mu , const int nu , double _Complex *Uptr )
{
  Matrix<3,3> tmp( Uptr );
  Site xm = x-mu;
  tmp = S(x,mu,nu,-1) * Hermite( S(x,mu,nu,1) );
  tmp = tmp * Hermite( S(xm,mu,nu,1) );
  tmp = tmp * S(xm,mu,nu,-1);
  return __real__( Trace(tmp) );
}
double W7_local( const Gauge& V , const SingleStaple& S , const Site& x , double _Complex *Uptr )
{
  double ans = 0.;
  for( const auto& m : HPE_comb::C2 )
    ans += W7_local_elem( V , S , x , m[0] , m[1] , Uptr );
  return -128. * ans;
}

double W8_local_elem( const Gauge& V , const SingleStaple& S , const Site& x , const int mu , const int nu , const int rho , double _Complex *Uptr )
{
  Matrix<3,3> tmp( Uptr );
  double ans = 0.;
  for( int sn : {-1,1} ){
    for( int sr : {-1,1} ){
      tmp = S(x,mu,nu,sn) * Hermite( V(x,mu) );
      tmp = tmp * S(x,mu,rho,sr);
      tmp = tmp * Hermite( V(x,mu) );
      ans += __real__( Trace(tmp) );
    }
  }
  return ans;
}
double W8_local( const Gauge& V , const SingleStaple& S , const Site& x , double _Complex *Uptr )
{
  double ans = 0.;
  for( const auto& m : HPE_comb::C1C2 )
    ans += W8_local_elem( V , S , x , m[0] , m[1] , m[2] , Uptr );
  return 32. * ans;
}

double W9_local_elem( const Gauge& V , const SingleStaple& S , const Site& x , const int mu , const int nu , const int rho , double _Complex *Uptr )
{
  Matrix<3,3> tmp( Uptr );
  Site xm = x - mu;
  double ans = 0.;
  for( int sn : {-1,1} ){
    for( int sr : {-1,1} ){
      tmp = V(x,mu) * Hermite( S(x,mu,nu,sn) );
      tmp = tmp * Hermite( S(xm,mu,rho,sr) );
      tmp = tmp * V(xm,mu);
      ans += __real__( Trace(tmp) );
    }
  }
  return ans;
}
double W9_local( const Gauge& V , const SingleStaple& S , const Site& x , double _Complex *Uptr )
{
  double ans = 0.;
  for( const auto& m : HPE_comb::P3 )
    ans += W9_local_elem( V , S , x , m[0] , m[1] , m[2] , Uptr );
  return 32. * ans;
}

double W10_local( const Gauge& V , const SingleStaple& S , const Site& x , double _Complex *Uptr )
{
  return 0.;
}

double W11_local_elem( const Gauge& V , const SingleStaple& S , const Site& x , const int mu , const int nu , const int rho , double _Complex *Uptr )
{
  Matrix<3,3> tmp( Uptr );
  Site xm = x + mu;
  Site xn = x + nu;
  double ans = 0.;
  for ( int s : {-1,1} ){
    tmp = V(x,mu) * S(xm,nu,mu,1);
    tmp = tmp * Hermite( V(xn,mu) );
    tmp = tmp * Hermite( S(x,nu,rho,s) );
    ans += __real__( Trace(tmp) );
    tmp = V(x,mu) * S(xm,nu,rho,s);
    tmp = tmp * Hermite( V(xn,mu) );
    tmp = tmp * Hermite( S(x,nu,mu,-1) );
    ans += __real__( Trace(tmp) );
  }
  return ans;
}
double W11_local( const Gauge& V , const SingleStaple& S , const Site& x , double _Complex *Uptr )
{
  double ans = 0.;
  for( const auto& m : HPE_comb::P3 )
    ans += W11_local_elem( V , S , x , m[0] , m[1] , m[2] , Uptr );
  return -64. * ans;
}

double W12_local_elem( const Gauge& V , const SingleStaple& S , const Site& x , const int mu , const int nu , const int rho , double _Complex *Uptr )
{
  Matrix<3,3> tmp( Uptr );
  Site xm = x + mu;
  Site xn = x + nu;
  double ans = 0.;
  for( int s : {-1,1} ){
    tmp = V(x,mu) * S(xm,nu,rho,s);
    tmp = tmp * Hermite( V(xn,mu) );
    tmp = tmp * Hermite( S(x,nu,rho,s) );
    ans += __real__( Trace(tmp) );
  }
  return ans;
}
double W12_local( const Gauge& V , const SingleStaple& S , const Site& x , double _Complex *Uptr )
{
  double ans = 0.;
  for( const auto& m : HPE_comb::C2C1 )
    ans += W12_local_elem( V , S , x , m[0] , m[1] , m[2] , Uptr );
  return -32. * ans;
}

double W13_local_elem( const Gauge& V , const SingleStaple& S , const Site& x , const int mu , const int nu , const int rho , double _Complex *Uptr )
{
  Matrix<3,3> tmp( Uptr );
  Site xm = x + mu;
  Site xn = x + nu;
  double ans = 0.;
  for( int s : {-1,1} ){
    tmp = V(x,mu) * S(xm,nu,rho,s);
    tmp = tmp * Hermite( V(xn,mu) );
    tmp = tmp * Hermite( S(x,nu,rho,-s) );
    ans += __real__( Trace(tmp) );
  }
  return ans;
}
double W13_local( const Gauge& V , const SingleStaple& S , const Site& x , double _Complex *Uptr )
{
  double ans = 0.;
  for( const auto& m : HPE_comb::P3 )
    ans += W13_local_elem( V , S , x , m[0] , m[1] , m[2] , Uptr );
  return -32. * ans;
}

double W14_local_elem( const Gauge& V , const SingleStaple& S , const Site& x , const int mu , const int nu , const int rho , double _Complex *Uptr )
{
  Matrix<3,3> tmp( Uptr );
  Site xm = x - mu;
  Site xn = x - nu;
  double ans = 0.;
  for( int sr : {-1,1} ){
    for( int sn : {-1,1} ){
      tmp = V(x,mu) * Hermite( S(x,mu,nu,sn) );
      tmp = tmp * Hermite( S(xm,mu,nu,sn) );
      tmp = tmp * S(xm,mu,rho,sr);
      ans += __real__( Trace(tmp) );
      tmp = S(x,mu,rho,sr) * Hermite( S(x,mu,nu,sn) );
      tmp = tmp * Hermite( S(xm,mu,nu,sn) );
      tmp = tmp * V(xm,mu);
      ans += __real__( Trace(tmp) );
    }
  }
  return ans;
}
double W14_local( const Gauge& V , const SingleStaple& S , const Site& x , double _Complex *Uptr )
{
  double ans = 0.;
  for( const auto& m : HPE_comb::P3 )
    ans += W14_local_elem( V , S , x , m[0] , m[1] , m[2] , Uptr );
  return -32. * ans;
}

double W15_local_elem( const Gauge& V , const SingleStaple& S , const Site& x , const int mu , const int nu , const int rho , double _Complex *Uptr )
{
  Matrix<3,3> tmp( Uptr );
  Site xm = x + mu;
  Site xn = x + nu;
  double ans = 0.;
  for( int sr : {-1,1} ){
    tmp = Hermite( V(x,nu) ) * V(x,mu);
    tmp = tmp * S(xm,nu,rho,sr);
    tmp = tmp * Hermite( S(xn,mu,rho,-sr) );
    ans += __real__( Trace(tmp) );
    tmp = V(x,mu) * V(xm,nu);
    tmp = tmp * Hermite( S(xn,mu,rho,sr) );
    tmp = tmp * Hermite( S(x,nu,rho,-sr) );
    ans += __real__( Trace(tmp) );
    tmp = V(xm,nu) * Hermite( V(xn,mu) );
    tmp = tmp * Hermite( S(x,nu,rho,sr));
    tmp = tmp * S(x,mu,rho,-sr);
    ans += __real__( Trace(tmp) );
    tmp = Hermite( V(xn,mu) ) * Hermite( V(x,nu));
    tmp = tmp * S(x,mu,rho,sr);
    tmp = tmp * S(xm,nu,rho,-sr);
    ans += __real__( Trace(tmp) );
  }
  return ans;
}
double W15_local( const Gauge& V , const SingleStaple& S , const Site& x , double _Complex *Uptr )
{
  double ans = 0.;
  for( const auto& m : HPE_comb::C2C1 )
    ans += W15_local_elem( V , S , x , m[0] , m[1] , m[2] , Uptr );
  return -32. * ans;
}

double W16_local_elem( const Gauge& V , const SingleStaple& S , const Site& x , const int mu , const int nu , const int rho , double _Complex *Uptr )
{
  Matrix<3,3> tmp( Uptr );
  Site xm = x - mu;
  Site xn = x - nu;
  double ans = 0.;
  for( int sn : {-1,1} ){
    for( int sr : {-1,1} ){
      tmp = S(x,mu,nu,sn) * Hermite( S(x,mu,rho,sr) );
      tmp = tmp * Hermite( S(xm,mu,rho,sr) );
      tmp = tmp * S(xm,mu,nu,sn);
      ans += __real__( Trace(tmp) );
    }
  }
  return ans;
}
double W16_local( const Gauge& V , const SingleStaple& S , const Site& x , double _Complex *Uptr )
{
  double ans = 0.;
  for( const auto& m : HPE_comb::C1C2 )
    ans += W16_local_elem( V , S , x , m[0] , m[1] , m[2] , Uptr );
  return -64. * ans;
}

double W17_local_elem( const Gauge& V , const SingleStaple& S , const Site& x , const int mu , const int nu , const int rho , double _Complex *Uptr )
{
  Matrix<3,3> tmp( Uptr );
  double ans = 0.;
  Site xn = x + nu;
  Site xr = x + rho;
  for( int s : {-1,1} ){
    tmp = S(x,nu,mu,s) * S(xn,rho,mu,s);
    tmp = tmp * Hermite( S(xr,nu,mu,-s) );
    tmp = tmp * Hermite( S(x,rho,mu,-s) );
    ans += __real__( Trace(tmp) );
    tmp = S(xn,rho,mu,s) * Hermite( S(xr,nu,mu,s) );
    tmp = tmp * Hermite( S(x,rho,mu,-s) );
    tmp = tmp * S(x,nu,mu,-s);
    ans += __real__( Trace(tmp) );
  }
  return ans;
}
double W17_local( const Gauge& V , const SingleStaple& S , const Site& x , double _Complex *Uptr )
{
  double ans = 0.;
  for( const auto& m : HPE_comb::C1C2 )
    ans += W17_local_elem( V , S , x , m[0] , m[1] , m[2] , Uptr );
  return -64. * ans;
}

double W18_local_elem( const Gauge& V , const SingleStaple& S , const Site& x , const int mu , const int nu , const int rho , double _Complex *Uptr )
{
  Matrix<3,3> tmp( Uptr );
  double ans = 0.;
  Site xn = x + nu;
  Site xr = x + rho;
  for( int s : {-1,1} ){
    tmp = S(x,nu,mu,s) * V(xn,rho);
    tmp = tmp * Hermite( S(xr,nu,mu,-s) );
    tmp = tmp * Hermite( S(x,rho,mu,-s) );
    ans += __real__( Trace(tmp) );
    tmp = S(xn,rho,mu,s) * Hermite( V(xr,nu) );
    tmp = tmp * Hermite( S(x,rho,mu,-s) );
    tmp = tmp * S(x,nu,mu,-s);
    ans += __real__( Trace(tmp) );
    tmp = Hermite( S(xr,nu,mu,s) ) * Hermite( V(x,rho) );
    tmp = tmp * S(x,nu,mu,-s);
    tmp = tmp * S(xn,rho,mu,-s);
    ans += __real__( Trace(tmp) );
    tmp = Hermite( S(x,rho,mu,s) ) * V(x,nu);
    tmp = tmp * S(xn,rho,mu,-s);
    tmp = tmp * Hermite( S(xr,nu,mu,-s) );
    ans += __real__( Trace(tmp) );
  }
  return ans;
}
double W18_local( const Gauge& V , const SingleStaple& S , const Site& x , double _Complex *Uptr )
{
  double ans = 0.;
  for( const auto& m : HPE_comb::P3 )
    ans += W18_local_elem( V , S , x , m[0] , m[1] , m[2] , Uptr );
  return -32. * ans;
}

bool apneq(double x, double y, double eps = 1e-12) {
    return std::abs(x - y) > eps;
}
double W19_local_elem( const Gauge& V , const SingleStaple& S , const Site& x , const int mu , const int nu , const int rho , const int sig , double _Complex *Uptr )
{
  Matrix<3,3> tmp1( Uptr ) , tmp2( Uptr+9 ), tmp3( Uptr+18 );
  Site xm = x - mu;
  Site xr = x - rho;
  double ans = 0.;
  tmp1 = S(x,mu,nu,1) + S(x,mu,nu,-1);
  tmp2 = tmp1 * Hermite( V(x,mu) );
  tmp1 = Hermite( S(xm,mu,nu,1) ) + Hermite( S(xm,mu,nu,-1) );
  tmp2 += tmp1 * V(xm,mu);
  tmp1 = S(x,rho,sig,1) + S(x,rho,sig,-1);
  tmp3 = tmp1 * Hermite( V(x,rho) );
  tmp1 = Hermite( S(xr,rho,sig,1) ) + Hermite( S(xr,rho,sig,-1) );
  tmp3 += tmp1 * V(xr,rho);
  tmp3 = tmp3 + Hermite( tmp3 );
  tmp1 = tmp2 * tmp3;
  ans += __real__( Trace(tmp1) );
  return ans;
}
double W19_local( const Gauge& V , const SingleStaple& S , const Site& x , double _Complex *Uptr )
{
  double ans = 0.;
  for( const auto& m : HPE_comb::C22 )
    ans += W19_local_elem( V , S , x , m[0] , m[1] , m[2] , m[3] , Uptr );
  return 16. * ans;
}

double W20_local_elem( const Gauge& V , const SingleStaple& S , const Site& x , const int mu , const int nu , const int rho , const int sig , double _Complex *Uptr )
{
  Matrix<3,3> tmp( Uptr );
  Site xm = x + mu;
  Site xn = x + nu;
  double ans = 0.;
  for( int sr : {-1,1} ){
    for( int ss : {-1,1} ){
      tmp = V(x,mu) * S(xm,nu,rho,sr);
      tmp = tmp * Hermite( V(xn,mu) );
      tmp = tmp * Hermite( S(x,nu,sig,ss));
      ans += __real__( Trace(tmp) );
    }
  }
  return ans;
}
double W20_local( const Gauge& V , const SingleStaple& S , const Site& x , double _Complex *Uptr )
{
  double ans = 0.;
  for( const auto& m : HPE_comb::P4 )
    ans += W20_local_elem( V , S , x , m[0] , m[1] , m[2] , m[3] , Uptr );
  return -32. * ans;
}

double W21_local_elem( const Gauge& V , const SingleStaple& S , const Site& x , const int mu , const int nu , const int rho , const int sig , double _Complex *Uptr )
{
  Matrix<3,3> tmp( Uptr );
  Site xm = x + mu;
  Site xn = x + nu;
  double ans = 0.;
  for( int sr : {-1,1} ){
    for( int ss : {-1,1} ){
      tmp = V(x,mu) * V(xm,nu);
      tmp = tmp * Hermite( S(xn,mu,rho,sr) );
      tmp = tmp * Hermite( S(x,nu,sig,ss));
      ans += __real__( Trace(tmp) );
      tmp = V(xm,nu) * Hermite( V(xn,mu) );
      tmp = tmp * Hermite( S(x,nu,rho,sr));
      tmp = tmp * S(x,mu,sig,ss);
      ans += __real__( Trace(tmp) );
      tmp = Hermite( V(xn,mu) ) * Hermite( V(x,nu));
      tmp = tmp * S(x,mu,rho,sr);
      tmp = tmp * S(xm,nu,sig,ss);
      ans += __real__( Trace(tmp) );
      tmp = Hermite( V(x,nu) ) * V(x,mu);
      tmp = tmp * S(xm,nu,rho,sr);
      tmp = tmp * Hermite( S(xn,mu,sig,ss) );
      ans += __real__( Trace(tmp) );
    }
  }
  return ans;
}
double W21_local( const Gauge& V , const SingleStaple& S , const Site& x , double _Complex *Uptr )
{
  double ans = 0.;
  for( const auto& m : HPE_comb::C2C1C1 )
    ans += W21_local_elem( V , S , x , m[0] , m[1] , m[2] , m[3] , Uptr );
  return -16. * ans;
}

double W22_local_elem( const Gauge& V , const SingleStaple& S , const Site& x , const int mu , const int nu , const int rho , const int sig , double _Complex *Uptr )
{
  Matrix<3,3> tmp( Uptr );
  double ans = 0.;
  Site xm = x - mu;
  Site xn = x - nu;
  for( int sr : {-1,1} ){
    for( int ss : {-1,1} ){
      tmp = Hermite( S(x,nu,rho,sr) ) * S(x,mu,rho,sr);
      tmp = tmp * Hermite( S(x,mu,sig,ss) );
      tmp = tmp * S(x,nu,sig,ss);
      ans += __real__( Trace(tmp) );
      tmp = S(xm,mu,rho,sr) * S(x,nu,rho,sr);
      tmp = tmp * Hermite( S(x,nu,sig,ss) );
      tmp = tmp * Hermite( S(xm,mu,sig,ss) );
      ans += __real__( Trace(tmp) );
      tmp = S(xn,nu,rho,sr) * Hermite( S(xm,mu,rho,sr) );
      tmp = tmp * S(xm,mu,sig,ss);
      tmp = tmp * Hermite( S(xn,nu,sig,ss) );
      ans += __real__( Trace(tmp) );
      tmp = Hermite( S(x,mu,rho,sr) ) * Hermite( S(xn,nu,rho,sr) );
      tmp = tmp * S(xn,nu,sig,ss);
      tmp = tmp * S(x,mu,sig,ss);
      ans += __real__( Trace(tmp) );
    }
  }
  return ans;
}
double W22_local( const Gauge& V , const SingleStaple& S , const Site& x , double _Complex *Uptr )
{
  double ans = 0.;
  for( const auto& m : HPE_comb::C22 )
    ans += W22_local_elem( V , S , x , m[0] , m[1] , m[2] , m[3] , Uptr );
 return -32. * ans;
}


double W23_local_elem( const Gauge& V , const SingleStaple& S , const Site& x , const int mu , const int nu , const int rho , const int sig , double _Complex *Uptr )
{
  Matrix<3,3> tmp( Uptr );
  double ans = 0.;
  Site xm = x + mu;
  Site xn = x + nu;
  for( int sr : {-1,1} ){
    for( int ss : {-1,1} ){
      tmp = S(x,mu,rho,sr) * S(xm,nu,rho,sr);
      tmp = tmp * Hermite( S(xn,mu,sig,ss) );
      tmp = tmp * Hermite( S(x,nu,sig,ss) );
      ans += __real__( Trace(tmp) );
      tmp = S(x,mu,sig,ss) * S(xm,nu,rho,sr);
      tmp = tmp * Hermite( S(xn,mu,rho,sr) );
      tmp = tmp * Hermite( S(x,nu,sig,ss) );
      ans += __real__( Trace(tmp) );
    }
  }
  return ans;
}
double W23_local( const Gauge& V , const SingleStaple& S , const Site& x , double _Complex *Uptr )
{
  double ans = 0.;
  for( const auto& m : HPE_comb::C2C1C1 )
    ans += W23_local_elem( V , S , x , m[0] , m[1] , m[2] , m[3] , Uptr );
  return -32. * ans;
}

double W24_local_elem( const Gauge& V , const SingleStaple& S , const Site& x , const int mu , const int nu , const int rho , const int sig , double _Complex *Uptr )
{
  Matrix<3,3> tmp( Uptr );
  double ans = 0.;
  Site xm = x + mu;
  Site xn = x + nu;
  for( int sr : {-1,1} ){
    for( int ss : {-1,1} ){
      tmp = S(x,mu,sig,ss) * V(xm,nu);
      tmp = tmp * Hermite( S(xn,mu,rho,sr) );
      tmp = tmp * Hermite( S(x,nu,rho,sr) );
      ans += __real__( Trace(tmp) );
      tmp = S(x,mu,rho,sr) * S(xm,nu,sig,ss);
      tmp = tmp * Hermite( V(xn,mu) );
      tmp = tmp * Hermite( S(x,nu,rho,sr) );
      ans += __real__( Trace(tmp) );
      tmp = S(x,mu,rho,sr) * S(xm,nu,rho,sr);
      tmp = tmp * Hermite( S(xn,mu,sig,ss) );
      tmp = tmp * Hermite( V(x,nu) );
      ans += __real__( Trace(tmp) );
      tmp = V(x,mu) * S(xm,nu,rho,sr);
      tmp = tmp * Hermite( S(xn,mu,rho,sr) );
      tmp = tmp * Hermite( S(x,nu,sig,ss) );
      ans += __real__( Trace(tmp) );
    }
  }
  return ans/2.;//double counting
}
double W24_local( const Gauge& V , const SingleStaple& S , const Site& x , double _Complex *Uptr )
{
  double ans = 0.;
  //for( const auto& m : HPE_comb::C2C1C1 )
  for( const auto& m : HPE_comb::P4 )
    ans += W24_local_elem( V , S , x , m[0] , m[1] , m[2] , m[3] , Uptr );
  return -16. * ans;
}

double W25_local_elem( const Gauge& V , const SingleStaple& S , const Site& x , const int mu , const int nu , const int rho , const int sig , double _Complex *Uptr )
{
  Matrix<3,3> tmp( Uptr );
  double ans = 0.;
  Site x1 = x + mu;
  Site x2 = x1 + nu;
  Site x3 = x2 + rho;
  Site x7 = x + sig;
  Site x6 = x7 + rho;
  Site x5 = x6 + nu;
  tmp = V(x,mu) * V(x1,nu);
  tmp = tmp * V(x2,rho);
  tmp = tmp * V(x3,sig);
  tmp = tmp * Hermite( V(x5,mu) );
  tmp = tmp * Hermite( V(x6,nu) );
  tmp = tmp * Hermite( V(x7,rho) );
  tmp = tmp * Hermite( V(x,sig) );
  ans += __real__( Trace(tmp) );
  x2 = x1 - nu;
  x3 = x2 + rho;
  x5 = x6 - nu;
  tmp = V(x,mu) * Hermite( V(x2,nu) );
  tmp = tmp * V(x2,rho);
  tmp = tmp * V(x3,sig);
  tmp = tmp * Hermite( V(x5,mu) );
  tmp = tmp * V(x5,nu);
  tmp = tmp * Hermite( V(x7,rho) );
  tmp = tmp * Hermite( V(x,sig) );
  ans += __real__( Trace(tmp) );
  return ans;
}
double W25_local( const Gauge& V , const SingleStaple& S , const Site& x , double _Complex *Uptr )
{
  double ans = 0.;
  for( const auto& m : HPE_comb::P4R )
    ans += W25_local_elem( V , S , x , m[0] , m[1] , m[2] , m[3] , Uptr );
  //for( const auto& m : HPE_comb::P4 )
  //  ans += W25_local_elem( V , S , x , m[0] , m[1] , m[2] , m[3] , Uptr )/2.;
  return -16. * ans;
}

double W_local( const Gauge& V , const SingleStaple& S , const Site& x , double _Complex *Uptr , int N ){
  using Wfunc = double(*)( const Gauge& , const SingleStaple& , const Site& , double _Complex* );
  Wfunc Wlist[] = { 
    W1_local , W2_local , W3_local , W4_local , W5_local ,
    W6_local , W7_local , W8_local , W9_local , W10_local ,
    W11_local , W12_local , W13_local , W14_local , W15_local ,
    W16_local , W17_local , W18_local , W19_local , W20_local ,
    W21_local , W22_local , W23_local , W24_local , W25_local 
  };
  return Wlist[N]( V , S , x , Uptr );
}

double NNLO_W( const Gauge& V , const SingleStaple& S , int N ){
  double ans = 0.;
#pragma omp parallel reduction( + : ans )
  {
    double _Complex ptr[27];
    Site x(V.lattice());
    ForAllSitesOMP(x){
      ans += W_local( V , S , x , ptr , N );
    }
  }
  mpi.add(ans);
  return -2. * ans / V.lattice().size();
}

double NNLO_Wall( const Gauge& V , const SingleStaple& S ){
  double ans = 0.;
#pragma omp parallel reduction( + : ans )
  {
    double _Complex ptr[27];
    Site x(V.lattice());
    ForAllSitesOMP(x){
      for( int i=0 ; i<25 ; i++ )
        ans += W_local( V , S , x , ptr , i );
    }
  }
  mpi.add(ans);
  return -2. * ans / V.lattice().size();
}

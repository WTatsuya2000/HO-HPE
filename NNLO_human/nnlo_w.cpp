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
  static constexpr int C2P2[12][4] = {
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
  static constexpr int remnant[4][4][2] = {
    // mu, nu -> (rho, sig)
    {{-1,-1},{2,3},{1,3},{1,2}},
    {{2,3},{-1,-1},{0,3},{0,2}},
    {{1,3},{0,3},{-1,-1},{0,1}},
    {{1,2},{0,2},{0,1},{-1,-1}}
  };
};
template< typename m1 , typename m2 >
double RTP( const m1& A , const m2& B )
{
  return __real__( TraceProd( A , B ) );
}

double W1_11_12_13_20_local( const Gauge& V , const SingleStaple& S , const LocalElements& LE , const Site& x , double _Complex *Uptr )
{
  Matrix<3,3> tmp1( Uptr ) , tmp2( Uptr+9 ) , tmp3(Uptr+18);
  double ans = 0.;
  for( const auto& m : HPE_comb::P2 ){
    const int mu=m[0] , nu=m[1];
    const int rho = HPE_comb::remnant[mu][nu][0] , sig = HPE_comb::remnant[mu][nu][1];
    Site xm = x+mu;
    Site xn = x+nu;

    tmp1  = S(x,nu,rho,0) + S(x,nu,sig,0);
    tmp1 += 2.* S(x,nu,mu,-1);
    tmp2  = S(xm,nu,rho,1) + S(xm,nu,sig,1);
    tmp2 += 2.* S(xm,nu,mu,1);
    tmp2  = V(x,mu) * tmp2;
    tmp2  = tmp2 * Hermite(V(xn,mu));
    ans += -32.* RTP( tmp1 , Hermite(tmp2) );

    tmp1  = S(x,nu,rho,1) + S(x,nu,sig,1);
    tmp1 += 2.* S(x,nu,mu,-1);
    tmp3 = tmp1 + S(x,nu,rho,-1);
    tmp2 = V(x,mu) * S(xm,nu,sig,-1);
    tmp2 = tmp2 * Hermite(V(xn,mu));
    ans += -32.* RTP( tmp3 , Hermite(tmp2) );
    tmp3 = tmp1 + S(x,nu,sig,-1);
    tmp2 = V(x,mu) * S(xm,nu,rho,-1);
    tmp2 = tmp2 * Hermite(V(xn,mu));
    ans += -32.* RTP( tmp3 , Hermite(tmp2) );
  }
  return ans;
}

double W2_5_local( const Gauge& V , const SingleStaple& S , const LocalElements& LE , const Site& x , double _Complex *Uptr )
{
  Matrix<3,3> tmp( Uptr );
  double ans = 0.;
  for( const auto& m : HPE_comb::C2 ){
    const int mu=m[0] , nu=m[1];
    tmp = LE.P(mu,nu,-1,1) + LE.P(mu,nu,1,-1);
    tmp += .5 * LE.P(mu,nu,1,1);
    ans += RTP( LE.P(mu,nu,1,1) , tmp );
  }
  return 32. * ans;
}

double W3_local( const Gauge& V , const SingleStaple& S , const LocalElements& LE , const Site& x , double _Complex *Uptr )
{
  double ans = 0.;
  for( const auto& m : HPE_comb::C2 ){
    const int mu=m[0] , nu=m[1];
    Site xn = x-nu;
    ans += RTP( LE.P(mu,nu,1,1) , LE.P(mu,nu,-1,-1) );
    ans += RTP( LE.P(mu,nu,-1,1) , LE.P(mu,nu,1,-1) );
  }
  return -32. * ans;
}

double W4_local( const Gauge& V , const SingleStaple& S , const LocalElements& LE , const Site& x , double _Complex *Uptr )
{
  double ans = 0.;
  for( const auto& m : HPE_comb::C2 ){
    const int mu=m[0] , nu=m[1];
    Site xn = x-nu;
    ans += RTP( LE.P(mu,nu,1,1) , LE.P(nu,mu,-1,-1) );
    ans += RTP( LE.P(mu,nu,-1,1) , LE.P(nu,mu,-1,1) );
  }
  return 64. * ans;
}

double W6_7_local( const Gauge& V , const SingleStaple& S , const LocalElements& LE , const Site& x , double _Complex *Uptr )
{
  Matrix<3,3> tmp1( Uptr ) , tmp2( Uptr+9 ), tmp3( Uptr+18 );
  double ans = 0.;
  for( const auto& m : HPE_comb::C2 ){
    const int mu=m[0] , nu=m[1];
    tmp1 = LE.P(mu,nu,1,1) * LE.P(nu,mu,1,-1);
    tmp2 = LE.P(mu,nu,-1,-1) * LE.P(nu,mu,-1,1);
    tmp3 = LE.P(mu,nu,1,1) + LE.P(nu,mu,1,-1);
    tmp3 += tmp1;
    ans += RTP( tmp3 , tmp2 );
    tmp3 = LE.P(mu,nu,-1,-1) + LE.P(nu,mu,-1,1);
    tmp3 += tmp2;
    ans += RTP( tmp1 , tmp3 );
  }
  return -64. * ans;
}

double W8_local( const Gauge& V , const SingleStaple& S , const LocalElements& LE , const Site& x , double _Complex *Uptr )
{
  double ans = 0.;
  for( const auto& m : HPE_comb::C2C1 ){
    const int mu=m[0] , nu=m[1] , rho=m[2];
    ans += RTP( LE.P(rho,mu,1,0) , LE.P(rho,nu,1,0) );
  }
  return 32. * ans;
}

double W9_local( const Gauge& V , const SingleStaple& S , const LocalElements& LE , const Site& x , double _Complex *Uptr )
{
  double ans = 0.;
  for( const auto& m : HPE_comb::P3 ){
    const int mu=m[0] , nu=m[1] , rho=m[2];
    ans += RTP( LE.P(mu,nu,1,0) , LE.P(rho,mu,0,-1) );
  }
  return 32. * ans;
}

double W14_local( const Gauge& V , const SingleStaple& S , const LocalElements& LE , const Site& x , double _Complex *Uptr )
{
  Matrix<3,3> tmp( Uptr ) , tmp1( Uptr+9 );
  double ans = 0.;
  for( const auto& m : HPE_comb::P3 ){
    const int mu=m[0] , nu=m[1] , rho=m[2];
    Site xm = x - mu;
    tmp  = LE.P(mu,nu,1,1) * LE.P(nu,mu,1,-1);
    tmp += LE.P(mu,nu,1,-1) * LE.P(nu,mu,-1,-1);
    tmp1 = LE.P(rho,mu,0,1) + LE.P(mu,rho,-1,0);
    ans += RTP( tmp , tmp1 );
  }
  return -32. * ans;
}

double W15_21_local( const Gauge& V , const SingleStaple& S , const LocalElements& LE , const Site& x , double _Complex *Uptr )
{
  Matrix<3,3> tmp1( Uptr ) , tmp2( Uptr+9 );
  double ans = 0.;
  for( const auto& m : HPE_comb::C2){
    const int mu=m[0] , nu=m[1];
    const int
      rho = HPE_comb::remnant[mu][nu][0] ,
      sig = HPE_comb::remnant[mu][nu][1];
    for( int sn : {-1,1} ){
      for( int s : {-1,1} ){
        tmp1  = LE.P(sig,mu,0,1);
        tmp1 += 2.* LE.P(rho,mu,-s,1);
        tmp2  = tmp1 * LE.P(mu,nu,1,sn);
        tmp1  = LE.P(sig,mu,0,-1);
        tmp1 += 2.* LE.P(rho,mu,-s,-1);
        tmp2 += tmp1 * LE.P(mu,nu,-1,sn);
        ans += RTP( tmp2 , LE.P(nu,rho,sn,s) );
    
        tmp1  = LE.P(rho,mu,0,1);
        tmp1 += 2.* LE.P(sig,mu,-s,1);
        tmp2  = tmp1 * LE.P(mu,nu,1,sn);
        tmp1  = LE.P(rho,mu,0,-1);
        tmp1 += 2.* LE.P(sig,mu,-s,-1);
        tmp2 += tmp1 * LE.P(mu,nu,-1,sn);
        ans += RTP( tmp2 , LE.P(nu,sig,sn,s) );
      }
    }
  }
  return -16. * ans;
}

double W16_17_local( const Gauge& V , const SingleStaple& S , const LocalElements& LE , const Site& x , double _Complex *Uptr )
{
  Matrix<3,3> tmp( Uptr );
  double ans = 0.;
  for( const auto& m : HPE_comb::C2P2 ){
    const int mu=m[0] , nu=m[1] , rho=m[2];
    tmp = LE.Sig(nu,mu,rho,-1,-1,-1) + LE.Sig(mu,nu,rho,-1,-1,-1);
    ans += RTP( LE.Sig(mu,nu,rho, 1, 1,1) , tmp );
    tmp = LE.Sig(nu,mu,rho,-1, 1,-1) + LE.Sig(mu,nu,rho, 1,-1,-1);
    ans += RTP( LE.Sig(mu,nu,rho,-1, 1,1) , tmp );
    tmp = LE.Sig(nu,mu,rho, 1, 1,-1) + LE.Sig(mu,nu,rho, 1, 1,-1);    
    ans += RTP( LE.Sig(mu,nu,rho,-1,-1,1) , tmp );
    tmp = LE.Sig(nu,mu,rho, 1,-1,-1) + LE.Sig(mu,nu,rho,-1, 1,-1);
    ans += RTP( LE.Sig(mu,nu,rho, 1,-1,1) , tmp );
  }
  return -64. * ans;
}

double W18_local( const Gauge& V , const SingleStaple& S , const LocalElements& LE , const Site& x , double _Complex *Uptr )
{
  Matrix<3,3> tmp( Uptr ) , tmp1( Uptr+9 );
  double ans = 0.;
  for( const auto& m : HPE_comb::C2){
    const int mu=m[0] , nu=m[1];
    const int
      rho = HPE_comb::remnant[mu][nu][0] ,
      sig = HPE_comb::remnant[mu][nu][1];
    Site xm = x + mu;
    Site xn = x + nu;
    tmp  = S(x,mu,nu,-1) * V(xm,nu);
    tmp += V(x,mu) * S(xm,nu,mu,1);
    tmp1 = LE.Sig(mu,nu,rho,-1,-1,0) + LE.Sig(mu,nu,sig,-1,-1,0);
    ans += RTP( tmp , tmp1 );

    tmp  = S(xm,nu,mu,1) * Hermite(V(xn,mu));
    tmp += V(xm,nu) * Hermite(S(xn,mu,nu,1));
    tmp1 = LE.Sig(mu,nu,rho,1,-1,0) + LE.Sig(mu,nu,sig,1,-1,0);
    ans += RTP( tmp , tmp1 );

    tmp  = Hermite(S(xn,mu,nu,1)) * Hermite(V(x,nu));
    tmp += Hermite(V(xn,mu)) * Hermite(S(x,nu,mu,-1));
    tmp1 = LE.Sig(mu,nu,rho,1,1,0) + LE.Sig(mu,nu,sig,1,1,0);
    ans += RTP( tmp , tmp1 );

    tmp  = Hermite(V(x,nu)) * S(x,mu,nu,-1);
    tmp += Hermite(S(x,nu,mu,-1)) * V(x,mu);
    tmp1 = LE.Sig(mu,nu,rho,-1,1,0) + LE.Sig(mu,nu,sig,-1,1,0);
    ans += RTP( tmp , tmp1 );
  }
  return -32. * ans;
}

double W19_local( const Gauge& V , const SingleStaple& S , const LocalElements& LE , const Site& x , double _Complex *Uptr )
{
  Matrix<3,3> tmp( Uptr );
  double ans = 0.;
  for( const auto& m : HPE_comb::C22 ){
    const int mu= m[0] , nu=m[1] , rho=m[2] , sig=m[3];
    tmp = LE.P(rho,sig,0,0) + LE.P(sig,rho,0,0);
    ans += RTP( LE.P(mu,nu,0,0) , tmp );
  }
  return 16. * ans;
}

double W22_local( const Gauge& V , const SingleStaple& S , const LocalElements& LE , const Site& x , double _Complex *Uptr )
{
  double ans = 0.;
  for( const auto& m : HPE_comb::C22 ){
    const int mu= m[0] , nu=m[1] , rho=m[2] , sig=m[3];
    ans += RTP( LE.Sig(nu,mu,sig,1,-1,0) , LE.Sig(mu,nu,rho,1,-1,0) );
    ans += RTP( LE.Sig(nu,mu,sig,-1,-1,0) , LE.Sig(mu,nu,rho,1,1,0) );
    ans += RTP( LE.Sig(nu,mu,sig,-1,1,0) , LE.Sig(mu,nu,rho,-1,1,0) );
    ans += RTP( LE.Sig(nu,mu,sig,1,1,0) , LE.Sig(mu,nu,rho,-1,-1,0) );
  }
 return -32. * ans;
}

double W23_local( const Gauge& V , const SingleStaple& S , const LocalElements& LE , const Site& x , double _Complex *Uptr )
{
  double ans = 0.;
  for( const auto& m : HPE_comb::C2P2 ){
    const int mu= m[0] , nu=m[1] , rho=m[2] , sig=m[3];
    ans += RTP( LE.Sig(mu,nu,rho,1,1,0) , LE.Sig(mu,nu,sig,-1,-1,0) );
    ans += RTP( LE.Sig(mu,nu,rho,-1,1,0) , LE.Sig(mu,nu,sig,1,-1,0) );
  }
  return -32. * ans;
}

double W24_local( const Gauge& V , const SingleStaple& S , const LocalElements& LE , const Site& x , double _Complex *Uptr )
{
  Matrix<3,3> tmp( Uptr );
  double ans = 0.;
  for( const auto& m : HPE_comb::C2P2 ){
    const int mu= m[2] , nu=m[1] , rho=m[0] , sig=m[3];
    Site xm = x + mu;
    Site xn = x + nu;
    tmp = LE.Sig(mu,nu,rho,-1,-1,0) * S(x,mu,sig,0);
    tmp += Hermite(S(xn,mu,sig,0)) * LE.Sig(mu,nu,rho,1,-1,0);
    ans += RTP( tmp , V(xm,nu) );
    tmp  = LE.Sig(mu,nu,rho,1,1,0) * Hermite(S(xn,mu,sig,0));
    tmp += S(x,mu,sig,0) * LE.Sig(mu,nu,rho,-1,1,0);
    ans += RTP( tmp , Hermite(V(x,nu)) );
  }
  return -16. * ans;
}

double W25_local( const Gauge& V , const SingleStaple& S , const LocalElements& LE , const Site& x , double _Complex *Uptr )
{
  Matrix<3,3> tmp( Uptr ) , tmp1( Uptr+9 );
  double ans = 0.;
  for( const auto& m : HPE_comb::P4R ){
    const int mu= m[0] , nu=m[1] , rho=m[2] , sig=m[3];
    Site xn = x+nu;
    //(mu+)(nu+)(rho+)(sig+)(mu-)(nu-)(rho-)(sig-)
    tmp = S(xn,rho,sig,-1) * Hermite(S(x+rho,nu,mu,-1));
    tmp = tmp * Hermite(S(x,rho,mu,-1));
    Site xr = x-rho;
    //(mu+)(nu+)(rho-)(sig+)(mu-)(nu-)(rho+)(sig-)
    tmp1 = Hermite(S(xn-rho,rho,sig,-1)) * Hermite(S(xr,nu,mu,-1));
    tmp += tmp1 * S(xr,rho,mu,-1);
    tmp = tmp * LE.P(mu,sig,-1,-1);
    ans += RTP( tmp , S(x,nu,sig,-1) );
  }
  return -16. * ans;
}

double NNLO_W( const Gauge& V ){
  using Wfunc = double(*)( const Gauge& , const SingleStaple& , const LocalElements& , const Site& , double _Complex* );
  Wfunc Wlist[] = { 
    W1_11_12_13_20_local ,
    //W2_local , W3_local , W4_local , W5_local ,
    W2_5_local , W3_local , W4_local ,
    W6_7_local , W8_local , W9_local , W14_local ,
    W15_21_local , W16_17_local , W18_local , 
    W19_local , W22_local , W23_local , W24_local ,W25_local
  };
  const SingleStaple& S(V);  
  double ans = 0.;
#pragma omp parallel reduction( + : ans )
  {
    double _Complex ptr[27];
    Site x(V.lattice());
    ForAllSitesOMP(x){
      LocalElements LE( V , S , x );
      for( auto& func : Wlist )
        ans += func( V , S , LE , x , ptr );
    }
  }
  mpi.add(ans);
  return -2. * ans / V.lattice().size();
}


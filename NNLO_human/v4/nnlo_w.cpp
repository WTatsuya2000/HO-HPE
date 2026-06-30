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
  static constexpr int C2C2[6][4] = {
    {0,1,2,3} , {0,2,1,3} , {0,3,1,2} , {1,2,0,3} , {1,3,0,2} , {2,3,0,1} 
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

double W1_local( const Gauge& V , const SingleStaple& S , const LocalElements& LE , const Site& x , double _Complex *Uptr )
{
  Matrix<3,3> tmp( Uptr );
  double ans = 0.;
  for( const auto& m : HPE_comb::P2 ){
    const int mu=m[0] , nu=m[1];
    Site xm = x+mu;
    Site xn = x+nu;
    tmp = Hermite( S(x,nu,mu,-1) ) * V(x,mu);
    tmp = tmp * S(xm,nu,mu,1);
    tmp = tmp * Hermite( V(xn,mu) );
    ans += __real__( Trace(tmp) );
  }
  return -128. * ans;
}

double W2_local( const Gauge& V , const SingleStaple& S , const LocalElements& LE , const Site& x , double _Complex *Uptr )
{
  Matrix<3,3> tmp( Uptr );
  double ans = 0.;
  for( const auto& m : HPE_comb::C2 ){
    const int mu=m[0] , nu=m[1];
    tmp = LE.P(mu,nu,-1,1) + LE.P(mu,nu,1,-1);
    tmp = LE.P(mu,nu,1,1) * tmp;
    ans += __real__( Trace(tmp) );
  }
  return 32. * ans;
}

double W3_local( const Gauge& V , const SingleStaple& S , const LocalElements& LE , const Site& x , double _Complex *Uptr )
{
  Matrix<3,3> tmp( Uptr );
  double ans = 0.;
  for( const auto& m : HPE_comb::C2 ){
    const int mu=m[0] , nu=m[1];
    Site xn = x-nu;
    tmp = LE.P(mu,nu,1,1) * LE.P(mu,nu,-1,-1);
    ans += __real__( Trace(tmp) );
    tmp = LE.P(mu,nu,-1,1) * LE.P(mu,nu,1,-1);
    ans += __real__( Trace(tmp) );
  }
  return -32. * ans;
}

double W4_local( const Gauge& V , const SingleStaple& S , const LocalElements& LE , const Site& x , double _Complex *Uptr )
{
  Matrix<3,3> tmp( Uptr );
  double ans = 0.;
  for( const auto& m : HPE_comb::C2 ){
    const int mu=m[0] , nu=m[1];
    Site xn = x-nu;
    tmp = LE.P(mu,nu,1,1) * LE.P(nu,mu,-1,-1);
    ans += __real__( Trace(tmp) );
    tmp = LE.P(mu,nu,-1,1) * LE.P(nu,mu,-1,1);
    ans += __real__( Trace(tmp) );
  }
  return 64. * ans;
}

double W5_local( const Gauge& V , const SingleStaple& S , const LocalElements& LE , const Site& x , double _Complex *Uptr )
{
  Matrix<3,3> tmp( Uptr );
  double ans = 0.;
  for( const auto& m : HPE_comb::C2 ){
    const int mu=m[0] , nu=m[1];
    tmp = LE.P(mu,nu,1,1) * LE.P(mu,nu,1,1);
    ans += __real__( Trace(tmp) );
  }
  return 32. * ans / 2.;
}

double W6_local( const Gauge& V , const SingleStaple& S , const LocalElements& LE , const Site& x , double _Complex *Uptr )
{
  Matrix<3,3> tmp( Uptr ) , tmp1( Uptr+9 );
  double ans = 0.;
  for( const auto& m : HPE_comb::C2 ){
    const int mu=m[0] , nu=m[1];
    for( int s : {-1,1} ){
      tmp1 = LE.P(mu,nu,1,s) * LE.P(nu,mu,s,-1);
      tmp = LE.P(mu,nu,-1,-s) + LE.P(nu,mu,-s,1);
      tmp = tmp * tmp1;
      ans += __real__( Trace(tmp) );
    }
  }
  return -64. * ans;
}

double W7_local( const Gauge& V , const SingleStaple& S , const LocalElements& LE , const Site& x , double _Complex *Uptr )
{
  Matrix<3,3> tmp( Uptr );
  double ans = 0.;
  for( const auto& m : HPE_comb::C2 ){
    const int mu=m[0] , nu=m[1];
    Site xm = x-mu;
    tmp = S(x,mu,nu,-1) * Hermite( S(x,mu,nu,1) );
    tmp = tmp * Hermite( S(xm,mu,nu,1) );
    tmp = tmp * S(xm,mu,nu,-1);
    ans += __real__( Trace(tmp) );
  }
  return -128. * ans;
}

double W8_local( const Gauge& V , const SingleStaple& S , const LocalElements& LE , const Site& x , double _Complex *Uptr )
{
  Matrix<3,3> tmp( Uptr ) , tmp1( Uptr+9);
  double ans = 0.;
  for( const auto& m : HPE_comb::C1C2 ){
    const int mu=m[0] , nu=m[1] , rho=m[2];
    tmp  = LE.P(mu,nu,1,1) + LE.P(mu,nu,1,-1);
    tmp1 = LE.P(mu,rho,1,1) + LE.P(mu,rho,1,-1);
    tmp = tmp * tmp1;
    ans += __real__( Trace(tmp) );
  }
  return 32. * ans;
}

double W9_local( const Gauge& V , const SingleStaple& S , const LocalElements& LE , const Site& x , double _Complex *Uptr )
{
  Matrix<3,3> tmp( Uptr ) , tmp1( Uptr+9 );
  double ans = 0.;
  for( const auto& m : HPE_comb::P3 ){
    const int mu=m[0] , nu=m[1] , rho=m[2];
    tmp  = LE.P(mu,nu,1,1) + LE.P(mu,nu,1,-1);
    tmp1 = LE.P(rho,mu,1,-1) + LE.P(rho,mu,-1,-1);
    tmp  = tmp * tmp1;
    ans += __real__( Trace(tmp) );
  }
  return 32. * ans;
}

double W10_local( const Gauge& V , const SingleStaple& S , const LocalElements& LE , const Site& x , double _Complex *Uptr )
{
  return 0.;
}

double W11_local( const Gauge& V , const SingleStaple& S , const LocalElements& LE , const Site& x , double _Complex *Uptr )
{
  Matrix<3,3> tmp( Uptr );
  double ans = 0.;
  for( const auto& m : HPE_comb::P2 ){
    const int mu=m[0] , nu=m[1];
    const int
      rho = HPE_comb::remnant[mu][nu][0] ,
      sig = HPE_comb::remnant[mu][nu][1];
    Site xm = x + mu;
    Site xn = x + nu;
    tmp = S(x,nu,rho,0) + S(x,nu,sig,0);
    tmp = Hermite(tmp) * V(x,mu);
    tmp = tmp * S(xm,nu,mu,1);
    tmp = tmp * Hermite( V(xn,mu) );
    ans += __real__( Trace(tmp) );
    tmp = S(xm,nu,rho,0) + S(xm,nu,sig,0);
    tmp = tmp * Hermite( V(xn,mu) );
    tmp = tmp * Hermite( S(x,nu,mu,-1) );
    tmp = tmp * V(x,mu);
    ans += __real__( Trace(tmp) );
  }
  /*    
  for( const auto& m : HPE_comb::P3 ){
    const int mu=m[0] , nu=m[1] , rho=m[2];
    Site xm = x + mu;
    Site xn = x + nu;
    tmp = Hermite( S(x,nu,rho,0) ) * V(x,mu);
    tmp = tmp * S(xm,nu,mu,1);
    tmp = tmp * Hermite( V(xn,mu) );
    ans += __real__( Trace(tmp) );
    tmp = V(x,mu) * S(xm,nu,rho,0);
    tmp = tmp * Hermite( V(xn,mu) );
    tmp = tmp * Hermite( S(x,nu,mu,-1) );
    ans += __real__( Trace(tmp) );
  }
  */
  return -64. * ans;
}

double W12_local( const Gauge& V , const SingleStaple& S , const LocalElements& LE , const Site& x , double _Complex *Uptr )
{
  Matrix<3,3> tmp( Uptr );
  double ans = 0.;
  for( const auto& m : HPE_comb::C2C1 ){
    const int mu=m[0] , nu=m[1] , rho=m[2];
    Site xm = x + mu;
    Site xn = x + nu;
    for( int s : {-1,1} ){
      tmp = V(x,mu) * S(xm,nu,rho,s);
      tmp = tmp * Hermite( V(xn,mu) );
      tmp = tmp * Hermite( S(x,nu,rho,s) );
      ans += __real__( Trace(tmp) );
    }
  }
  return -32. * ans;
}

double W13_local( const Gauge& V , const SingleStaple& S , const LocalElements& LE , const Site& x , double _Complex *Uptr )
{
  Matrix<3,3> tmp( Uptr );
  double ans = 0.;
  for( const auto& m : HPE_comb::P3 ){
    const int mu=m[0] , nu=m[1] , rho=m[2];
    Site xm = x + mu;
    Site xn = x + nu;
    for( int s : {-1,1} ){
      tmp = V(x,mu) * S(xm,nu,rho,s);
      tmp = tmp * Hermite( V(xn,mu) );
      tmp = tmp * Hermite( S(x,nu,rho,-s) );
      ans += __real__( Trace(tmp) );
    }
  }
  return -32. * ans;
}

double W14_local( const Gauge& V , const SingleStaple& S , const LocalElements& LE , const Site& x , double _Complex *Uptr )
{
  Matrix<3,3> tmp( Uptr );
  double ans = 0.;
  for( const auto& m : HPE_comb::P3 ){
    const int mu=m[0] , nu=m[1] , rho=m[2];
    Site xm = x - mu;
    tmp = LE.P(mu,nu,1,1) * Hermite( S(xm,mu,nu,1) );
    tmp += LE.P(mu,nu,1,-1) * Hermite( S(xm,mu,nu,-1) );
    tmp = tmp * S(xm,mu,rho,0);
    ans += __real__( Trace(tmp) );
    tmp = LE.P(mu,nu,-1,1) * S(x,mu,nu,1);
    tmp += LE.P(mu,nu,-1,-1) * S(x,mu,nu,-1);
    tmp = tmp * Hermite( S(x,mu,rho,0) );
    ans += __real__( Trace(tmp) );
  }
  return -32. * ans;
}

double W15_local( const Gauge& V , const SingleStaple& S , const LocalElements& LE , const Site& x , double _Complex *Uptr )
{
  Matrix<3,3> tmp( Uptr );
  double ans = 0.;
  for( const auto& m : HPE_comb::C2C1 ){
    const int mu=m[0] , nu=m[1] , rho=m[2];
    for( int sm : {-1,1} ){
      for( int sn : {-1,1} ){
	tmp = LE.P(nu,rho,sn,-1) * LE.P(rho,mu,1,sm);
	tmp += LE.P(nu,rho,sn,1) * LE.P(rho,mu,-1,sm);
	tmp = tmp * LE.P(mu,nu,sm,sn);
	ans += __real__( Trace(tmp) );
      }
    }
  }
  return -32. * ans;
}

double W16_local( const Gauge& V , const SingleStaple& S , const LocalElements& LE , const Site& x , double _Complex *Uptr )
{
  Matrix<3,3> tmp( Uptr );
  double ans = 0.;
  for( const auto& m : HPE_comb::C2P2 ){
    const int mu=m[0] , nu=m[1] , rho=m[2];
    tmp = LE.Sig(mu,nu,rho, 1, 1,1) * LE.Sig(nu,mu,rho,-1,-1,-1);
    tmp+= LE.Sig(mu,nu,rho,-1, 1,1) * LE.Sig(nu,mu,rho,-1, 1,-1);
    tmp+= LE.Sig(mu,nu,rho,-1,-1,1) * LE.Sig(nu,mu,rho, 1, 1,-1);
    tmp+= LE.Sig(mu,nu,rho, 1,-1,1) * LE.Sig(nu,mu,rho, 1,-1,-1);
    ans += __real__( Trace(tmp) );
  }
  return -64. * ans;
}

double W17_local( const Gauge& V , const SingleStaple& S , const LocalElements& LE , const Site& x , double _Complex *Uptr )
{
  Matrix<3,3> tmp( Uptr );
  double ans = 0.;
  for( const auto& m : HPE_comb::C2P2 ){
    const int mu=m[0] , nu=m[1] , rho=m[2];
    tmp = LE.Sig(mu,nu,rho, 1, 1,1) * LE.Sig(mu,nu,rho,-1,-1,-1);
    tmp+= LE.Sig(mu,nu,rho, 1,-1,1) * LE.Sig(mu,nu,rho,-1, 1,-1);
    tmp+= LE.Sig(mu,nu,rho,-1, 1,1) * LE.Sig(mu,nu,rho, 1,-1,-1);
    tmp+= LE.Sig(mu,nu,rho,-1,-1,1) * LE.Sig(mu,nu,rho, 1, 1,-1);
    ans += -64.* __real__( Trace(tmp) );
  }
  return ans;
}

double W18_local( const Gauge& V , const SingleStaple& S , const LocalElements& LE , const Site& x , double _Complex *Uptr )
{
  Matrix<3,3> tmp( Uptr ) , tmp1( Uptr+9 );
  double ans = 0.;
  for( const auto& m : HPE_comb::C2C2 ){
    const int mu= m[0] , nu=m[1] , rho=m[2] , sig=m[3];
    Site xm = x + mu;
    Site xn = x + nu;
    tmp  = S(x,mu,nu,-1) * V(xm,nu);
    tmp += V(x,mu) * S(xm,nu,mu,1);
    tmp1 = LE.Sig(mu,nu,rho,-1,-1,1) + LE.Sig(mu,nu,rho,-1,-1,-1);
    tmp1+= LE.Sig(mu,nu,sig,-1,-1,1);
    tmp1+= LE.Sig(mu,nu,sig,-1,-1,-1);
    tmp  = tmp * tmp1;
    ans += __real__( Trace(tmp) );

    tmp  = S(xm,nu,mu,1) * Hermite(V(xn,mu));
    tmp += V(xm,nu) * Hermite(S(xn,mu,nu,1));
    tmp1 = LE.Sig(mu,nu,rho,1,-1,1) + LE.Sig(mu,nu,rho,1,-1,-1);
    tmp1+= LE.Sig(mu,nu,sig,1,-1,1);
    tmp1+= LE.Sig(mu,nu,sig,1,-1,-1);
    tmp  = tmp * tmp1;
    ans += __real__( Trace(tmp) );

    tmp  = Hermite(S(xn,mu,nu,1)) * Hermite(V(x,nu));
    tmp += Hermite(V(xn,mu)) * Hermite(S(x,nu,mu,-1));
    tmp1 = LE.Sig(mu,nu,rho,1,1,1) + LE.Sig(mu,nu,rho,1,1,-1);
    tmp1+= LE.Sig(mu,nu,sig,1,1,1);
    tmp1+= LE.Sig(mu,nu,sig,1,1,-1);
    tmp  = tmp * tmp1;
    ans += __real__( Trace(tmp) );

    tmp  = Hermite(V(x,nu)) * S(x,mu,nu,-1);
    tmp += Hermite(S(x,nu,mu,-1)) * V(x,mu);
    tmp1 = LE.Sig(mu,nu,rho,-1,1,1) + LE.Sig(mu,nu,rho,-1,1,-1);
    tmp1+= LE.Sig(mu,nu,sig,-1,1,1);
    tmp1+= LE.Sig(mu,nu,sig,-1,1,-1);
    tmp  = tmp * tmp1;
    ans += __real__( Trace(tmp) );
  }
  return -32. * ans;
}

double W19_local( const Gauge& V , const SingleStaple& S , const LocalElements& LE , const Site& x , double _Complex *Uptr )
{
  Matrix<3,3> tmp( Uptr ) , tmp1( Uptr+9 );
  double ans = 0.;
  for( const auto& m : HPE_comb::C22 ){
    const int mu= m[0] , nu=m[1] , rho=m[2] , sig=m[3];
    Site xm = x - mu;
    Site xr = x - rho;
    tmp  = LE.P(mu,nu,1,1) + LE.P(mu,nu,1,-1);
    tmp += LE.P(mu,nu,-1,1);
    tmp += LE.P(mu,nu,-1,-1);
    tmp1 = LE.P(rho,sig,1,1) + LE.P(rho,sig,1,-1);
    tmp1+= LE.P(rho,sig,-1,1);
    tmp1+= LE.P(rho,sig,-1,-1);
    tmp1 += Hermite(tmp1);
    tmp  = tmp * tmp1;
    ans += __real__( Trace(tmp) );
  }
  return 16. * ans;
}

double W20_local( const Gauge& V , const SingleStaple& S , const LocalElements& LE , const Site& x , double _Complex *Uptr )
{
  Matrix<3,3> tmp( Uptr );
  double ans = 0.;
  for( const auto& m : HPE_comb::P4 ){
    const int mu= m[0] , nu=m[1] , rho=m[2] , sig=m[3];
  Site xm = x + mu;
  Site xn = x + nu;
  tmp = V(x,mu) * S(xm,nu,rho,0);
  tmp = tmp * Hermite( V(xn,mu) );
  tmp = tmp * Hermite( S(x,nu,sig,0));
  ans += __real__( Trace(tmp) );
  }
  return -32. * ans;
}

double W21_local( const Gauge& V , const SingleStaple& S , const LocalElements& LE , const Site& x , double _Complex *Uptr )
{
  Matrix<3,3> tmp( Uptr ) , tmp1( Uptr+9 );
  double ans = 0.;
  for( int mu=0 ; mu<4 ; mu++ ){
    for( int nu=mu+1 ; nu<4 ; nu++ ){
      const int
	rho = HPE_comb::remnant[mu][nu][0] ,
	sig = HPE_comb::remnant[mu][nu][1];
      Site xm = x + mu;
      Site xn = x + nu;
      tmp  = Hermite(S(xn,mu,rho,0)) * Hermite(S(x,nu,sig,0));
      tmp += Hermite(S(xn,mu,sig,0)) * Hermite(S(x,nu,rho,0));
      tmp  = tmp * V(x,mu);
      tmp  = tmp * V(xm,nu);
      ans += __real__( Trace(tmp) );
      tmp  = Hermite(S(x,nu,rho,0)) * S(x,mu,sig,0);
      tmp += Hermite(S(x,nu,sig,0)) * S(x,mu,rho,0);
      tmp  = tmp * V(xm,nu);
      tmp  = tmp * Hermite(V(xn,mu));
      ans += __real__( Trace(tmp) );
      tmp  = S(x,mu,rho,0) * S(xm,nu,sig,0);
      tmp += S(x,mu,sig,0) * S(xm,nu,rho,0);
      tmp = tmp * Hermite(V(xn,mu));
      tmp = tmp * Hermite(V(x,nu));
      ans += __real__( Trace(tmp) );
      tmp = S(xm,nu,rho,0) * Hermite(S(xn,mu,sig,0));
      tmp += S(xm,nu,sig,0) * Hermite(S(xn,mu,rho,0));
      tmp = tmp * Hermite(V(x,nu));
      tmp = tmp * V(x,mu);
      ans += __real__( Trace(tmp) );
    }
  }
  return -16. * ans;
}

double W22_local( const Gauge& V , const SingleStaple& S , const LocalElements& LE , const Site& x , double _Complex *Uptr )
{
  Matrix<3,3> tmp( Uptr ) , tmp1( Uptr+9 );
  double ans = 0.;
  for( const auto& m : HPE_comb::C22 ){
    const int mu= m[0] , nu=m[1] , rho=m[2] , sig=m[3];
  Site xm = x - mu;
  Site xn = x - nu;
  tmp1 = LE.Sig(mu,nu,rho,1,-1,1) + LE.Sig(mu,nu,rho,1,-1,-1);
  tmp = LE.Sig(nu,mu,sig,1,-1,1) + LE.Sig(nu,mu,sig,1,-1,-1);
  tmp = tmp * tmp1;
  ans += __real__( Trace(tmp) );
  tmp1 = LE.Sig(mu,nu,rho,1,1,1) + LE.Sig(mu,nu,rho,1,1,-1);
  tmp = LE.Sig(nu,mu,sig,-1,-1,1) + LE.Sig(nu,mu,sig,-1,-1,-1);
  tmp = tmp * tmp1;
  ans += __real__( Trace(tmp) );
  tmp1 = LE.Sig(mu,nu,rho,-1,1,1) + LE.Sig(mu,nu,rho,-1,1,-1);
  tmp = LE.Sig(nu,mu,sig,-1,1,1) + LE.Sig(nu,mu,sig,-1,1,-1);
  tmp = tmp * tmp1;
  ans += __real__( Trace(tmp) );
  tmp1 = LE.Sig(mu,nu,rho,-1,-1,1) + LE.Sig(mu,nu,rho,-1,-1,-1);
  tmp = LE.Sig(nu,mu,sig,1,1,1) + LE.Sig(nu,mu,sig,1,1,-1);
  tmp = tmp * tmp1;
  ans += __real__( Trace(tmp) );
  }
 return -32. * ans;
}

double W23_local( const Gauge& V , const SingleStaple& S , const LocalElements& LE , const Site& x , double _Complex *Uptr )
{
  Matrix<3,3> tmp( Uptr ) , tmp1( Uptr+9 );
  double ans = 0.;
  for( const auto& m : HPE_comb::C2P2 ){
    const int mu= m[0] , nu=m[1] , rho=m[2] , sig=m[3];
  Site xm = x + mu;
  Site xn = x + nu;
  tmp1 = LE.Sig(mu,nu,rho,1,1,1) + LE.Sig(mu,nu,rho,1,1,-1);
  tmp = LE.Sig(mu,nu,sig,-1,-1,1) + LE.Sig(mu,nu,sig,-1,-1,-1);
  tmp = tmp * tmp1;
  ans += __real__( Trace(tmp) );
  tmp1 = LE.Sig(mu,nu,rho,-1,1,1) + LE.Sig(mu,nu,rho,-1,1,-1);
  tmp = LE.Sig(mu,nu,sig,1,-1,1) + LE.Sig(mu,nu,sig,1,-1,-1);
  tmp = tmp * tmp1;
  ans += __real__( Trace(tmp) );
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
    tmp = LE.Sig(mu,nu,rho,-1,-1,1) + LE.Sig(mu,nu,rho,-1,-1,-1);
    tmp = tmp * S(x,mu,sig,0);
    tmp = tmp * V(xm,nu);
    ans += __real__( Trace(tmp) );

    tmp = LE.Sig(mu,nu,rho,1,-1,1) + LE.Sig(mu,nu,rho,1,-1,-1);
    tmp = tmp * V(xm,nu);
    tmp = tmp * Hermite( S(xn,mu,sig,0) );
    ans += __real__( Trace(tmp) );

    tmp = LE.Sig(mu,nu,rho,1,1,1) + LE.Sig(mu,nu,rho,1,1,-1);
    tmp  = tmp * Hermite( S(xn,mu,sig,0) );
    tmp  = tmp * Hermite( V(x,nu) );
    ans += __real__( Trace(tmp) );

    tmp = LE.Sig(mu,nu,rho,-1,1,1) + LE.Sig(mu,nu,rho,-1,1,-1);
    tmp  = tmp * Hermite( V(x,nu) );
    tmp  = tmp * S(x,mu,sig,0);
    ans += __real__( Trace(tmp) );
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
    tmp = tmp * S(x,nu,sig,-1);
    ans += __real__( Trace(tmp) );
  }
  return -16. * ans;
}

double W_local( const Gauge& V , const SingleStaple& S , const LocalElements& LE , const Site& x , double _Complex *Uptr , int i ){
  using Wfunc = double(*)( const Gauge& , const SingleStaple& , const LocalElements& , const Site& , double _Complex* );
  Wfunc Wlist[] = { 
    W1_local , W2_local , W3_local , W4_local , W5_local ,
    W6_local , W7_local , W8_local , W9_local , W10_local ,
    W11_local , W12_local , W13_local , W14_local , W15_local ,
    W16_local , W17_local , W18_local , W19_local , W20_local ,
    W21_local , W22_local , W23_local , W24_local , W25_local 
  };
  return Wlist[i]( V , S , LE , x , Uptr );
}

double NNLO_W( const Gauge& V , const SingleStaple& S , int i ){
  double ans = 0.;
#pragma omp parallel reduction( + : ans )
  {
    double _Complex ptr[27];
    Site x(V.lattice());
    ForAllSitesOMP(x){
      LocalElements LE( V , S , x );
      ans += W_local( V , S , LE , x , ptr , i );
    }
  }
  mpi.add(ans);
  return -2. * ans / V.lattice().size();
}

double NNLO_W( const Gauge& V , const SingleStaple& S ){
  using Wfunc = double(*)( const Gauge& , const SingleStaple& , const LocalElements& , const Site& , double _Complex* );
  Wfunc Wlist[] = { 
    W1_local , W2_local , W3_local , W4_local , W5_local ,
    W6_local , W7_local , W8_local , W9_local , 
    W11_local , W12_local , W13_local , W14_local , W15_local ,
    W16_local , W17_local , W18_local , W19_local , W20_local ,
    W21_local , W22_local , W23_local , W24_local , W25_local 
  };
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


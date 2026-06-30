double _Complex P1_local( const Gauge& V , const SingleStaple& S , const Site& x , double _Complex *Uptr )
{
  Matrix<3,3> tmp( Uptr );
  tmp = V(x,0);
  auto xt = x;
  for( int i=1 ; i<4 ; i++ ){
    xt = xt+0;
    tmp = tmp * V(xt,0);
  }
  tmp = tmp * tmp;
  return -512./2. * Trace(tmp);
  //This function is evaluated only at tau=0 timeslice. Degeneracy factor 2 for two loops
}

double _Complex P2_local( const Gauge& V , const SingleStaple& S , const Site& x , double _Complex *Uptr )
{
  Matrix<3,3> tmp( Uptr );
  tmp = 0.;
  for( int mu=1 ; mu<4 ; mu++ )
    tmp += Hermite(S(x,0,mu,0));
  tmp = V(x,0) * tmp;
  auto xt = x;
  tmp = tmp * V(x,0);
  for( int i=1 ; i<4 ; i++ ){
    xt = xt+0;
    tmp = tmp * V(xt,0);
  }
  return -128. * Trace(tmp);
}

double _Complex P3_local( const Gauge& V , const SingleStaple& S , const Site& x , double _Complex *Uptr )
{
  Matrix<3,3> tmp( Uptr ) , tmp1( Uptr+9 );
  tmp = 0.;
  for( int mu=1 ; mu<4 ; mu++ ){
    auto xm = x - mu;
    tmp1 = V(x,mu) * S(x+mu,0,mu,1);
    tmp += tmp1 * Hermite( V(x+0,mu) );
    tmp1 = Hermite( V(xm,mu) ) * S(xm,0,mu,-1);
    tmp += tmp1 * V(xm+0,mu);
  }
  auto xt = x;
  for( int i=1 ; i<4 ; i++ ){
    xt = xt+0;
    tmp = tmp * V(xt,0);
  }
  return 128. * Trace(tmp);
}

double _Complex P4_local( const Gauge& V , const SingleStaple& S , const Site& x , double _Complex *Uptr )
{
  auto xt = x+0;
  auto xt2 = xt+0;
  Matrix<3,3> tmp( Uptr ) , tmp1( Uptr+9 );
  tmp = 0.;
  for( int mu=1 ; mu<4 ; mu++ ){
    tmp1 = V(x,mu) * S(x+mu,0,mu,1);
    tmp1 = tmp1 * S(xt+mu,0,mu,1);
    tmp += tmp1 * Hermite( V(xt2,mu) );
  }
  tmp = tmp * V(xt2,0);
  tmp = tmp * V(x-0,0);
  return 128. * Trace(tmp);
}

double _Complex P5_local( const Gauge& V , const SingleStaple& S , const Site& x , double _Complex *Uptr )
{
  auto xt = x+0;
  auto xt2 = xt+0;
  Matrix<3,3> tmp( Uptr ) , tmp1( Uptr+9 );
  tmp = 0.;
  for( int mu=1 ; mu<4 ; mu++ ){
    tmp1 = V(x,mu) * V(x+mu,0);
    tmp1 = tmp1 * S(xt+mu,0,mu,1);
    tmp += tmp1 * Hermite( V(xt2,mu) );
    auto xm = x-mu;
    tmp1 = Hermite( V(xm,mu) ) * V(xm,0);
    tmp1 = tmp1 * S(xm+0,0,mu,-1);
    tmp += tmp1 * V(xt2-mu,mu);
  }
  for( int i=2 ; i<4 ; i++ ){
    xt = xt+0;
    tmp = tmp * V(xt,0);
  }
  return 64. * Trace(tmp);
}

double _Complex P6_local( const Gauge& V , const SingleStaple& S , const Site& x , double _Complex *Uptr )
{
  auto xt = x+0;
  auto xt2 = xt+0;
  Matrix<3,3> tmp( Uptr ) , tmp1( Uptr+9 );
  tmp = 0.;
  for( int mu=1 ; mu<4 ; mu++ ){
    tmp1 = V(x,mu) * S(x+mu,0,mu,1);
    tmp1 = tmp1 * V(xt+mu,0);
    tmp += tmp1 * Hermite(V(xt2,mu));
    auto xm = x-mu;
    tmp1 = Hermite( V(xm,mu) ) * S(xm,0,mu,-1);
    tmp1 = tmp1 * V(xm+0,0);
    tmp += tmp1 * V(xt2-mu,mu);
  }
  for( int i=2 ; i<4 ; i++ ){
    xt = xt+0;
    tmp = tmp * V(xt,0);
  }
  return 64. * Trace(tmp);
}

double _Complex P4_5_6_local( const Gauge& V , const SingleStaple& S , const Site& x , double _Complex *Uptr )
{
  auto xt = x+0;
  auto xt2 = xt+0;
  Matrix<3,3> tmp( Uptr ) , tmp1( Uptr+9 );
  tmp = 0.;
  for( int mu=1 ; mu<4 ; mu++ ){
    auto xm = x + mu;
    auto xtm = xt + mu;
    tmp1 = S(xm,0,mu,1) * S(xtm,0,mu,1);
    tmp1 += V(xm,0) * S(xtm,0,mu,1);
    tmp1 += S(xm,0,mu,1) * V(xtm,0);
    tmp1 = V(x,mu) * tmp1;
    tmp += tmp1 * Hermite( V(xt2,mu) );
  }
  for( int mu=1 ; mu<4 ; mu++ ){
    auto xm = x - mu;
    auto xtm = xt - mu;
    tmp1 = S(xm,0,mu,-1) * S(xtm,0,mu,-1);
    tmp1 += V(xm,0) * S(xtm,0,mu,-1);
    tmp1 += S(xm,0,mu,-1) * V(xtm,0);
    tmp1 = Hermite( V(xm,mu) ) * tmp1;
    tmp += tmp1 * V(xt2-mu,mu);
  }

  tmp = tmp * V(xt2,0);
  tmp = tmp * V(x-0,0);
  return 64. * Trace(tmp);
}

double _Complex P7_local( const Gauge& V , const SingleStaple& S , const Site& x , double _Complex *Uptr )
{
  auto xt = x+0;
  Matrix<3,3> tmp( Uptr );
  tmp = 0.;
  for( int mu=1 ; mu<4 ; mu++ ){
    tmp += S(x,0,mu,1) * S(xt,0,mu,-1);
    tmp += S(x,0,mu,-1) * S(xt,0,mu,1);
  }
  for( int i=2 ; i<4 ; i++ ){
    xt = xt+0;
    tmp = tmp * V(xt,0);
  }
  return 64. * Trace(tmp);
}

double _Complex P8_local( const Gauge& V , const SingleStaple& S , const Site& x , double _Complex *Uptr )
{
  auto xt = x+0;
  auto xt2 = xt+0;
  Matrix<3,3> tmp( Uptr ) , tmp1( Uptr+9 );
  tmp = 0.;
  for( int mu=1 ; mu<4 ; mu++ ){
    tmp1 = S(x,0,mu,1) * V(xt,0);
    tmp += tmp1 * S(xt2,0,mu,-1);
  }
  tmp = tmp * V(x-0,0);
  return 32. * Trace(tmp);
}

double _Complex P9_local( const Gauge& V , const SingleStaple& S , const Site& x , double _Complex *Uptr )
{
  auto xt = x+0;
  auto xt2 = xt+0;
  Matrix<3,3> tmp( Uptr ) , tmp1( Uptr+9 );
  tmp = 0.;
  for( int mu=1 ; mu<4 ; mu++ ){
    tmp1 = S(x,0,mu,1) * V(xt,0);
    tmp += tmp1 * S(xt2,0,mu,1);
    tmp1 = S(x,0,mu,-1) * V(xt,0);
    tmp += tmp1 * S(xt2,0,mu,-1);
  }
  tmp = tmp * V(x-0,0);
  return 32. * Trace(tmp);
  //This function is evaluated only at tau=0 timeslice.
}

double _Complex P10_local( const Gauge& V , const SingleStaple& S , const Site& x , double _Complex *Uptr )
{
  Matrix<3,3> tmp( Uptr );
  tmp = 0.;
  for( int mu=1 ; mu<4 ; mu++ ){
    auto xm = x - mu;
    for( int nu=mu+1 ; nu<4 ; nu++ ){
      tmp += S(x,mu,nu,0) * Hermite( V(x,mu) );
      tmp += Hermite( S(xm,mu,nu,0) ) * V(xm,mu);
    }
  }
  tmp += Hermite(tmp);
  auto xt = x;
  tmp = tmp * V(x,0);
  for( int i=1 ; i<4 ; i++ ){
    xt = xt+0;
    tmp = tmp * V(xt,0);
  }
  return -64. * Trace(tmp);
}

double _Complex P11_local( const Gauge& V , const SingleStaple& S , const Site& x , double _Complex *Uptr )
{
  return 0.;
}
double _Complex P12_local( const Gauge& V , const SingleStaple& S , const Site& x , double _Complex *Uptr )
{
  return 0.;
}
double _Complex P13_local( const Gauge& V , const SingleStaple& S , const Site& x , double _Complex *Uptr )
{
  return 0.;
}
double _Complex P14_local( const Gauge& V , const SingleStaple& S , const Site& x , double _Complex *Uptr )
{
  return 0.;
}
double _Complex P15_local( const Gauge& V , const SingleStaple& S , const Site& x , double _Complex *Uptr )
{
  return 0.;
}

double _Complex P16_local( const Gauge& V , const SingleStaple& S , const Site& x , double _Complex *Uptr )
{
  Matrix<3,3> tmp( Uptr ) , tmp1( Uptr+9 );
  tmp = 0.;
  auto xt = x+0;
  auto xt2 = xt+0;
  for( int mu=1 ; mu<4 ; mu++ ){
    auto xm = x - mu;
    auto xtm = xt - mu;
    auto xt2m = xt2 - mu;
    for( int nu=1 ; nu<4 ; nu++ ){
      if( mu==nu ) continue;
      tmp1 = S(xt+mu,0,nu,1) * Hermite( S(xt2,mu,nu,1) );
      tmp1 += S(xt+mu,0,nu,-1) * Hermite( S(xt2,mu,nu,-1) );
      tmp1 = V(x+mu,0) * tmp1;
      tmp += V(x,mu) * tmp1;
      tmp1 = S(xtm,0,nu,1) * S(xt2m,mu,nu,1);
      tmp1 += S(xtm,0,nu,-1) * S(xt2m,mu,nu,-1);
      tmp1 = V(xm,0) * tmp1;
      tmp += Hermite( V(xm,mu) ) * tmp1;
    }
  }
  tmp = tmp * V(xt2,0);
  tmp = tmp * V(x-0,0);
  return 32. * Trace(tmp);
}

double _Complex P17_local( const Gauge& V , const SingleStaple& S , const Site& x , double _Complex *Uptr )
{
  Matrix<3,3> tmp( Uptr ) , tmp1( Uptr+9 );
  tmp = 0.;
  auto xt = x+0;
  auto xt2 = xt+0;
  for( int mu=1 ; mu<4 ; mu++ ){
    for( int nu=1 ; nu<4 ; nu++ ){
      if( mu==nu ) continue;
      auto xtn = xt - nu;
      auto xt2n = xt2 - nu;
      tmp1 = S(x,0,mu,1) * S(xt,nu,mu,1);
      tmp1 += S(x,0,mu,-1) * S(xt,nu,mu,-1);
      tmp1 = tmp1 * V(xt+nu,0);
      tmp += tmp1 * Hermite(V(xt2,nu));
      tmp1 = S(x,0,mu,1) * Hermite(S(xtn,nu,mu,1) );
      tmp1 += S(x,0,mu,-1) * Hermite(S(xtn,nu,mu,-1));
      tmp1 = tmp1 * V(xtn,0);
      tmp += tmp1 * V(xt2n,nu);
    }
  }
  tmp = tmp * V(xt2,0);
  tmp = tmp * V(x-0,0);
  return 32. * Trace(tmp);
}

double _Complex P18_local( const Gauge& V , const SingleStaple& S , const Site& x , double _Complex *Uptr )
{
  Matrix<3,3> tmp( Uptr ) , tmp1( Uptr+9 );
  tmp = 0.;
  auto xt = x+0;
  auto xt2 = xt+0;
  for( int mu=1 ; mu<4 ; mu++ ){
    for( int nu=1 ; nu<4 ; nu++ ){
      if( mu==nu ) continue;
      auto xn = x - nu;
      tmp1 = S(x,nu,mu,1) * S(x+nu,0,mu,1);
      tmp1 += S(x,nu,mu,-1) * S(x+nu,0,mu,-1);
      tmp1 = tmp1 * V(xt+nu,0);
      tmp += tmp1 * Hermite(V(xt2,nu));
      tmp1 = Hermite( S(xn,nu,mu,1) ) * S(xn,0,mu,1);
      tmp1 += Hermite( S(xn,nu,mu,-1) ) * S(xn,0,mu,-1);
      tmp1 = tmp1 * V(xt-nu,0);
      tmp += tmp1 * V(xt2-nu,nu);
    }
  }
  tmp = tmp * V(xt2,0);
  tmp = tmp * V(x-0,0);
  return 32. * Trace(tmp);
}

double _Complex P19_local( const Gauge& V , const SingleStaple& S , const Site& x , double _Complex *Uptr )
{
  Matrix<3,3> tmp( Uptr ) , tmp1( Uptr+9 );
  tmp = 0.;
  auto xt = x+0;
  auto xt2 = xt+0;
  auto xt3 = xt2+0;
  for( int mu=1 ; mu<4 ; mu++ ){
    for( int nu=1 ; nu<4 ; nu++ ){
      if( mu==nu ) continue;
      auto xtm = xt+mu;
      auto xt2n = xt2+nu;
      tmp1 = V(x,mu) * V(x+mu,0);
      tmp1 = tmp1 * V(xtm,nu);
      tmp1 = tmp1 * V(xtm+nu,0);
      tmp1 = tmp1 * Hermite(V(xt2n,mu));
      tmp1 = tmp1 * V(xt2n,0);
      tmp1 = tmp1 * Hermite( V(xt3,nu) );
      tmp += tmp1 * V(xt3,0);
    }
  }
  return 32. * Trace(tmp);
}

double _Complex P20_local( const Gauge& V , const SingleStaple& S , const Site& x , double _Complex *Uptr )
{
  Matrix<3,3> tmp( Uptr ) , tmp1( Uptr+9 );
  tmp = 0.;
  auto xt = x+0;
  for( int mu=1 ; mu<4 ; mu++ ){
    auto xm = x - mu;
    for( int nu=mu+1 ; nu<4 ; nu++ ){
      auto xn = x - nu;
      tmp1 = V(x,mu) * S(x+mu,0,nu,0);
      tmp += tmp1 * Hermite( V(xt,mu) );
      tmp1 = Hermite( V(xm,mu) ) * S(xm,0,nu,0);
      tmp += tmp1 * V(xt-mu,mu);
      tmp1 = V(x,nu) * S(x+nu,0,mu,0);
      tmp += tmp1 * Hermite( V(xt,nu) );
      tmp1 = Hermite( V(xn,nu) ) * S(xn,0,mu,0);
      tmp += tmp1 * V(xt-nu,nu);
    }
  }
  tmp = tmp * V(xt,0);
  tmp = tmp * V(xt+0,0);
  tmp = tmp * V(x-0,0);
  return 64. * Trace(tmp);
}

double _Complex P3_20_local( const Gauge& V , const SingleStaple& S , const Site& x , double _Complex *Uptr )
{
  Matrix<3,3> tmp( Uptr ) , tmp1( Uptr+9 );
  tmp = 0.;
  auto xt = x+0;
  for( int mu=1 ; mu<4 ; mu++ ){
    const int nu=mu%3+1 , rho=(mu+1)%3+1;
    auto xm = x + mu;
    tmp1 = S(xm,0,nu,0) + S(xm,0,rho,0);
    tmp1+= 2.* S(xm,0,mu,1);
    tmp1 = V(x,mu) * tmp1;
    tmp += tmp1 * Hermite( V(xt,mu) );
    xm = x - mu;
    tmp1 = S(xm,0,nu,0) + S(xm,0,rho,0);
    tmp1+= 2.* S(xm,0,mu,-1);
    tmp1 = Hermite(V(xm,mu)) * tmp1;
    tmp += tmp1 * V(xt-mu,mu);
  }
  tmp = tmp * V(xt,0);
  tmp = tmp * V(xt+0,0);
  tmp = tmp * V(x-0,0);
  return 64. * Trace(tmp);
}

double _Complex P21_local( const Gauge& V , const SingleStaple& S , const Site& x , double _Complex *Uptr )
{
  Matrix<3,3> tmp( Uptr ) , tmp1( Uptr+9 );
  tmp = 0.;
  auto xt = x+0;
  auto xt2 = xt+0;
  for( int mu=1 ; mu<4 ; mu++ ){
    auto xm = x - mu;
    auto xtm = xt - mu;
    auto xt2m = xt2 - mu;
    for( int nu=mu+1 ; nu<4 ; nu++ ){
      tmp1 = V(x,mu) * S(x+mu,0,nu,1);
      tmp1 = tmp1 * S(xt+mu,0,nu,1);
      tmp += tmp1 * Hermite( V(xt2,mu) );
      tmp1 = V(x,mu) * S(x+mu,0,nu,-1);
      tmp1 = tmp1 * S(xt+mu,0,nu,-1);
      tmp += tmp1 * Hermite( V(xt2,mu) );
      tmp1 = Hermite( V(xm,mu) ) * S(xm,0,nu,1);
      tmp1 = tmp1 * S(xtm,0,nu,1);
      tmp += tmp1 * V(xt2m,mu);
      tmp1 = Hermite( V(xm,mu) ) * S(xm,0,nu,-1);
      tmp1 = tmp1 * S(xtm,0,nu,-1);
      tmp += tmp1 * V(xt2m,mu);
    }
  }
  tmp = tmp * V(xt2,0);
  tmp = tmp * V(x-0,0);
  return 64. * Trace(tmp);
}

double _Complex P22_local( const Gauge& V , const SingleStaple& S , const Site& x , double _Complex *Uptr )
{
  Matrix<3,3> tmp( Uptr ) , tmp1( Uptr+9 );
  tmp = 0.;
  auto xt = x+0;
  auto xt2 = xt+0;
  for( int mu=1 ; mu<4 ; mu++ ){
    auto xm = x - mu;
    for( int nu=mu+1 ; nu<4 ; nu++ ){
      auto xn = x - nu;
      tmp1 = V(x,mu) * V(x+mu,0);
      tmp1 = tmp1 * S(xt+mu,0,nu,0);
      tmp += tmp1 * Hermite( V(xt2,mu) );
      tmp1 = Hermite( V(xm,mu) )* V(xm,0);
      tmp1 = tmp1 * S(xt-mu,0,nu,0);
      tmp += tmp1 * V(xt2-mu,mu);
      tmp1 = V(x,nu) * V(x+nu,0);
      tmp1 = tmp1 * S(xt+nu,0,mu,0);
      tmp += tmp1 * Hermite( V(xt2,nu) );
      tmp1 = Hermite( V(xn,nu) )* V(xn,0);
      tmp1 = tmp1 * S(xt-nu,0,mu,0);
      tmp += tmp1 * V(xt2-nu,nu);
    }
  }
  tmp = tmp * V(xt+0,0);
  tmp = tmp * V(x-0,0);
  return 32. * Trace(tmp);
}

double _Complex P23_local( const Gauge& V , const SingleStaple& S , const Site& x , double _Complex *Uptr )
{
  Matrix<3,3> tmp( Uptr ) , tmp1( Uptr+9 );
  tmp = 0.;
  auto xt = x+0;
  auto xt2 = xt+0;
  for( int mu=1 ; mu<4 ; mu++ ){
    auto xm = x - mu;
    for( int nu=mu+1 ; nu<4 ; nu++ ){
      auto xn = x - nu;
      tmp1 = V(x,mu) * S(x+mu,0,nu,0);
      tmp1 = tmp1 * V(xt+mu,0);
      tmp += tmp1 * Hermite( V(xt2,mu) );
      tmp1 = Hermite( V(xm,mu) )* S(xm,0,nu,0);
      tmp1 = tmp1 * V(xt-mu,0);
      tmp += tmp1 * V(xt2-mu,mu);
      tmp1 = V(x,nu) * S(x+nu,0,mu,0);
      tmp1 = tmp1 * V(xt+nu,0);
      tmp += tmp1 * Hermite( V(xt2,nu) );
      tmp1 = Hermite( V(xn,nu) )* S(xn,0,mu,0);
      tmp1 = tmp1 * V(xt-nu,0);
      tmp += tmp1 * V(xt2-nu,nu);
    }
  }
  tmp = tmp * V(xt+0,0);
  tmp = tmp * V(x-0,0);
  return 32. * Trace(tmp);
}

double _Complex P4_5_6_21_22_23local( const Gauge& V , const SingleStaple& S , const Site& x , double _Complex *Uptr )
{
  auto xt = x+0;
  auto xt2 = xt+0;
  Matrix<3,3> tmp( Uptr ) , tmp1( Uptr+9 ) , tmp2( Uptr+18 );
  tmp = 0.;
  for( int mu=1 ; mu<4 ; mu++ ){
    const int nu=mu%3+1 , rho=(mu+1)%3+1;
    auto xm = x + mu;
    auto xtm = xt + mu;
    tmp1 = S(xm,0,mu,1) * S(xtm,0,mu,1);
    tmp1 = 2.* tmp1;
    tmp1+= S(xm,0,nu,1) * S(xtm,0,nu,1);
    tmp1+= S(xm,0,nu,-1) * S(xtm,0,nu,-1);
    tmp1+= S(xm,0,rho,1) * S(xtm,0,rho,1);
    tmp1+= S(xm,0,rho,-1) * S(xtm,0,rho,-1);
    tmp2 = S(xtm,0,nu,0) + S(xtm,0,rho,0);
    tmp2+= 2.* S(xtm,0,mu,1);
    tmp1+= V(xm,0) * tmp2;
    tmp2 = S(xm,0,nu,0) + S(xm,0,rho,0);
    tmp2+= 2.* S(xm,0,mu,1);
    tmp1+= tmp2 * V(xtm,0);
    tmp1 = V(x,mu) * tmp1;
    tmp += tmp1 * Hermite( V(xt2,mu) );

    xm = x - mu;
    xtm = xt - mu;
    tmp1 = S(xm,0,mu,-1) * S(xtm,0,mu,-1);
    tmp1 = 2.* tmp1;
    tmp1+= S(xm,0,nu,1) * S(xtm,0,nu,1);
    tmp1+= S(xm,0,nu,-1) * S(xtm,0,nu,-1);
    tmp1+= S(xm,0,rho,1) * S(xtm,0,rho,1);
    tmp1+= S(xm,0,rho,-1) * S(xtm,0,rho,-1);
    tmp2 = S(xtm,0,nu,0) + S(xtm,0,rho,0);
    tmp2+= 2.* S(xtm,0,mu,-1);
    tmp1+= V(xm,0) * tmp2;
    tmp2 = S(xm,0,nu,0) + S(xm,0,rho,0);
    tmp2+= 2.* S(xm,0,mu,-1);
    tmp1+= tmp2 * V(xtm,0);
    tmp1 = Hermite( V(xm,mu) ) * tmp1;
    tmp += tmp1 * V(xt2-mu,mu);
  }
  tmp = tmp * V(xt2,0);
  tmp = tmp * V(x-0,0);
  return 32. * Trace(tmp);
}

double _Complex P24_local( const Gauge& V , const SingleStaple& S , const Site& x , double _Complex *Uptr )
{
  Matrix<3,3> tmp( Uptr );
  tmp = 0.;
  auto xt = x+0;
  for( int mu=1 ; mu<4 ; mu++ ){
    for( int nu=mu+1 ; nu<4 ; nu++ ){
      tmp += S(x,0,mu,0) * S(xt,0,nu,0);
      tmp += S(x,0,nu,0) * S(xt,0,mu,0);
    }
  }
  tmp = tmp * V(xt+0,0);
  tmp = tmp * V(x-0,0);
  return 32. * Trace(tmp);
}

double _Complex P7_24_local( const Gauge& V , const SingleStaple& S , const Site& x , double _Complex *Uptr )
{
  Matrix<3,3> tmp( Uptr ) , tmp1( Uptr+9 ) , tmp2( Uptr+18 );
  tmp = 0.;
  auto xt = x+0;
  for( int mu=1 ; mu<4 ; mu++ ){
    const int nu=mu%3+1 , rho=(mu+1)%3+1;
    tmp1 = S(xt,0,nu,0) + S(xt,0,rho,0);
    tmp2 = tmp1;
    tmp1+= 2.* S(xt,0,mu,-1);
    tmp += S(x,0,mu,1) * tmp1;
    tmp2+= 2.* S(xt,0,mu,1);
    tmp += S(x,0,mu,-1) * tmp2;
  }
  tmp = tmp * V(xt+0,0);
  tmp = tmp * V(x-0,0);
  return 32. * Trace(tmp);
}

double _Complex P25_local( const Gauge& V , const SingleStaple& S , const Site& x , double _Complex *Uptr )
{
  Matrix<3,3> tmp( Uptr ) , tmp1( Uptr+9 );
  tmp = 0.;
  auto xt = x+0;
  auto xt2 = xt+0;
  for( int mu=1 ; mu<4 ; mu++ ){
    for( int nu=mu+1 ; nu<4 ; nu++ ){
      tmp1 = S(x,0,mu,0) * V(xt,0);
      tmp += tmp1 * S(xt2,0,nu,0);
    }
  }
  tmp = tmp * V(x-0,0);
  return 32. * Trace(tmp);
}

double _Complex P4_local( const Gauge& V , const SingleStaple& S , const Site& x , double _Complex *Uptr , int i ){
  using Pfunc = double _Complex(*)( const Gauge& , const SingleStaple& , const Site& , double _Complex* );
  Pfunc Plist[] = { 
    P1_local , P2_local , P3_local , P4_local , P5_local ,
    P6_local , P7_local , P8_local , P9_local , P10_local ,
    P11_local , P12_local , P13_local , P14_local , P15_local ,
    P16_local , P17_local , P18_local , P19_local , P20_local ,
    P21_local , P22_local , P23_local , P24_local , P25_local ,
    P4_5_6_local , P3_20_local , P7_24_local , P4_5_6_21_22_23local
  };
  return Plist[i]( V , S , x , Uptr );
}

double _Complex NNLO_P4( const Gauge& V , const SingleStaple& S , int i ){
  double _Complex ans = 0.;
#pragma omp parallel reduction( + : ans )
  {
    double _Complex ptr[27];
    Site x(V.lattice());
    if( i==0 or i==8 ){
      ForAllSitesOMP(x){
	if( x(0) == 0 )
	ans += P4_local( V , S , x , ptr , i );
      }
    }
    else{ 
      ForAllSitesOMP(x){
	ans += P4_local( V , S , x , ptr , i );
      }
    }
  }
  double _Complex tmp;
  MPI_Allreduce( (void*)&ans , (void*)&tmp , 1 ,
		 MPI_DOUBLE_COMPLEX , MPI_SUM , MPI_COMM_WORLD );
  return 2. * tmp / V.lattice().size();
}

double _Complex NNLO_P4( const Gauge& V , const SingleStaple& S )
{
  using Pfunc = double _Complex(*)( const Gauge& , const SingleStaple& , const Site& , double _Complex* );
  Pfunc Plist_t_ne0[] = { 
    P2_local , P3_20_local , P4_5_6_21_22_23local ,
    P7_24_local , P8_local , P10_local ,
    P16_local , P17_local , P18_local , P19_local , 
    P25_local
  };
  Pfunc Plist_t0[] = { P1_local , P9_local };
  double _Complex ans = 0.;
#pragma omp parallel reduction( + : ans )
  {
    double _Complex ptr[27];
    Site x(V.lattice());
    ForAllSitesOMP(x){
      if( x(0) == 0 )
	for( auto func : Plist_t0 )
	  ans += func( V , S , x , ptr );
      for( auto func : Plist_t_ne0 )
	ans += func( V , S , x , ptr );
    }
  }
  double _Complex tmp;
  MPI_Allreduce( (void*)&ans , (void*)&tmp , 1 ,
		 MPI_DOUBLE_COMPLEX , MPI_SUM , MPI_COMM_WORLD );
  return 2. * tmp / V.lattice().size();
}

double _Complex NNLO_P4( const Gauge& V )
{
  using Pfunc = double _Complex(*)( const Gauge& , const SingleStaple& , const Site& , double _Complex* );
  Pfunc Plist_t_ne0[] = { 
    P2_local , P3_20_local , P4_5_6_21_22_23local ,
    P7_24_local , P8_local , P10_local ,
    P16_local , P17_local , P18_local , P19_local , 
    P25_local
  };
  Pfunc Plist_t0[] = { P1_local , P9_local };
  const SingleStaple& S(V);
  double _Complex ans = 0.;
#pragma omp parallel reduction( + : ans )
  {
    double _Complex ptr[27];
    Site x(V.lattice());
    ForAllSitesOMP(x){
      if( x(0) == 0 )
	for( auto func : Plist_t0 )
	  ans += func( V , S , x , ptr );
      for( auto func : Plist_t_ne0 )
	ans += func( V , S , x , ptr );
    }
  }
  double _Complex tmp;
  MPI_Allreduce( (void*)&ans , (void*)&tmp , 1 ,
		 MPI_DOUBLE_COMPLEX , MPI_SUM , MPI_COMM_WORLD );
  return 2. * tmp / V.lattice().size();
}

double _Complex NNLO_P4L1( const Gauge& V , const SingleStaple& S )
{
  using Pfunc = double _Complex(*)( const Gauge& , const SingleStaple& , const Site& , double _Complex* );
  Pfunc Plist_t_ne0[] = { 
    P2_local , P3_20_local , P4_5_6_21_22_23local ,
    P7_24_local , P8_local , P10_local ,
    P16_local , P17_local , P18_local , P19_local , 
    P25_local
  };
  Pfunc Plist_t0[] = { P9_local };
  double _Complex ans = 0.;
#pragma omp parallel reduction( + : ans )
  {
    double _Complex ptr[27];
    Site x(V.lattice());
    ForAllSitesOMP(x){
      if( x(0) == 0 )
	for( auto func : Plist_t0 )
	  ans += func( V , S , x , ptr );
      for( auto func : Plist_t_ne0 )
	ans += func( V , S , x , ptr );
    }
  }
  double _Complex tmp;
  MPI_Allreduce( (void*)&ans , (void*)&tmp , 1 ,
		 MPI_DOUBLE_COMPLEX , MPI_SUM , MPI_COMM_WORLD );
  return 2. * tmp / V.lattice().size();
}

double _Complex NNLO_P4L2( const Gauge& V , const SingleStaple& S )
{
  double _Complex ans = 0.;
#pragma omp parallel reduction( + : ans )
  {
    double _Complex ptr[27];
    Site x(V.lattice());
    ForAllSitesOMP(x){
      if( x(0) == 0 )
	ans += P1_local(V,S,x,ptr);
    }
  }
  double _Complex tmp;
  MPI_Allreduce( (void*)&ans , (void*)&tmp , 1 ,
		 MPI_DOUBLE_COMPLEX , MPI_SUM , MPI_COMM_WORLD );
  return 2. * tmp / V.lattice().size();
}


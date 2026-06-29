
void set_cold(gauge_field &U) {
  mdp << "Creating a cold gauge configuration" << '\n';
  site x(U.lattice());
  int mu;
  forallsites(x) 
    for( mu=0 ; mu<U.ndim ; mu++ ) 
      U(x,mu) = mdp_identity(U.nc);
  mpi.barrier();
  mdp << "communicating" << endl;
  mpi.barrier();
  U.update();
}

bool heatbath_SU2( mdp_prng &random,
		   double beta_eff , double _Complex *a )
{
  const double e0 = __real__(a[0]) + __real__(a[3]);
  const double e1 = __imag__(a[1]) + __imag__(a[2]);
  const double e2 = __real__(a[1]) - __real__(a[2]);
  const double e3 = __imag__(a[0]) - __imag__(a[3]);

#if defined(BLUEGENE)// this routine is little bit faster on bluegene
  const double dk = 1./sqrt(e0*e0+e1*e1+e2*e2+e3*e3);
  const double p0 = beta_eff/dk;
  double _Complex u0 = e0 - e3 * _Complex_I;
  double _Complex u2 = e2 - e1 * _Complex_I;
  double _Complex u1 =-e2 - e1 * _Complex_I;
  double _Complex u3 = e0 + e3 * _Complex_I;
  u0 *= dk; u1 *= dk; u2 *= dk; u3 *= dk;
#else
  const double dk = sqrt( e0*e0 + e1*e1 + e2*e2 + e3*e3 );
  const double p0 = dk * beta_eff;
  double _Complex u0 = e0 - e3 * _Complex_I;
  double _Complex u2 = e2 - e1 * _Complex_I;
  double _Complex u1 =-e2 - e1 * _Complex_I;
  double _Complex u3 = e0 + e3 * _Complex_I;
  u0 /= dk; u2 /= dk; u1 /= dk; u3 /= dk;
#endif

  int hit=0; // 5-hit ; by M.K. Apr.2008
  double r4 , delta;
  do{
    if( hit++ > 5 ) return false; // 5-hit ;
    const double r1 = random.plain();
    const double r2 = random.plain();
    const double r3 = cos( 2.*Pi*random.plain() );
    delta = ( log(r2) + r3*r3 * log(r1) ) / (-p0);
    r4 = random.plain();
  } while (r4*r4 > (1.0-0.5*delta));

  const double a0 = 1. - delta;
  const double cos_theta = 2. * random.plain() - 1.;
  const double sin_theta = sqrt( 1. - cos_theta*cos_theta );
  const double sin_alpha = sqrt( 1. - a0*a0 );

  const double a3 = sin_alpha*cos_theta;
  const double _Complex v0 = a0 + a3 * _Complex_I;
  const double _Complex v3 = a0 - a3 * _Complex_I;

  const double phi = 2. * Pi * random.plain();
  const double a1 = cos(phi);
  const double a2 = sin(phi);
  double _Complex v1 = a2 + a1 * _Complex_I;
  double _Complex v2 =-a2 + a1 * _Complex_I;
  v1 *= sin_alpha*sin_theta;
  v2 *= sin_alpha*sin_theta;

  a[0] = v0 * u0 + v1 * u2;
  a[1] = v0 * u1 + v1 * u3;
  a[2] = v2 * u0 + v3 * u2;
  a[3] = v2 * u1 + v3 * u3;

  return true;
}

bool HeatBathSU3_Polyakov( gauge_field &U , double beta , double lambda ,
			   int n_iter = 1 )
{
  Gauge V(U);
  const double beta_eff = beta / U.nc;

  const int order[6][3] = { { 0,1,2 } , { 0,2,1 } , { 1,0,2 } ,
			    { 1,2,0 } , { 2,0,1 } , { 2,1,0 } };
  const int element[3][2] = { { 1,2 } , { 0,2 } , { 0,1 } };

  for( int iter=0 ; iter<n_iter ; iter++ ){
    for( int parity=0 ; parity<2 ; parity++ ) 
      for( int mu=0 ; mu<U.ndim ; mu++ ){
#pragma omp parallel
	{
	mdp_prng& random = mdp_randomOMP[omp_get_thread_num()];
	Site x( U.lattice() );
	double _Complex tmp_ptr1[9] , tmp_ptr2[9];
	Const_Matrix<3,3,false,false,1> Staple( tmp_ptr1 );
	Const_Matrix<3,3,false,false,1> M( tmp_ptr2 );

	ForAllSitesOfParityOMP( x , parity ) {
	  Staple_Polyakov( V , x,mu , lambda , Staple );
		
	  const int n = (int)( 6.* random.plain() );
	  for( int m=0 ; m<U.nc ; m++ ){
	    const int i = element[ order[n][m] ][0];
	    const int j = element[ order[n][m] ][1];

	    M = V(x,mu) * Staple;
	    double _Complex a[4];
	    a[0] = M(i,i); 
	    a[1] = M(i,j);
	    a[2] = M(j,i);
	    a[3] = M(j,j);

	    if( heatbath_SU2( random , beta_eff , a ) )
	      for( int k=0 ; k<U.nc ; k++ ) {
		const double _Complex tmpUik = a[0] * V(x,mu)(i,k) + a[1] * V(x,mu)(j,k);
		V(x,mu)(j,k) = a[2] * V(x,mu)(i,k) + a[3] * V(x,mu)(j,k);
		V(x,mu)(i,k) = tmpUik;
	      }
	  }
		
	}// loop for allsites 
	}
	Gauge_Comm_Wilson2<1>( V , mu , parity );
	    
      }// loop for mu
  }// loop for iter
  return true;
}

bool OverRelaxationSU3_Polyakov( gauge_field &U, double lambda , int n_iter )
{
  Gauge V(U);

  const int order[6][3] = { { 0,1,2 } , { 0,2,1 } , { 1,0,2 } ,
			    { 1,2,0 } , { 2,0,1 } , { 2,1,0 } };
  const int element[3][2] = { { 1,2 } , { 2,0 } , { 0,1 } };

  if( U.nc != 3 ) mdp << "overrelaxation is written for Nc=3" << endl;

  const int n = (int)( 6.* mdp_random.plain() );

  for(  int iter=0 ; iter<n_iter ; iter++  ) 
    for(  int parity=0 ; parity<2 ; parity++  ) {
      for(  int mu=0 ; mu<U.ndim ; mu++  ){
#pragma omp parallel
	{
	Site x(U.lattice());
	double _Complex tmp_ptr1[9] , tmp_ptr2[9];
	Const_Matrix<3,3,false,false,1> Staple( tmp_ptr1 );
	Const_Matrix<3,3,false,false,1> UV( tmp_ptr2 );

	ForAllSitesOfParityOMP( x , parity ){
	  Staple_Polyakov( V , x,mu , lambda , Staple );

	  for(  int m=0 ; m<3 ; m++  ){

	    UV = V( x,mu ) * Staple;
	    const int i = element[ order[n][m] ][0];
	    const int j = element[ order[n][m] ][1];

	    double _Complex w1[2] , ww[2];
	    w1[0] = UV(i,i) + ~UV(j,j);
 	    w1[1] = UV(i,j) - ~UV(j,i);

	    const double det = 1./ __real__( w1[0]*~w1[0] + w1[1]*~w1[1] );
	    ww[0] = w1[0] * w1[0] - w1[1] * ~w1[1];
	    ww[1] = (  w1[0] + ~w1[0]  ) * w1[1];

	    double _Complex tmp0[3] , tmp1[3];
	    double _Complex* const Ui = V(x,mu).address+i*3;
	    double _Complex* const Uj = V(x,mu).address+j*3;
	    for(  int r=0 ; r<3 ; r++  ){
	      tmp0[r] = Ui[r] * ~ww[0] - ww[1] * Uj[r];
	      tmp1[r] = Ui[r] * ~ww[1] + ww[0] * Uj[r];
	    }
	    for(  int r=0 ; r<3 ; r++  ) Ui[r] = tmp0[r] * det;
	    for(  int r=0 ; r<3 ; r++  ) Uj[r] = tmp1[r] * det;

	  }// loop for m
	}// loop for x
	}
	Gauge_Comm_Wilson2<1>( V , mu , parity );
      }// loop for mu
    }
  return true;
}

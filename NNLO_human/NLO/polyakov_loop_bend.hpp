//
// evaluate the polyakov loop expectation value
// by M.K.

// This algorithm works only when
// (1) temporal direction is not separated by MPI
// (2) 
template< int N >
inline
double _Complex Polyakov_bend_local( const Gauge& V , Site x )
{
  const int Nt = V.lattice().size(0);
  double _Complex tmp_ptr1[9] , tmp_ptr2[9];
  Const_Matrix<3,3,false,false,1> tmp1( tmp_ptr1 );
  Const_Matrix<3,3,false,false,1> tmp2( tmp_ptr2 );

  tmp2 = 0.;
  for( int nu=1 ; nu<4 ; nu++ ){
    Site y = x + 0;
    tmp1 = V(x,nu) * V(x+nu,0);
    for( int i=0 ; i<N-1 ; i++ ){
      tmp1 = tmp1 * V(y+nu,0);
      y = y+0;
    }
    Site xnd = x-nu;
    tmp2 += tmp1 * Hermite(V(y,nu));
    y = x + 0 - nu;
    tmp1 = Hermite(V(xnd,nu)) * V(xnd,0);
    for( int i=0 ; i<N-1 ; i++ ){
      tmp1 = tmp1 * V(y,0);
      y = y+0;
    }
    tmp2 += tmp1 * V(y,nu);
  }
  for( int t=1 ; t<N ; t++ )
    x = x+0;
  for( int t=N ; t<Nt ; t++ ){
    x = x+0;
    tmp2 = tmp2 * V(x,0);
  }
  return ( tmp2(0,0) + tmp2(1,1) + tmp2(2,2) ) / 3. / 6.;
}
template< int N >
double _Complex Polyakov_bend( const Gauge& V )
{
  const int vol = V.lattice().nx[0] * V.lattice().nx[1] * V.lattice().nx[2] * V.lattice().nx[3];

  double _Complex ans = 0.;
#pragma omp parallel reduction( + : ans )
  {
    Site x( V.lattice() );

    ForAllSitesOMP( x ){
      ans += Polyakov_bend_local<N>( V , x );
    }
    ans /= (double)vol;
  }
  double _Complex tmp;
  MPI_Allreduce( (void*)&ans , (void*)&tmp , 1 ,
		 MPI_DOUBLE_COMPLEX , MPI_SUM , MPI_COMM_WORLD );
  return tmp;
}

inline double _Complex Polyakov_bendN( const Gauge& V , int n )
{
  switch(n){
  case 1:
    return Polyakov_bend<1>(V);
  case 2:
    return Polyakov_bend<2>(V);
  case 3:
    return Polyakov_bend<3>(V);
  case 4:
    return Polyakov_bend<4>(V);
  case 5:
    return Polyakov_bend<5>(V);
  case 6:
    return Polyakov_bend<6>(V);
  default:
    throw std::out_of_range("Polyakov_bendN: n must be 1..6");
  }
}

inline
double _Complex Polyakov_bend1_local( const Gauge& V , Site x )
{
  const int Nt = V.lattice().size(0);
  double _Complex tmp_ptr1[9] , tmp_ptr2[9];
  Const_Matrix<3,3,false,false,1> tmp1( tmp_ptr1 );
  Const_Matrix<3,3,false,false,1> tmp2( tmp_ptr2 );

  tmp2 = 0.;
  for( int nu=1 ; nu<4 ; nu++ ){
    Site xnd = x-nu;
    tmp1 = V(x,nu) * V(x+nu,0);
    tmp2 += tmp1 * Hermite(V(x+0,nu));
    tmp1 = Hermite(V(xnd,nu)) * V(xnd,0);
    tmp2 += tmp1 * V(xnd+0,nu);
  }
  for( int t=1 ; t<Nt ; t++ ){
    x = x+0;
    tmp2 = tmp2 * V(x,0);
  }
  return ( tmp2(0,0) + tmp2(1,1) + tmp2(2,2) ) / 3. / 6.;
}
double _Complex Polyakov_bend1( const Gauge& V )
{
  const int vol = V.lattice().nx[0] * V.lattice().nx[1] * V.lattice().nx[2] * V.lattice().nx[3];

  double _Complex ans = 0.;
#pragma omp parallel reduction( + : ans )
  {
    Site x( V.lattice() );

    ForAllSitesOMP( x ){
      ans += Polyakov_bend1_local( V , x );
    }
    ans /= (double)vol;
  }
  double _Complex tmp;
  MPI_Allreduce( (void*)&ans , (void*)&tmp , 1 ,
		 MPI_DOUBLE_COMPLEX , MPI_SUM , MPI_COMM_WORLD );
  return tmp;
}

inline
double _Complex Polyakov_bend2_local( const Gauge& V , Site x )
{
  const int Nt = V.lattice().size(0);
  double _Complex tmp_ptr1[9] , tmp_ptr2[9];
  Const_Matrix<3,3,false,false,1> tmp1( tmp_ptr1 );
  Const_Matrix<3,3,false,false,1> tmp2( tmp_ptr2 );

  tmp2 = 0.;
  for( int nu=1 ; nu<4 ; nu++ ){
    Site y = x + 0;
    Site xnd = x-nu;
    tmp1 = V(x,nu) * V(x+nu,0);
    tmp1 = tmp1 * V(y+nu,0);
    tmp2 += tmp1 * Hermite(V(y+0,nu));
    y = y - nu;
    tmp1 = Hermite(V(xnd,nu)) * V(xnd,0);
    tmp1 = tmp1 * V(y,0);
    tmp2 += tmp1 * V(y+0,nu);
  }
  x = x+0;
  for( int t=2 ; t<Nt ; t++ ){
    x = x+0;
    tmp2 = tmp2 * V(x,0);
  }
  return ( tmp2(0,0) + tmp2(1,1) + tmp2(2,2) ) / 3. / 6.;
}
double _Complex Polyakov_bend2( const Gauge& V )
{
  const int vol = V.lattice().nx[0] * V.lattice().nx[1] * V.lattice().nx[2] * V.lattice().nx[3];

  double _Complex ans = 0.;
#pragma omp parallel reduction( + : ans )
  {
    Site x( V.lattice() );

    ForAllSitesOMP( x ){
      ans += Polyakov_bend2_local( V , x );
    }
    ans /= (double)vol;
  }
  double _Complex tmp;
  MPI_Allreduce( (void*)&ans , (void*)&tmp , 1 ,
		 MPI_DOUBLE_COMPLEX , MPI_SUM , MPI_COMM_WORLD );
  return tmp;
}

inline
double _Complex Polyakov_bend_half_local( const Gauge& V , Site x )
{
  const int Nt = V.lattice().size(0);
  const int Nt2 = Nt/2;
  double _Complex tmp_ptr1[9] , tmp_ptr2[9];
  Const_Matrix<3,3,false,false,1> tmp1( tmp_ptr1 );
  Const_Matrix<3,3,false,false,1> tmp2( tmp_ptr2 );

  tmp2 = 0.;
  for( int nu=1 ; nu<4 ; nu++ ){
    Site y = x;
    tmp1 = V(x,nu) * V(x+nu,0);
    for( int t=1 ; t<Nt2 ; t++ ){
      y = y+0;
      tmp1 = tmp1 * V(y+nu,0);
    }
    tmp2 += tmp1 * Hermite(V(y+0,nu));
  }
  for( int t=1 ; t<Nt2 ; t++ )
    x = x+0;
  for( int t=0 ; t<Nt2 ; t++ ){
    x = x+0;
    tmp2 = tmp2 * V(x,0);
  }
  return ( tmp2(0,0) + tmp2(1,1) + tmp2(2,2) ) / 3. / 3.;
}
double _Complex Polyakov_bend_half( const Gauge& V )
{
  const int vol = V.lattice().nx[0] * V.lattice().nx[1] * V.lattice().nx[2] * V.lattice().nx[3];

  double _Complex ans = 0.;
#pragma omp parallel reduction( + : ans )
  {
    Site x( V.lattice() );

    ForAllSitesOMP( x ){
      ans += Polyakov_bend_half_local( V , x );
    }
    ans /= (double)vol;
  }
  double _Complex tmp;
  MPI_Allreduce( (void*)&ans , (void*)&tmp , 1 ,
		 MPI_DOUBLE_COMPLEX , MPI_SUM , MPI_COMM_WORLD );
  return tmp;
}


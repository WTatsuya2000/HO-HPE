//
// evaluate the polyakov loop expectation value
// by M.K.

// This algorithm works only when
// (1) temporal direction is not separated by MPI
// (2) 
inline
double _Complex Polyakov_local( const Gauge& V , Site x , 
				Const_Matrix<3,3,false,false,1> tmp )
{
  const int Nt = V.lattice().size(0);

  tmp = V(x,0);
  for( int t=1 ; t<Nt ; t++ ){
    x = x+0;
    tmp = tmp * V(x,0);
  }
  return ( tmp(0,0) + tmp(1,1) + tmp(2,2) ) / 3.;
}

double _Complex Polyakov_loop( const Gauge& V )
{
  const int vol = V.lattice().nx[1] * V.lattice().nx[2] * V.lattice().nx[3];

  double _Complex ans = 0.;
#pragma omp parallel reduction( + : ans )
  {
    Site x( V.lattice() );
    double _Complex p_tmp[9];
    Const_Matrix<3,3,false,false,1> tmp(p_tmp);

    ForAllSitesOMP( x ){
      if( x(0)==0 )
	ans += Polyakov_local( V , x , tmp );
    }
    ans /= (double)vol;
  }
  double _Complex tmp;
  MPI_Allreduce( (void*)&ans , (void*)&tmp , 1 ,
		 MPI_DOUBLE_COMPLEX , MPI_SUM , MPI_COMM_WORLD );
  return tmp;
}

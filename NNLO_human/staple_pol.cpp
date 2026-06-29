/*
inline void staple_H( const Gauge& U , const Site& x , int mu , Const_Matrix<3,3,false,false,1>& ans )
{
  double _Complex tmp_ptr[9];
  Const_Matrix<3,3,false,false,1> tmp( tmp_ptr );
  
  ans = 0.;

  for( int nu=0; nu<4; nu++)if(nu!=mu){
      Site xnd = x-nu;
      tmp = U(x+mu,nu) * Hermite(U(x+nu,mu));
      ans += tmp * Hermite(U(x,nu));
      tmp = Hermite(U(xnd+mu,nu)) * Hermite(U(xnd,mu));
      ans += tmp * U(xnd,nu);
    }
}
*/
inline void Staple_Polyakov( const Gauge& U , const Site& x , int mu ,
			     double lambda , 
			     Const_Matrix<3,3,false,false,1>& ans )
{
  double _Complex tmp_ptr[9];
  Const_Matrix<3,3,false,false,1> tmp( tmp_ptr );
  
  ans = 0.;

  for( int nu=0; nu<4; nu++)if(nu!=mu){
      Site xnd = x-nu;
      tmp = U(x+mu,nu) * Hermite(U(x+nu,mu));
      ans += tmp * Hermite(U(x,nu));
      tmp = Hermite(U(xnd+mu,nu)) * Hermite(U(xnd,mu));
      ans += tmp * U(xnd,nu);
    }
  if( mu==0 ){
    const int Nt = U.lattice().size(0);
    Site y = x+0;
    tmp = U( y,0 );
    for( int t=1 ; t<Nt-1 ; t++ ){
      y = y+0;
      tmp = tmp * U( y,0);
    }
    ans += lambda * tmp; // \beta( S_g - \lambda \Omega )
  }
}


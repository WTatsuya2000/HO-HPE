class SingleStaple: public mdp_complex_field {
  const int Nelem = 36;
public:
  double _Complex* m;

  SingleStaple( const Gauge& V ){
    allocate_field( V.lattice() , 9*Nelem );
    m = (double _Complex*)physical_address();
    CalcStaple( V );
  }
    
    // mu direction , on the (mu,nu) plane, direction sign(=+1,-1; 0:sum)
    int index( int mu , int nu , int sign ) const {
        static constexpr int index[4][4] = { {-1,0,1,2} , {3,-1,4,5} , {6,7,-1,8} , {9,10,11,-1} };
        return 3*index[mu][nu] + (sign+1);
    }

    Matrix<3,3> operator()( const Site& x , int mu , int nu , int dir ){
        return Matrix<3,3>( m + x.idx*9*Nelem + 9*index(mu,nu,dir) );
    }
    Matrix<3,3> operator()( const Site& x , int mu , int nu , int dir ) const {
        return Matrix<3,3>( m + x.idx*9*Nelem + 9*index(mu,nu,dir) );
    }
    
    void CalcStaple( const Gauge& V ){
        #pragma omp parallel 
        {
        double _Complex ptr[9];
        Matrix<3,3> tmp(ptr);
        Site x( V.lattice() );
        ForAllSitesOMP(x){
            for( int mu : {0,1,2,3} ){
                for( int nu : {0,1,2,3} ){
                    if( mu == nu ) continue;
                    Site xn = x - nu;
                    tmp = V(x,nu) * V(x+nu,mu);
                    (*this)(x,mu,nu,1) = tmp * Hermite( V(x+mu,nu) );
                    tmp = Hermite( V(xn,nu) ) * V(xn,mu);
                    (*this)(x,mu,nu,-1) = tmp * V(xn+mu,nu);

		    (*this)(x,mu,nu,0) = (*this)(x,mu,nu,1) + (*this)(x,mu,nu,-1);
                }
            }
        }
        }
        this->update();
    }
};

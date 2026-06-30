#include <array>

struct NNNLO_WPath
{
  double coeff;
  int traj[10];
};

struct NNNLO_WNode
{
  int dir;
  double coeff;
  std::vector<int> child;
};

struct NNNLO_WPair
{
  int prefix;
  int suffix;
  double coeff;
};

std::array<int,10> NNNLO_cyclic_min( const std::vector<int>& traj )
{
  std::array<int,10> ans;
  for( int n=0 ; n<10 ; n++ )
    ans[n] = traj[n];
  for( int shift=1 ; shift<10 ; shift++ ){
    std::array<int,10> tmp;
    for( int n=0 ; n<10 ; n++ )
      tmp[n] = traj[(shift+n)%10];
    if( tmp < ans )
      ans = tmp;
  }
  return ans;
}

std::vector<NNNLO_WPath>& NNNLO_W_paths()
{
  static std::vector<NNNLO_WPath> paths;
  if( paths.empty() ){
    for( int i=0 ; i<W10traj.size() ; i++ ){
      if( W10D[i]==0 ) continue;
      const double coeff = W10D[i] / (double)cyclic_symmetry(W10traj[i]);
      const auto trajlist = MakeTrajList(W10traj[i]);
      for( const auto& tr : trajlist ){
	const auto traj = NNNLO_cyclic_min( tr );
	NNNLO_WPath path;
	path.coeff = coeff;
	for( int n=0 ; n<10 ; n++ )
	  path.traj[n] = traj[n];
	paths.push_back( path );
      }
    }
  }
  return paths;
}

std::vector<NNNLO_WNode>& NNNLO_W_tree()
{
  static std::vector<NNNLO_WNode> tree;
  if( tree.empty() ){
    tree.push_back( {0,0.,std::vector<int>()} );
    const auto& paths = NNNLO_W_paths();
    for( const auto& path : paths ){
      int node = 0;
      for( int depth=0 ; depth<10 ; depth++ ){
	const int dir = path.traj[depth];
	int next = -1;
	for( int c : tree[node].child ){
	  if( tree[c].dir == dir ){
	    next = c;
	    break;
	  }
	}
	if( next<0 ){
	  next = tree.size();
	  tree.push_back( {dir,0.,std::vector<int>()} );
	  tree[node].child.push_back( next );
	}
	node = next;
      }
      tree[node].coeff += path.coeff;
    }
  }
  return tree;
}

std::vector<NNNLO_WNode>& NNNLO_W_prefix_tree()
{
  static std::vector<NNNLO_WNode> tree;
  if( tree.empty() ){
    tree.push_back( {0,0.,std::vector<int>()} );
    const auto& paths = NNNLO_W_paths();
    for( const auto& path : paths ){
      int node = 0;
      for( int depth=0 ; depth<5 ; depth++ ){
	const int dir = path.traj[depth];
	int next = -1;
	for( int c : tree[node].child ){
	  if( tree[c].dir == dir ){
	    next = c;
	    break;
	  }
	}
	if( next<0 ){
	  next = tree.size();
	  tree.push_back( {dir,0.,std::vector<int>()} );
	  tree[node].child.push_back( next );
	}
	node = next;
      }
    }
  }
  return tree;
}

std::vector<NNNLO_WNode>& NNNLO_W_suffix_tree()
{
  static std::vector<NNNLO_WNode> tree;
  if( tree.empty() ){
    tree.push_back( {0,0.,std::vector<int>()} );
    const auto& paths = NNNLO_W_paths();
    for( const auto& path : paths ){
      int node = 0;
      for( int depth=9 ; depth>=5 ; depth-- ){
	const int dir = -path.traj[depth];
	int next = -1;
	for( int c : tree[node].child ){
	  if( tree[c].dir == dir ){
	    next = c;
	    break;
	  }
	}
	if( next<0 ){
	  next = tree.size();
	  tree.push_back( {dir,0.,std::vector<int>()} );
	  tree[node].child.push_back( next );
	}
	node = next;
      }
    }
  }
  return tree;
}

int NNNLO_W_half_leaf( const std::vector<NNNLO_WNode>& tree ,
		       const int *traj , bool suffix )
{
  int node = 0;
  if( suffix ){
    for( int depth=9 ; depth>=5 ; depth-- ){
      const int dir = -traj[depth];
      for( int c : tree[node].child ){
	if( tree[c].dir == dir ){
	  node = c;
	  break;
	}
      }
    }
  }
  else{
    for( int depth=0 ; depth<5 ; depth++ ){
      const int dir = traj[depth];
      for( int c : tree[node].child ){
	if( tree[c].dir == dir ){
	  node = c;
	  break;
	}
      }
    }
  }
  return node;
}

std::vector<NNNLO_WPair>& NNNLO_W_pairs()
{
  static std::vector<NNNLO_WPair> pairs;
  if( pairs.empty() ){
    const auto& ptree = NNNLO_W_prefix_tree();
    const auto& stree = NNNLO_W_suffix_tree();
    const auto& paths = NNNLO_W_paths();
    for( const auto& path : paths ){
      NNNLO_WPair pair;
      pair.prefix = NNNLO_W_half_leaf( ptree , path.traj , false );
      pair.suffix = NNNLO_W_half_leaf( stree , path.traj , true );
      pair.coeff = path.coeff;
      pairs.push_back( pair );
    }
  }
  return pairs;
}

long long NNNLO_W_v1_matrix_multiplications_per_site()
{
  return (long long)NNNLO_W_tree().size() - 1LL;
}

long long NNNLO_W_v2_matrix_multiplications_per_site()
{
  return (long long)NNNLO_W_prefix_tree().size() - 1LL
    + (long long)NNNLO_W_suffix_tree().size() - 1LL
    + (long long)NNNLO_W_pairs().size();
}

long long NNNLO_W_matrix_multiplications_per_site()
{
  return (long long)NNNLO_W_prefix_tree().size() - 1LL
    + (long long)NNNLO_W_suffix_tree().size() - 1LL;
}

long long NNLO_W10_safe_matrix_multiplications_per_site()
{
  long long ans = 0;
  for( int i=0 ; i<W10traj.size() ; i++ ){
    if( W10D[i]==0 ) continue;
    ans += 28LL * (long long)MakeTrajList(W10traj[i]).size();
  }
  return ans;
}

double NNNLO_W_tree_local( const Gauge& V , const Site& x ,
			   const std::vector<NNNLO_WNode>& tree ,
			   int node , int depth , double _Complex *Uptr )
{
  Matrix<3,3> tmp( Uptr + 9*depth );
  double ans = tree[node].coeff * __real__( Trace(tmp) );
  for( int c : tree[node].child ){
    const int dir = tree[c].dir;
    const int mu = std::abs(dir) % 4;
    Site xh = x;
    Matrix<3,3> tmp1( Uptr + 9*(depth+1) );
    if( dir>0 ){
      tmp1 = tmp * V(x,mu);
      xh = x + mu;
    }
    else{
      xh = x - mu;
      tmp1 = tmp * Hermite(V(xh,mu));
    }
    ans += NNNLO_W_tree_local( V , xh , tree , c , depth+1 , Uptr );
  }
  return ans;
}

double NNNLO_W_v1( const Gauge& V )
{
  double ans = 0.;
  const auto& tree = NNNLO_W_tree();
#pragma omp parallel reduction( + : ans )
  {
    double _Complex ptr[99];
    Site x(V.lattice());
    ForAllSitesOMP(x){
      Matrix<3,3> tmp( ptr );
      tmp = 1.;
      ans += NNNLO_W_tree_local( V , x , tree , 0 , 0 , ptr );
    }
  }
  mpi.add(ans);
  return -2. * ans / V.lattice().size();
}

void NNNLO_W_eval_half_tree( const Gauge& V , const Site& x ,
			     const std::vector<NNNLO_WNode>& tree ,
			     int node , double _Complex *Uptr )
{
  Matrix<3,3> tmp( Uptr + 9*node );
  for( int c : tree[node].child ){
    const int dir = tree[c].dir;
    const int mu = std::abs(dir) % 4;
    Matrix<3,3> tmp1( Uptr + 9*c );
    if( dir>0 )
      tmp1 = tmp * V(x,mu);
    else
      tmp1 = tmp * Hermite(V(x-mu,mu));
    Site xh = (dir>0) ? x+mu : x-mu;
    NNNLO_W_eval_half_tree( V , xh , tree , c , Uptr );
  }
}

double NNNLO_W_trace_AHB( const Matrix<3,3>& A , const Matrix<3,3>& B )
{
  double ans = 0.;
  for( int i=0 ; i<9 ; i++ )
    ans += __real__( A.address[i] * conj(B.address[i]) );
  return ans;
}

double NNNLO_W_v2( const Gauge& V )
{
  double ans = 0.;
  const auto& ptree = NNNLO_W_prefix_tree();
  const auto& stree = NNNLO_W_suffix_tree();
  const auto& pairs = NNNLO_W_pairs();
#pragma omp parallel reduction( + : ans )
  {
    std::vector<double _Complex> pptr( 9*ptree.size() );
    std::vector<double _Complex> sptr( 9*stree.size() );
    double _Complex tptr[9];
    Site x(V.lattice());
    ForAllSitesOMP(x){
      Matrix<3,3> p0( pptr.data() ) , s0( sptr.data() );
      p0 = 1.;
      s0 = 1.;
      NNNLO_W_eval_half_tree( V , x , ptree , 0 , pptr.data() );
      NNNLO_W_eval_half_tree( V , x , stree , 0 , sptr.data() );
      for( const auto& pair : pairs ){
	Matrix<3,3> tmp( tptr );
	tmp = Matrix<3,3>(pptr.data() + 9*pair.prefix)
	  * Hermite( Matrix<3,3>(sptr.data() + 9*pair.suffix) );
	ans += pair.coeff * __real__( Trace(tmp) );
      }
    }
  }
  mpi.add(ans);
  return -2. * ans / V.lattice().size();
}

double NNNLO_W( const Gauge& V )
{
  double ans = 0.;
  const auto& ptree = NNNLO_W_prefix_tree();
  const auto& stree = NNNLO_W_suffix_tree();
  const auto& pairs = NNNLO_W_pairs();
#pragma omp parallel reduction( + : ans )
  {
    std::vector<double _Complex> pptr( 9*ptree.size() );
    std::vector<double _Complex> sptr( 9*stree.size() );
    Site x(V.lattice());
    ForAllSitesOMP(x){
      Matrix<3,3> p0( pptr.data() ) , s0( sptr.data() );
      p0 = 1.;
      s0 = 1.;
      NNNLO_W_eval_half_tree( V , x , ptree , 0 , pptr.data() );
      NNNLO_W_eval_half_tree( V , x , stree , 0 , sptr.data() );
      for( const auto& pair : pairs )
	ans += pair.coeff * NNNLO_W_trace_AHB(
	  Matrix<3,3>(pptr.data() + 9*pair.prefix) ,
	  Matrix<3,3>(sptr.data() + 9*pair.suffix) );
    }
  }
  mpi.add(ans);
  return -2. * ans / V.lattice().size();
}

double NNNLO_W10_7_try_local( const Gauge& V , const SingleStaple& S ,
			      const LocalElements& LE , const Site& x ,
			      double _Complex *Uptr )
{
  Matrix<3,3> tmp( Uptr ) , st1( Uptr+9 ) , st2( Uptr+18 );
  double ans = 0.;
  for( const auto& m : HPE_comb::C2 ){
    const int mu=m[0] , nu=m[1];
    for( int sm : {-1,1} ){
      for( int sn : {-1,1} ){
	if( sm==1 ){
	  st1 = S(x,mu,nu,-sn);
	  st2 = S(x,mu,nu, sn);
	}
	else{
	  st1 = Hermite( S(x-mu,mu,nu,-sn) );
	  st2 = Hermite( S(x-mu,mu,nu, sn) );
	}
	tmp = LE.P(mu,nu,-sm,sn) * st1;
	tmp = tmp * Hermite( st2 );
	ans += __real__( Trace(tmp) );
      }
    }
  }
  return ans;
}

double NNNLO_W10_7_try( const Gauge& V )
{
  const SingleStaple& S(V);
  double ans = 0.;
#pragma omp parallel reduction( + : ans )
  {
    double _Complex ptr[27];
    Site x(V.lattice());
    ForAllSitesOMP(x){
      LocalElements LE( V , S , x );
      ans += NNNLO_W10_7_try_local( V , S , LE , x , ptr );
    }
  }
  mpi.add(ans);
  return -2. * ans / V.lattice().size() * W10D[6];
}

double NNNLO_W10_safe_i( const Gauge& V , int i )
{
  double tmp = __real__(TrajSum_safe( V , MakeTrajList(W10traj[i])) );
  tmp *= W10D[i] / cyclic_symmetry(W10traj[i]);
  return tmp;
}

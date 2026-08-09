#include <array>
#include <vector>
#include <fstream>
#include <algorithm>
#include <cstdio>

// ===========================================================================
//  hpe_pool.cpp : length-generic HPE measurement engine (lengths 6 ... 12),
//                 i.e. NLO ... N4LO, based on a shared trie ("pool").
//
//  Wilson type   : symmetry orbit = MakeTrajList  (24 perm x 16 sign),  -2/V
//  Polyakov type : symmetry orbit = MakeTrajListP (6 perm x 8 sign, time
//                  axis fixed), classified by the temporal winding number.
//                  Coefficients carry the antiperiodic-BC sign (-1)^(w-1).
//
//  Each trajectory is cut at a chosen starting point into a forward and a
//  backward subtrajectory of equal length L/2.  The backward half is reversed
//  and negated so that both halves are ordinary forward link products, and the
//  loop trace is recovered as Tr( P . S^dagger ).
//
//  The two halves are therefore sequences drawn from the same alphabet with
//  the same length, so they are stored in ONE trie rather than two.  Whenever
//  a direction sequence occurs both as a forward and as a backward half it is
//  then evaluated only once, and matrix products shared by several
//  subtrajectories are shared across the whole pool.  Sectors may in addition
//  share a single pool with each other (HPEP_pool_mode below), which lets a
//  short sector ride for free on the internal nodes of a longer one.
//
//  Cold (U=1) reference values, reproduced by loop_classification/:
//    W6  = 8448        W8  = 245952      W10 = 7372800   W12 = 225232896
//    L1(4,6)  = 1728
//    L1(4,8)  = 45792        L2(4,8)  = -384
//    L1(4,10) = 645120       L2(4,10) = -64512
//    L1(4,12) = -26224128    L2(4,12) = -6842880    L3(4,12) = 4096
//
//  External dependencies (supplied by the production driver, see README):
//    lattice library ......... Gauge, Site, Matrix<3,3>, Hermite,
//                              ForAllSitesOMP, mpi.add, mpicout
//    NNLO_human/nnlo_safe.cpp  cyclic_symmetry, MakeTrajList, MakeTrajListP
//
//  Defining HPE_COST_ONLY before including this file drops everything that
//  touches the lattice library and keeps only table reading, trie building and
//  the cost counters.  That build needs nothing but the symmetry helpers above
//  and is how the matrix-multiplication counts quoted in the paper are
//  obtained (see src/hpe_cost.cpp).
// ===========================================================================

#ifndef HPE_POOL_CPP
#define HPE_POOL_CPP

// Directory holding the coefficient tables, relative to the working directory.
#ifndef HPE_TABLE_DIR
#define HPE_TABLE_DIR "trajectory_list/tables"
#endif

struct HPEPNode { int dir; std::vector<int> child; };
struct HPEPPair { int prefix; int suffix; double coeff; };   // leaves of the same trie

// ---- choice of the starting point (anchor) --------------------------------
// For a cut after k links, net(d) = (#(+d) in the forward half) - (#(-d) in it)
// is the displacement between the two endpoints of the forward subtrajectory.
// Since the trajectory is closed, the backward half has the opposite net, so a
// cut that makes |net| small makes the two halves geometrically similar and
// maximises the sharing in the pool.
//
//   0 = reach_min  minimise the spatial bounding-box reach.  This minimises the
//                  MPI halo depth rather than the cost, and is the right choice
//                  for winding Polyakov loops.
//   1 = split_min  minimise sum_d net(d)^2.
//   2 = auto       split_min for Wilson, reach_min for Polyakov.
//   3 = paper      minimise sum_d |net(d)|, i.e. the rule stated in the paper.
//
// Measured cost in matrix multiplications per site (Nt = 4, Lmax = 12), one
// pool per sector, for the two scoring rules:
//                reach_min   split_min
//   W10             8460        4815
//   W12            70307       34217
//   L1(4,12)       30293       38834
// A winding Polyakov loop has a structurally nonzero net in the time direction,
// so trying to balance it is counterproductive; that is what mode 2 avoids.
int HPEP_anchor_mode = 1;

// How sectors are grouped into pools:
//   0 = one pool per sector
//   1 = all Wilson sectors share pool 0, each Polyakov sector gets its own
//   2 = every sector shares pool 0
// Total cost over all 12 sectors at Nt = 4, Lmax = 12 (MM/site):
//                  reach_min   split_min   auto
//   per sector       115396      85608     75157
//   Wilson shared    105645      80006     69555
//   all shared        72703      51596     53166
int HPEP_pool_mode = 2;

// Merge pairs that end up on the same (prefix, suffix) leaves by adding their
// coefficients, and drop pairs whose coefficients cancel.
bool HPEP_merge_pairs = true;

// Group pairs sharing (prefix, coefficient) so that the suffixes are summed
// before a single trace is taken.  This trades traces for gathers: it removes
// flops but adds a variable-length inner loop.  Measured on the production
// machine (Intel icx / EPYC, 16 threads) it was a net loss (1.09 s -> 1.56 s),
// while GCC gained about 10%.  Off by default.
bool HPEP_use_star = false;

struct HPEPPool { std::vector<HPEPNode> tree; };
std::vector<HPEPPool> gp_pools;

struct HPEPSector
{
  int half = 0;
  int maxreach = 0;                    // spatial reach = MPI halo this sector needs
  double factor = 0.;                  // -2 (Wilson) / +2 (Polyakov)
  int pool = 0;                        // which shared pool it lives in
  std::vector<HPEPPair> pairs;         // leaf indices into gp_pools[pool].tree

  // Star-shaped CSR view of the pair list, used when HPEP_use_star is set.
  struct StarGroup { double coeff; int p9; int off; int len; };
  std::vector<StarGroup> star;
  std::vector<int>       sidx;         // concatenated suffix leaves (9 * node)
  long long npair_raw = 0;             // pair count before merging
};

// ---- geometric reach ------------------------------------------------------
// Walk the vertices of the trajectory from its starting point and take the
// per-axis bounding box; the reach is the symmetric halo width this anchoring
// needs.  SPATIAL axes only: the time direction is never MPI-split, so it
// requires no halo, and including it would inflate the reach of winding
// Polyakov loops (up to ~m*Nt) and push the trie through the spatial halo.
static int hpep_reach( const std::vector<int>& c )
{
  int m[4]={0,0,0,0}, M[4]={0,0,0,0}, p[4]={0,0,0,0};
  const int L = c.size();
  for( int k=0 ; k<L-1 ; k++ ){          // the last step closes the loop
    const int a = std::abs(c[k]) - 1;
    p[a] += (c[k]>0) ? 1 : -1;
    if( p[a] < m[a] ) m[a] = p[a];
    if( p[a] > M[a] ) M[a] = p[a];
  }
  int r = 0;
  for( int a=0 ; a<3 ; a++ ){
    const int e = std::max( M[a], -m[a] );
    if( e > r ) r = e;
  }
  return r;
}

static void hpep_net( const std::vector<int>& c , int half , int net[4] )
{
  net[0]=net[1]=net[2]=net[3]=0;
  for( int k=0 ; k<half ; k++ ){
    const int a = std::abs(c[k]) - 1;
    net[a] += (c[k]>0) ? 1 : -1;
  }
}

// Minimise sum_d net(d)^2; ties broken by maximising net(1), net(2), net(3),
// net(4) in that order, then by the smallest k.  Cutting at k and at k+L/2
// gives the same score with opposite net, so the tie-break fixes one of them.
static std::vector<int> hpep_split_min( const std::vector<int>& t )
{
  const int L = t.size() , half = L/2;
  std::vector<int> c( L );
  int bestNet[4] , bestScore = 0 , bestk = 0;
  hpep_net( t , half , bestNet );
  for( int d=0 ; d<4 ; d++ ) bestScore += bestNet[d]*bestNet[d];

  for( int s=1 ; s<L ; s++ ){
    for( int n=0 ; n<L ; n++ ) c[n] = t[(s+n)%L];
    int net[4]; hpep_net( c , half , net );
    int sc = 0; for( int d=0 ; d<4 ; d++ ) sc += net[d]*net[d];

    bool better = false;
    if( sc < bestScore ) better = true;
    else if( sc == bestScore ){
      for( int d=0 ; d<4 ; d++ )
	if( net[d] != bestNet[d] ){ better = ( net[d] > bestNet[d] ); break; }
    }
    if( better ){
      bestScore = sc; bestk = s;
      for( int d=0 ; d<4 ; d++ ) bestNet[d] = net[d];
    }
  }
  std::vector<int> best( L );
  for( int n=0 ; n<L ; n++ ) best[n] = t[(bestk+n)%L];
  return best;
}

// The rule as stated in the paper: among the cuts minimising sum_d |net(d)|,
// keep those whose first nonzero net, scanned in the order mu = 4, 1, 2, 3, is
// positive; remaining ties are resolved by the smallest k.
//
// The paper resolves the last tie by minimising the total number of matrix
// multiplications, which is a global combinatorial optimisation over all
// trajectories at once (simulated annealing plus coordinate descent in the
// W(10)-specific implementation).  That search is not reproduced here; the
// deterministic rule below is what the general engine uses.
static std::vector<int> hpep_paper_min( const std::vector<int>& t )
{
  static const int ORD[4] = { 3 , 0 , 1 , 2 };      // mu = 4, 1, 2, 3
  const int L = t.size() , half = L/2;
  std::vector<int> c( L );
  int bestScore = -1 , bestk = -1;

  for( int s=0 ; s<L ; s++ ){
    for( int n=0 ; n<L ; n++ ) c[n] = t[(s+n)%L];
    int net[4]; hpep_net( c , half , net );
    int sc = 0; for( int d=0 ; d<4 ; d++ ) sc += std::abs(net[d]);

    bool oriented = true;                           // all-zero net counts as oriented
    for( int i=0 ; i<4 ; i++ )
      if( net[ORD[i]] != 0 ){ oriented = ( net[ORD[i]] > 0 ); break; }
    if( !oriented ) continue;

    if( bestk < 0 || sc < bestScore ){ bestScore = sc; bestk = s; }
  }
  if( bestk < 0 ) bestk = 0;                        // cannot happen for a closed loop

  std::vector<int> best( L );
  for( int n=0 ; n<L ; n++ ) best[n] = t[(bestk+n)%L];
  return best;
}

static std::vector<int> hpep_reach_min( const std::vector<int>& t )
{
  const int L = t.size();
  std::vector<int> best( t );
  int bestR = hpep_reach( best );
  for( int s=1 ; s<L ; s++ ){
    std::vector<int> c( L );
    for( int n=0 ; n<L ; n++ ) c[n] = t[(s+n)%L];
    const int r = hpep_reach( c );
    if( r < bestR || ( r==bestR && c < best ) ){ bestR = r; best = c; }
  }
  return best;
}

static int hpep_trie_insert( std::vector<HPEPNode>& tree , const std::vector<int>& dirs )
{
  int node = 0;
  for( int d : dirs ){
    int nx = -1;
    for( int c : tree[node].child )
      if( tree[c].dir == d ){ nx = c; break; }
    if( nx<0 ){
      nx = tree.size();
      tree.push_back( {d,std::vector<int>()} );
      tree[node].child.push_back( nx );
    }
    node = nx;
  }
  return node;
}

HPEPSector hpep_build( const std::vector<std::vector<int>>& T ,
		       const std::vector<int>& D , bool poly , int pool )
{
  HPEPSector s;
  s.factor = poly ? 2. : -2.;
  s.pool   = pool;
  if( (int)gp_pools.size() <= pool ) gp_pools.resize( pool+1 );
  std::vector<HPEPNode>& tree = gp_pools[pool].tree;
  if( tree.empty() ) tree.push_back( {0,std::vector<int>()} );   // root = identity
  for( int i=0 ; i<(int)T.size() ; i++ ){
    if( D[i]==0 ) continue;
    const int L = T[i].size();
    s.half = L/2;
    const double coeff = D[i] / (double)cyclic_symmetry(T[i]);
    const auto orbit = poly ? MakeTrajListP(T[i]) : MakeTrajList(T[i]);
    for( const auto& tr : orbit ){
      std::vector<int> cm;
      switch( HPEP_anchor_mode ){
      case 1:  cm = hpep_split_min ( tr ); break;
      case 2:  cm = poly ? hpep_reach_min( tr ) : hpep_split_min( tr ); break;
      case 3:  cm = hpep_paper_min ( tr ); break;
      default: cm = hpep_reach_min ( tr ); break;
      }
      const int rr = hpep_reach( cm );
      if( rr > s.maxreach ) s.maxreach = rr;
      std::vector<int> pre( cm.begin() , cm.begin()+s.half );
      std::vector<int> suf( s.half );
      for( int n=0 ; n<s.half ; n++ ) suf[n] = -cm[L-1-n];
      HPEPPair pr;
      pr.prefix = hpep_trie_insert( tree , pre );   // same trie
      pr.suffix = hpep_trie_insert( tree , suf );   // same trie
      pr.coeff  = coeff;
      s.pairs.push_back( pr );
    }
  }

  // Pairs that landed on the same (prefix, suffix) are summed into one.  In a
  // shared pool distinct trajectories can reach the same leaf pair, so this is
  // a sort-and-collapse rather than a search.
  s.npair_raw = (long long)s.pairs.size();
  if( HPEP_merge_pairs && !s.pairs.empty() ){
    std::sort( s.pairs.begin() , s.pairs.end() ,
	       []( const HPEPPair& a , const HPEPPair& b ){
		 if( a.prefix != b.prefix ) return a.prefix < b.prefix;
		 return a.suffix < b.suffix;
	       } );
    std::vector<HPEPPair> out;
    out.reserve( s.pairs.size() );
    out.push_back( s.pairs[0] );
    for( size_t i=1 ; i<s.pairs.size() ; i++ ){
      HPEPPair& last = out.back();
      if( s.pairs[i].prefix==last.prefix && s.pairs[i].suffix==last.suffix )
	last.coeff += s.pairs[i].coeff;
      else
	out.push_back( s.pairs[i] );
    }
    std::vector<HPEPPair> out2;
    out2.reserve( out.size() );
    for( const auto& p : out ) if( p.coeff != 0. ) out2.push_back( p );
    s.pairs.swap( out2 );
  }

  // ---- star-shaped CSR ----
  {
    std::vector<HPEPPair> sp( s.pairs );
    std::sort( sp.begin() , sp.end() ,
	       []( const HPEPPair& a , const HPEPPair& b ){
		 if( a.prefix != b.prefix ) return a.prefix < b.prefix;
		 if( a.coeff  != b.coeff  ) return a.coeff  < b.coeff;
		 return a.suffix < b.suffix;
	       } );
    s.sidx.reserve( sp.size() );
    for( size_t i=0 ; i<sp.size() ; ){
      size_t j = i;
      while( j<sp.size() && sp[j].prefix==sp[i].prefix && sp[j].coeff==sp[i].coeff ) j++;
      HPEPSector::StarGroup g;
      g.coeff = sp[i].coeff;
      g.p9    = 9*sp[i].prefix;
      g.off   = (int)s.sidx.size();
      g.len   = (int)(j-i);
      for( size_t k=i ; k<j ; k++ ) s.sidx.push_back( 9*sp[k].suffix );
      s.star.push_back( g );
      i = j;
    }
  }
  return s;
}

// Table format: the first token is the number of entries N, then N lines of
//   D  L  d_1 d_2 ... d_L
// where D is the Dirac coefficient of the shape and d_k are signed axis
// indices (1..4 = x,y,z,t, negative = traversed backwards).  Shapes with D = 0
// are omitted from the files.
bool hpep_read_table( const char* path ,
		      std::vector<std::vector<int>>& T , std::vector<int>& D )
{
  std::ifstream f( path );
  if( !f ) return false;
  int N; f >> N;
  if( !f || N < 0 ) return false;
  T.assign( N , std::vector<int>() ); D.assign( N , 0 );
  for( int i=0 ; i<N ; i++ ){
    int d , L; f >> d >> L;
    if( !f ) return false;
    D[i]=d; T[i].resize(L);
    for( int k=0 ; k<L ; k++ ) f >> T[i][k];
  }
  return (bool)f;
}

// ---- registry -------------------------------------------------------------
struct HPEPPolyEntry { int n, m; HPEPSector sec; };
struct HPEPWilEntry  { int n;    HPEPSector sec; };
static std::vector<HPEPPolyEntry> gp_poly;
static std::vector<HPEPWilEntry>  gp_wil;

// Discover every sector that has a table on disk, assign pools, then build.
// Building in increasing pool order keeps the node indices of one pool
// contiguous, which helps locality during evaluation.
void HPEP_init( int Nt , int Lmax = 12 )
{
  gp_pools.clear(); gp_wil.clear(); gp_poly.clear();

  struct Spec { char path[256]; bool poly; int n, m; };
  std::vector<Spec> spec;
  std::vector<std::vector<int>> T; std::vector<int> D;

  for( int n=6 ; n<=Lmax ; n+=2 ){
    Spec sp; sp.poly=false; sp.n=n; sp.m=0;
    snprintf( sp.path , sizeof(sp.path) , "%s/w%d.dat" , HPE_TABLE_DIR , n );
    if( !hpep_read_table( sp.path , T , D ) ) continue;
    spec.push_back( sp );
  }
  for( int n=6 ; n<=Lmax ; n+=2 )
    for( int m=1 ; m<= n/Nt ; m++ ){
      if( (n-m*Nt)%2 != 0 ) continue;
      Spec sp; sp.poly=true; sp.n=n; sp.m=m;
      snprintf( sp.path , sizeof(sp.path) , "%s/p%d_nt%d_w%d.dat" , HPE_TABLE_DIR , n , Nt , m );
      if( !hpep_read_table( sp.path , T , D ) ) continue;
      spec.push_back( sp );
    }

  std::vector<int> poolof( spec.size() , 0 );
  int npool = 1;
  if( HPEP_pool_mode == 0 ){
    for( size_t i=0 ; i<spec.size() ; i++ ) poolof[i] = (int)i;
    npool = (int)spec.size();
  }
  else if( HPEP_pool_mode == 1 ){
    int next = 1;
    for( size_t i=0 ; i<spec.size() ; i++ ) poolof[i] = spec[i].poly ? next++ : 0;
    npool = next;
  }
  if( npool < 1 ) npool = 1;
  gp_pools.resize( npool );

  for( size_t i=0 ; i<spec.size() ; i++ ){
    if( !hpep_read_table( spec[i].path , T , D ) ) continue;
    HPEPSector sec = hpep_build( T , D , spec[i].poly , poolof[i] );
    if( spec[i].poly ) gp_poly.push_back( { spec[i].n , spec[i].m , sec } );
    else               gp_wil .push_back( { spec[i].n , sec } );
  }
}

// ---- cost accounting ------------------------------------------------------
// Every trie node other than the root is produced by exactly one matrix
// multiplication, so the number of matrix multiplications per lattice site
// needed to build a pool is (number of nodes - 1).
long long HPEP_pool_mm( int p ){ return (long long)gp_pools[p].tree.size() - 1LL; }
int       HPEP_npool(){ return (int)gp_pools.size(); }
long long HPEP_total_mm(){
  long long t = 0;
  for( const auto& p : gp_pools ) t += (long long)p.tree.size() - 1LL;
  return t;
}

// One SingleStaple evaluation costs 48 matrix multiplications per site; the
// paper quotes all costs in this unit.
#ifndef HPE_MM_PER_STAPLE
#define HPE_MM_PER_STAPLE 48.0
#endif
double HPEP_staples( long long mm ){ return (double)mm / HPE_MM_PER_STAPLE; }

// Closing a pair costs one Re Tr[X Y] (Wilson) or Tr[X Y] (Polyakov).  For
// SU(N_c) a real trace of a product is worth 1/(2 N_c) of a matrix
// multiplication, a complex one 1/N_c.
#ifndef HPE_NC
#define HPE_NC 3
#endif
double HPEP_trace_mm( long long npair , bool poly )
{
  return (double)npair / ( poly ? (double)HPE_NC : 2.0*(double)HPE_NC );
}

int    HPEP_nwilson(){ return (int)gp_wil.size(); }
int    HPEP_wilson_len( int i ){ return gp_wil[i].n; }
long long HPEP_wilson_mm( int i ){ return HPEP_pool_mm( gp_wil[i].sec.pool ); }
int       HPEP_wilson_pool( int i ){ return gp_wil[i].sec.pool; }
long long HPEP_wilson_npair( int i ){ return (long long)gp_wil[i].sec.pairs.size(); }
long long HPEP_wilson_npair_raw( int i ){ return gp_wil[i].sec.npair_raw; }
int       HPEP_wilson_reach( int i ){ return gp_wil[i].sec.maxreach; }

int    HPEP_npoly(){ return (int)gp_poly.size(); }
int    HPEP_poly_len ( int i ){ return gp_poly[i].n; }
int    HPEP_poly_wind( int i ){ return gp_poly[i].m; }
long long HPEP_poly_mm( int i ){ return HPEP_pool_mm( gp_poly[i].sec.pool ); }
int       HPEP_poly_pool( int i ){ return gp_poly[i].sec.pool; }
long long HPEP_poly_npair( int i ){ return (long long)gp_poly[i].sec.pairs.size(); }
long long HPEP_poly_npair_raw( int i ){ return gp_poly[i].sec.npair_raw; }
int       HPEP_poly_reach( int i ){ return gp_poly[i].sec.maxreach; }

int HPEP_maxreach(){
  int r = 0;
  for( const auto& e : gp_wil )  if( e.sec.maxreach > r ) r = e.sec.maxreach;
  for( const auto& e : gp_poly ) if( e.sec.maxreach > r ) r = e.sec.maxreach;
  return r;
}

#ifndef HPE_COST_ONLY
// ===========================================================================
//  Everything below needs the lattice library.
// ===========================================================================

// Walk the trie once from site x, storing the link product of every node.
static void hpep_eval_tree( const Gauge& V , const Site& x ,
			    const std::vector<HPEPNode>& tree ,
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
    hpep_eval_tree( V , xh , tree , c , Uptr );
  }
}

// Re Tr(A B^dagger) = sum_k [ Re A_k Re B_k + Im A_k Im B_k ]
static inline double hpep_trace_re( const double _Complex* a , const double _Complex* b )
{
  double s = 0.;
  for( int i=0 ; i<9 ; i++ )
    s += __real__(a[i])*__real__(b[i]) + __imag__(a[i])*__imag__(b[i]);
  return s;
}
static inline double _Complex hpep_trace( const double _Complex* a , const double _Complex* b )
{
  double _Complex s = 0.;
  for( int i=0 ; i<9 ; i++ ) s += a[i] * conj(b[i]);
  return s;
}

double hpep_eval_real( const Gauge& V , const HPEPSector& s )
{
  double ar = 0.;
  const auto& tree  = gp_pools[s.pool].tree;
  const auto& pairs = s.pairs;
  const HPEPPair* const PR = pairs.data();
  const int npair = (int)pairs.size();
#pragma omp parallel reduction( + : ar )
  {
    std::vector<double _Complex> uptr( 9*tree.size() );
    double _Complex* const U = uptr.data();
    Site x(V.lattice());
    ForAllSitesOMP(x){
      Matrix<3,3> u0( U );
      u0 = 1.;
      hpep_eval_tree( V , x , tree , 0 , U );
      for( int i=0 ; i<npair ; i++ )
	ar += PR[i].coeff * hpep_trace_re( U + 9*PR[i].prefix , U + 9*PR[i].suffix );
    }
  }
  mpi.add(ar);
  return s.factor * ar / V.lattice().size();
}

double _Complex hpep_eval( const Gauge& V , const HPEPSector& s )
{
  double ar = 0. , ai = 0.;
  const auto& tree  = gp_pools[s.pool].tree;
  const auto& pairs = s.pairs;
  const HPEPPair* const PR = pairs.data();
  const int npair = (int)pairs.size();
#pragma omp parallel reduction( + : ar , ai )
  {
    std::vector<double _Complex> uptr( 9*tree.size() );
    double _Complex* const U = uptr.data();
    Site x(V.lattice());
    ForAllSitesOMP(x){
      Matrix<3,3> u0( U );
      u0 = 1.;
      hpep_eval_tree( V , x , tree , 0 , U );
      for( int i=0 ; i<npair ; i++ ){
	const double _Complex t = hpep_trace( U + 9*PR[i].prefix , U + 9*PR[i].suffix );
	ar += PR[i].coeff * __real__(t);
	ai += PR[i].coeff * __imag__(t);
      }
    }
  }
  mpi.add(ar);
  mpi.add(ai);
  double _Complex z;
  __real__ z = s.factor * ar / V.lattice().size();
  __imag__ z = s.factor * ai / V.lattice().size();
  return z;
}

double          HPEP_wilson_val( const Gauge& V , int i ){ return hpep_eval_real( V , gp_wil[i].sec ); }
double _Complex HPEP_poly_val  ( const Gauge& V , int i ){ return hpep_eval( V , gp_poly[i].sec ); }

// ---- batched measurement --------------------------------------------------
// Sweep the sites once per pool, evaluate the trie once, then run the pair
// lists of every sector living in that pool.  With pool_mode = 2 the lattice is
// swept exactly once for all sectors.
void HPEP_measure_all( const Gauge& V ,
		       std::vector<double>& wilson ,
		       std::vector<double _Complex>& poly )
{
  const int nw = (int)gp_wil.size() , np = (int)gp_poly.size();
  wilson.assign( nw , 0. );
  poly  .assign( np , 0. );
  const double vol = V.lattice().size();

  for( int p=0 ; p<(int)gp_pools.size() ; p++ ){
    std::vector<int> wi , pi;
    for( int i=0 ; i<nw ; i++ ) if( gp_wil [i].sec.pool==p ) wi.push_back(i);
    for( int i=0 ; i<np ; i++ ) if( gp_poly[i].sec.pool==p ) pi.push_back(i);
    if( wi.empty() && pi.empty() ) continue;

    const std::vector<HPEPNode>& tree = gp_pools[p].tree;
    std::vector<double> war( wi.size() , 0. ) , par( pi.size() , 0. ) , pai( pi.size() , 0. );

#pragma omp parallel
    {
      std::vector<double _Complex> uptr( 9*tree.size() );
      double _Complex* const U = uptr.data();
      std::vector<double> lw( wi.size() , 0. ) , lpr( pi.size() , 0. ) , lpi( pi.size() , 0. );
      Site x(V.lattice());
      ForAllSitesOMP(x){
	Matrix<3,3> u0( U );
	u0 = 1.;
	hpep_eval_tree( V , x , tree , 0 , U );

	if( !HPEP_use_star ){
	  for( size_t k=0 ; k<wi.size() ; k++ ){
	    const std::vector<HPEPPair>& pr = gp_wil[wi[k]].sec.pairs;
	    const HPEPPair* const PR = pr.data();
	    const int n = (int)pr.size();
	    double a = 0.;
	    for( int i=0 ; i<n ; i++ )
	      a += PR[i].coeff * hpep_trace_re( U + 9*PR[i].prefix , U + 9*PR[i].suffix );
	    lw[k] += a;
	  }
	  for( size_t k=0 ; k<pi.size() ; k++ ){
	    const std::vector<HPEPPair>& pr = gp_poly[pi[k]].sec.pairs;
	    const HPEPPair* const PR = pr.data();
	    const int n = (int)pr.size();
	    double ar=0. , ai=0.;
	    for( int i=0 ; i<n ; i++ ){
	      const double _Complex t = hpep_trace( U + 9*PR[i].prefix , U + 9*PR[i].suffix );
	      ar += PR[i].coeff * __real__(t);
	      ai += PR[i].coeff * __imag__(t);
	    }
	    lpr[k] += ar; lpi[k] += ai;
	  }
	}
	else{
	  for( size_t k=0 ; k<wi.size() ; k++ ){
	    const auto& sec = gp_wil[wi[k]].sec;
	    const auto* const G = sec.star.data();
	    const int* const SI = sec.sidx.data();
	    const int ng = (int)sec.star.size();
	    double a = 0.;
	    double _Complex B[9];
	    for( int g=0 ; g<ng ; g++ ){
	      const int* si = SI + G[g].off;
	      const double _Complex* s0 = U + si[0];
	      for( int i=0 ; i<9 ; i++ ) B[i] = s0[i];
	      for( int m=1 ; m<G[g].len ; m++ ){
		const double _Complex* sm = U + si[m];
		for( int i=0 ; i<9 ; i++ ) B[i] += sm[i];
	      }
	      a += G[g].coeff * hpep_trace_re( U + G[g].p9 , B );
	    }
	    lw[k] += a;
	  }
	  for( size_t k=0 ; k<pi.size() ; k++ ){
	    const auto& sec = gp_poly[pi[k]].sec;
	    const auto* const G = sec.star.data();
	    const int* const SI = sec.sidx.data();
	    const int ng = (int)sec.star.size();
	    double ar=0. , ai=0.;
	    double _Complex B[9];
	    for( int g=0 ; g<ng ; g++ ){
	      const int* si = SI + G[g].off;
	      const double _Complex* s0 = U + si[0];
	      for( int i=0 ; i<9 ; i++ ) B[i] = s0[i];
	      for( int m=1 ; m<G[g].len ; m++ ){
		const double _Complex* sm = U + si[m];
		for( int i=0 ; i<9 ; i++ ) B[i] += sm[i];
	      }
	      const double _Complex t = hpep_trace( U + G[g].p9 , B );
	      ar += G[g].coeff * __real__(t);
	      ai += G[g].coeff * __imag__(t);
	    }
	    lpr[k] += ar; lpi[k] += ai;
	  }
	}
      }
#pragma omp critical
      {
	for( size_t k=0 ; k<wi.size() ; k++ ) war[k] += lw[k];
	for( size_t k=0 ; k<pi.size() ; k++ ){ par[k] += lpr[k]; pai[k] += lpi[k]; }
      }
    }

    for( size_t k=0 ; k<wi.size() ; k++ ){
      double a = war[k]; mpi.add(a);
      wilson[ wi[k] ] = gp_wil[wi[k]].sec.factor * a / vol;
    }
    for( size_t k=0 ; k<pi.size() ; k++ ){
      double ar = par[k] , ai = pai[k]; mpi.add(ar); mpi.add(ai);
      double _Complex z;
      __real__ z = gp_poly[pi[k]].sec.factor * ar / vol;
      __imag__ z = gp_poly[pi[k]].sec.factor * ai / vol;
      poly[ pi[k] ] = z;
    }
  }
}

#endif  // HPE_COST_ONLY
#endif  // HPE_POOL_CPP

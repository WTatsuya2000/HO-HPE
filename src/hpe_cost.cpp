// ===========================================================================
//  hpe_cost.cpp — count the cost of the trie algorithm without a lattice.
//
//  Building the trie needs only the trajectory tables and the symmetry orbits,
//  none of which touch a gauge field, so the cost of the algorithm can be
//  measured with no lattice library present.  This program does exactly that
//  and reproduces the matrix-multiplication counts quoted in the paper.
//
//  Build (from the repository root):
//    c++ -std=gnu++17 -O2 -fpermissive -I src -o hpe_cost src/hpe_cost.cpp
//
//  Run (from the repository root, so that trajectory_list/tables is found):
//    ./hpe_cost [Nt] [Lmax] [anchor] [pool]
//      Nt      temporal extent selecting the Polyakov tables   (default 4)
//      Lmax    maximum loop length, 12 = N4LO                  (default 12)
//      anchor  0 reach_min, 1 split_min, 2 auto, 3 paper       (default 1)
//      pool    0 per sector, 1 Wilson shared, 2 all shared     (default 0)
//
//  Columns:
//    MM/site   matrix multiplications per lattice site to build the pool
//    staples   the same in units of one SingleStaple evaluation (48 MM)
//    pairs     number of (prefix, suffix) traces closing the trajectories
//    +ReTr     staples including the trace cost, 1/(2 N_c) MM per real trace
//              and 1/N_c per complex one
//    reach     spatial halo depth this sector needs under MPI decomposition
//
//  With pool = 1 or 2 several sectors share one trie, so the per-sector MM
//  column repeats the cost of the shared pool; only the total is meaningful
//  there.  Per-term costs should be read off the default pool = 0.
// ===========================================================================

#define HPE_COST_ONLY
#include "lattice_stub.hpp"
#include "../NNLO_human/nnlo_safe.cpp"   // cyclic_symmetry, MakeTrajList, MakeTrajListP
#include "hpe_pool.cpp"

#include <cstdio>
#include <cstdlib>

static const char* anchor_name( int m )
{
  switch(m){
  case 0: return "reach_min";
  case 1: return "split_min";
  case 2: return "auto";
  case 3: return "paper";
  }
  return "?";
}

static const char* pool_name( int p )
{
  switch(p){
  case 0: return "per-sector";
  case 1: return "wilson-shared";
  case 2: return "all-shared";
  }
  return "?";
}

int main( int argc , char** argv )
{
  setvbuf( stdout , nullptr , _IOLBF , 0 );
  const int Nt     = (argc>1) ? atoi(argv[1]) : 4;
  const int Lmax   = (argc>2) ? atoi(argv[2]) : 12;
  const int anchor = (argc>3) ? atoi(argv[3]) : 1;
  const int pool   = (argc>4) ? atoi(argv[4]) : 0;

  HPEP_anchor_mode = anchor;
  HPEP_pool_mode   = pool;
  HPEP_init( Nt , Lmax );

  if( HPEP_nwilson()==0 && HPEP_npoly()==0 ){
    fprintf( stderr , "no tables found under %s -- run from the repository root\n" ,
             HPE_TABLE_DIR );
    return 1;
  }

  printf( "### anchor = %-9s  pool = %-13s  Nt = %d  Lmax = %d\n" ,
          anchor_name(anchor) , pool_name(pool) , Nt , Lmax );
  printf( "%-10s %5s %10s %10s %12s %10s %6s\n" ,
          "sector" , "len" , "MM/site" , "staples" , "pairs" , "+ReTr" , "reach" );

  for( int i=0 ; i<HPEP_nwilson() ; i++ ){
    char nm[16]; snprintf( nm , sizeof(nm) , "W%d" , HPEP_wilson_len(i) );
    const long long mm = HPEP_wilson_mm(i) , np = HPEP_wilson_npair(i);
    printf( "%-10s %5d %10lld %10.1f %12lld %10.1f %6d\n" , nm , HPEP_wilson_len(i) ,
            mm , HPEP_staples(mm) , np ,
            HPEP_staples(mm) + HPEP_trace_mm(np,false)/HPE_MM_PER_STAPLE ,
            HPEP_wilson_reach(i) );
  }
  for( int i=0 ; i<HPEP_npoly() ; i++ ){
    char nm[24];
    snprintf( nm , sizeof(nm) , "L%d(%d,%d)" , HPEP_poly_wind(i) , Nt , HPEP_poly_len(i) );
    const long long mm = HPEP_poly_mm(i) , np = HPEP_poly_npair(i);
    printf( "%-10s %5d %10lld %10.1f %12lld %10.1f %6d\n" , nm , HPEP_poly_len(i) ,
            mm , HPEP_staples(mm) , np ,
            HPEP_staples(mm) + HPEP_trace_mm(np,true)/HPE_MM_PER_STAPLE ,
            HPEP_poly_reach(i) );
  }
  printf( "pools = %d, total MM/site = %lld (%.1f staples), max reach = %d\n" ,
          HPEP_npool() , HPEP_total_mm() , HPEP_staples(HPEP_total_mm()) , HPEP_maxreach() );
  return 0;
}

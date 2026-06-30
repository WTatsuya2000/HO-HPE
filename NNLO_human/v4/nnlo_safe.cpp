#include <algorithm>
#include "nnlo_traj_list.cpp" 

int cyclic_symmetry(const std::vector<int>& traj)
{
    const int n = traj.size();
    int count = 0;
    for (int shift = 0; shift < n; ++shift) {
        bool same = true;
        for (int i = 0; i < n; ++i) {
            if (traj[i] != traj[(i + shift) % n]) {
                same = false;
                break;
            }
        }
        if (same) ++count;
    }
    return count;
}

bool cyclic_equal(const std::vector<int>& a,
                  const std::vector<int>& b)
{
    if (a.size() != b.size()) return false;
    std::vector<int> aa;
    aa.reserve(2 * a.size());
    aa.insert(aa.end(), a.begin(), a.end());
    aa.insert(aa.end(), a.begin(), a.end());
    return std::search(
        aa.begin(), aa.end(),
        b.begin(), b.end()
    ) != aa.end();
}

bool path_equal(const std::vector<int>& a,
                const std::vector<int>& b)
{
    if (cyclic_equal(a, b))
        return true;
    std::vector<int> rev;
    rev.reserve(a.size());
    for (auto it = a.rbegin(); it != a.rend(); ++it)
        rev.push_back(-*it);
    return cyclic_equal(rev, b);
}

std::vector<int> TrajConvert( const std::vector<int> traj ,
			      const int mulist[4] ,
			      const int signlist[4] )
{
  const int len = traj.size();
  std::vector<int> ans(len);
  for( int i=0 ; i<len ; i++ ){
    int n = std::abs(traj[i])-1;
    int s = (traj[i]>0) - (traj[i]<0);
    ans[i] = (mulist[n]+1) * s * signlist[n];
  }
  return ans;
}

std::vector<std::vector<int> > MakeTrajList( const std::vector<int> traj )
{
  static constexpr int P4[24][4] = {
    {0,1,2,3} , {0,1,3,2} , {0,2,1,3} , {0,2,3,1} , {0,3,1,2} , {0,3,2,1} ,
    {1,0,2,3} , {1,0,3,2} , {1,2,0,3} , {1,2,3,0} , {1,3,0,2} , {1,3,2,0} ,
    {2,0,1,3} , {2,0,3,1} , {2,1,0,3} , {2,1,3,0} , {2,3,0,1} , {2,3,1,0} ,
    {3,0,1,2} , {3,0,2,1} , {3,1,0,2} , {3,1,2,0} , {3,2,0,1} , {3,2,1,0}
  };
  static constexpr int slist4[16][4] = {
    {1,1,1,1},{1,1,1,-1},{1,1,-1,1},{1,1,-1,-1},{1,-1,1,1},{1,-1,1,-1},{1,-1,-1,1},{1,-1,-1,-1},
    {-1,1,1,1},{-1,1,1,-1},{-1,1,-1,1},{-1,1,-1,-1},{-1,-1,1,1},{-1,-1,1,-1},{-1,-1,-1,1},{-1,-1,-1,-1}
  };

  const int len = traj.size();
  std::vector<std::vector<int> > tlist;
  tlist.push_back( traj );
  for( const auto& m : P4 ){
    for( const auto& s : slist4 ){
      auto tmp = TrajConvert( traj , m , s );
      bool found = false;
      for( auto t : tlist )
	if( path_equal( t , tmp ) ){
	  found = true;
	  break;
	}
      if( !found )
	tlist.push_back( tmp );
    }
  }
  return tlist;
}  

std::vector<std::vector<int> > MakeTrajListP( const std::vector<int> traj )
{
  static constexpr int P3[6][4] = {
    {0,1,2,3},  {0,2,1,3},  {1,0,2,3},  {1,2,0,3},  {2,0,1,3},  {2,1,0,3}
  };
  static constexpr int slist3[8][4] = {
    {1,1,1,1},{1,1,-1,1},{1,-1,1,1},{1,-1,-1,1},
    {-1,1,1,1},{-1,1,-1,1},{-1,-1,1,1},{-1,-1,-1,1}
  };

  const int len = traj.size();
  std::vector<std::vector<int> > tlist;
  tlist.push_back( traj );
  for( const auto& m : P3 ){
    for( const auto& s : slist3 ){
      auto tmp = TrajConvert( traj , m , s );
      bool found = false;
      for( auto t : tlist )
	if( path_equal( t , tmp ) ){
	  found = true;
	  break;
	}
      if( !found )
	tlist.push_back( tmp );
    }
  }
  return tlist;
}  

double _Complex Safe_local1( const Gauge& V , const Site& x , const std::vector<int> traj )
{
  const int len = traj.size();
  double _Complex Uptr[9];
  Matrix<3,3> tmp( Uptr );
  double _Complex ans = 0.;
  //for(int y : traj) std::cout << y << ' ';
  //std::cout << '\n';

  for( int origin=0 ; origin<len ; origin++ ){
    auto xh = x;
    tmp = 1.;
    for( int n=0 ; n<len ; n++ ){
      const int index = ( origin + n ) % len;
      const int dir = traj[index];
      const int mu = std::abs(dir) % 4;//(1,2,3,4)-->(1,2,3,0)
      if( dir>0 ){
	tmp = tmp * V(xh,mu);
	xh = xh + mu;
      }
      else{
	xh = xh - mu;
	tmp = tmp * Hermite(V(xh,mu));
      }
    }
    ans += Trace(tmp) / (double)len;
  }
  return ans;
}

double _Complex TrajSum_safe( const Gauge& V , const vector<vector<int>> trajlist )
{
  double _Complex ans = 0.;
#pragma omp parallel reduction( + : ans )
  {
    Site x(V.lattice());
    ForAllSitesOMP(x){
      for( const auto& tr : trajlist )
        ans += Safe_local1( V , x , tr );
    }
  }
  mpi.add(ans);
  return -2. * ans / V.lattice().size();
}

std::array<double,25> NNLO_W_safe_decomp( const Gauge& V , bool disp=true )
{
  std::array<double,25> ans;
  for( int i=0 ; i<25 ; i++ ){
    ans[i] = __real__(TrajSum_safe( V , MakeTrajList(W8traj[i])) ) * W8D[i];
    ans[i] /= cyclic_symmetry(W8traj[i]);
    if(disp) mpicout << i+1 << " " << ans[i] << " / ";
  }
  if(disp) mpicout << endl;
  return ans;
}

double NNLO_W_safe( const Gauge& V , bool disp=true )
{
  double ans = 0.;
  for( int i=0 ; i<25 ; i++ ){
    double tmp = __real__(TrajSum_safe( V , MakeTrajList(W8traj[i])) );
    tmp *= W8D[i] / cyclic_symmetry(W8traj[i]);
    ans += tmp;
    if(disp) mpicout << i+1 << " " << tmp << " / ";
  }
  if(disp) mpicout << endl;
  return ans;
}

std::array<double _Complex,25> NNLO_P_safe_decomp( const Gauge& V , bool disp=true )
{
  std::array<double _Complex,25> ans;
  for( int i=0 ; i<25 ; i++ ){
    ans[i] = -TrajSum_safe( V , MakeTrajListP(P8Nt4traj[i]) ) * P8Nt4D[i];
    ans[i] /= cyclic_symmetry(P8Nt4traj[i]);
    if(disp) mpicout << i+1 << " " << ans[i] << " / ";
  }
  if(disp) mpicout << endl;
  return ans;
}

double _Complex NNLO_P_safe( const Gauge& V , bool disp=true )
{
  double _Complex ans;
  for( int i=0 ; i<25 ; i++ ){
    double _Complex tmp = -TrajSum_safe( V , MakeTrajListP(P8Nt4traj[i]) );
    tmp *= P8Nt4D[i] / cyclic_symmetry(P8Nt4traj[i]);
    ans += tmp;
    if(disp) mpicout << i+1 << " " << tmp << " / ";
  }
  if(disp) mpicout << endl;
  return ans;
}

void TrajList()
{
  mpicout << "W:" << endl;
  for( int i=0 ; i<25 ; i++ ){
    const auto& traj = W8traj[i];
    mpicout << i+1 << ": ";
    for(int t : traj)
      mpicout << t << ' ';
    mpicout << " : "
	    << MakeTrajList(traj).size() << " "
	    << cyclic_symmetry(traj) << " "
	    << W8D[i] << endl;
  }

  mpicout << "P:" << endl;
  for( int i=0 ; i<25 ; i++ ){
    const auto& traj = P8Nt4traj[i];
    mpicout << i+1 << ": ";
    for(int t : traj)
      mpicout << t << ' ';
    mpicout << " : "
	    << MakeTrajListP(traj).size() << " "
	    << cyclic_symmetry(traj) << " "
	    << P8Nt4D[i] << endl;
  }
}

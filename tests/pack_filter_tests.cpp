#include "PackDmtPayback.hpp"
#include "H1EdgeRisk.hpp"
#include <random>
#include <stdexcept>
static void require(bool b,const char* s){if(!b)throw std::runtime_error(s);}
static std::vector<std::vector<size_t>> membership(const std::vector<SelectedPV>& v){std::vector<std::vector<size_t>> o;for(auto& p:v)o.push_back(p.representatives);return o;}
int main(){
 std::mt19937 rng(20261023);size_t warning_cases=0,old_scans=0;
 for(size_t trial=0;trial<120;++trial){
  std::vector<std::vector<double>> points;std::vector<SelectedPV> supports;
  for(size_t group=0;group<2+trial%6;++group){
   std::vector<size_t> members;const size_t m=3+(trial+group)%4;
   for(size_t j=0;j<m;++j){members.push_back(points.size());points.push_back({double(group),double(j)*.02,double((trial/6)%3)*.003});}
   supports.push_back({members,.12});
  }
  std::shuffle(supports.begin(),supports.end(),rng);DistanceMatrix d(points);WindowState state(points.size());
  if(trial%3==0){auto pv=supports.back();supports.pop_back();std::unordered_set<size_t> absorbed(pv.representatives.begin(),pv.representatives.end());state.commitSelectedPVs({std::move(pv)},absorbed);}
  const auto before=state.activeLabels();auto risk=H1EdgeRisk::inspect(d,state,supports,.3,10.);
  if(!risk.warnings.empty())++warning_cases;
  auto drops=H1EdgeRisk::greedyDrops(risk,supports,state.totalLabelCount());auto fast=supports;std::vector<bool> kept(fast.size(),true);
  for(size_t id:drops)kept[id]=false;std::vector<SelectedPV> fast_kept;for(size_t i=0;i<fast.size();++i)if(kept[i])fast_kept.push_back(fast[i]);
  auto old=supports;std::vector<std::vector<size_t>> old_drops;auto current=risk;++old_scans;
  while(!current.warnings.empty()){
   std::vector<size_t> count(old.size(),0);for(auto& w:current.warnings)for(auto b:{w.lhs,w.rhs})if(current.blocks[b].proposed)++count.at(current.blocks[b].label-state.totalLabelCount());
   size_t id=old.size();double best=INFINITY;for(size_t i=0;i<old.size();++i)if(count[i]){const double c=double(old[i].representatives.size()-1)/count[i];if(c<best){best=c;id=i;}}
   require(id<old.size(),"old selector could not remove support");old_drops.push_back(old[id].representatives);old.erase(old.begin()+id);current=H1EdgeRisk::inspect(d,state,old,.3,10.);++old_scans;
  }
  require(membership(old)==membership(fast_kept),"one-scan retained order differs");
  require(old_drops.size()==drops.size(),"drop count differs");for(size_t i=0;i<drops.size();++i)require(old_drops[i]==supports[drops[i]].representatives,"drop order differs");
  require(state.activeLabels()==before,"graph filter mutated state");
 }
 // Third block provides fillers; pairwise warning deliberately remains conservative.
 std::vector<std::vector<double>> filler={{0,-.1},{0,0},{0,.1},{1,-.1},{1,0},{1,.1},{.5,0}};
 DistanceMatrix d(filler);WindowState state(7);std::vector<SelectedPV> proposed={{{0,1,2},.2},{{3,4,5},.2}};
 require(H1EdgeRisk::inspect(d,state,proposed,.3,1.1).warnings.size()==1,"third block changed pairwise detector");
 PackDmtPayback p;p.load(d,1.1,"measure","");require(p.candidates.size()<=8,"unbounded candidates");
 bool dup=false;require(PackDmtPayback::fullEndpoint(d,1.1,dup),"full endpoint not detected");
 require(PackDmtPayback::actualH(d,state,proposed,1.1,true)==PackDmtPayback::completeH(3),"actual partition count wrong");
 for(auto& c:p.candidates){require(c.q<1.1,"candidate at endpoint");require(c.Hend==PackDmtPayback::completeH(c.blocks),"full-endpoint shortcut differs");require(c.labels_reduced==7-c.blocks,"packing label accounting differs");}
 DistanceMatrix duplicate({{0,0},{0,0},{1,0}});require(!PackDmtPayback::fullEndpoint(duplicate,2,dup)&&dup,"duplicate contract not marked");
 require(warning_cases>0 && old_scans>120,"test did not exercise scan removal");
 std::cout<<"warning graph parity trials=120 warning_cases="<<warning_cases<<" old_scans="<<old_scans<<" fast_scans=120\n";
}

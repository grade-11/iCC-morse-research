#include "GeometryController.hpp"
#include "H1EdgeRisk.hpp"
#include <stdexcept>
#include <string>

static void require(bool condition,const std::string& message)
{if(!condition)throw std::runtime_error(message);}

int main()
{
    // Two disjoint bridges appear together, then mixed triangles join them.
    const std::vector<std::vector<double>> rectangle={{0,0},{0,.4},{2,0},{2,.4}};
    DistanceMatrix d(rectangle);WindowState state(4);
    const std::vector<SelectedPV> supports={{{0,1},.4},{{2,3},.4}};
    auto risk=H1EdgeRisk::inspect(d,state,supports,.6,3.);
    require(risk.warnings.size()==1,"missed disconnected cross-edge fiber");
    const auto& warning=risk.warnings.front();
    require(warning.epsilon==std::nextafter(2.,INFINITY),"strict epsilon event incorrect");
    require(warning.components==2 && warning.witness_edges.size()==2,"lost independent edge components");
    require(warning.resolved_at==std::nextafter(std::sqrt(4.16),INFINITY),"mixed triangle resolution incorrect");
    require(state.activeLabels()==std::vector<size_t>({0,1,2,3}) && state.pseudoVertices().empty(),"forecast mutated live state");
    require(H1EdgeRisk::inspect(d,state,supports,.6,2.).warnings.empty(),"included edge at excluded endpoint");
    require(H1EdgeRisk::inspect(d,state,supports,.6,std::nextafter(2.,INFINITY)).warnings.size()==1,"missed admitted final event");
    require(H1EdgeRisk::inspect(d,state,{{{0,1},.4}},.6,3.).warnings.empty(),"PV-singleton star is connected");

    // All four cross edges tie: inspect only the whole tie batch, never its prefixes.
    DistanceMatrix tied({{-.1,0,0},{.1,0,0},{0,-.1,2},{0,.1,2}});
    require(H1EdgeRisk::inspect(tied,state,{{{0,1},.2},{{2,3},.2}},.3,3.).warnings.empty(),"false alarm inside tie batch");

    auto padded=rectangle;
    for(size_t j=0;j<12;++j)padded.push_back({10.+.1*j,5.});
    DistanceMatrix padded_d(padded);WindowState padded_state(padded.size());
    auto padded_risk=H1EdgeRisk::inspect(padded_d,padded_state,supports,.6,20.);
    require(padded_risk.warnings.size()==1 && padded_risk.warnings.front().epsilon==warning.epsilon,"remote vertices changed local edge-fiber risk");

    auto scaled=rectangle;for(auto& p:scaled)for(auto& x:p)x*=8;
    DistanceMatrix scaled_d(scaled);
    const auto scaled_risk=H1EdgeRisk::inspect(scaled_d,state,{{{0,1},3.2},{{2,3},3.2}},4.8,24.);
    require(scaled_risk.warnings.size()==1 && scaled_risk.warnings.front().epsilon==warning.epsilon*8,"power-of-two units changed risk event");

    GeometryController controller;controller.load(padded_d,20.,.85,2);
    for(size_t i=1;i<controller.schedule.size();++i)require(controller.schedule[i-1]<controller.schedule[i],"schedule not strictly increasing");
    require(controller.schedule.size()<=3,"too many cuts");
    require(controller.snapshots.size()<=32 && controller.candidates.size()<=64,"unbounded feature grid");
    if(!controller.fallback)require(controller.schedule.back()==20.,"changed shared endpoint");
    std::cout<<"geometry and H1 risk regression tests passed\n";
}

#pragma once
#include <array>
#include <cstddef>
#include <cstdint>
#include <unordered_map>
#include <vector>
#include "DistanceMatrix.hpp"

// Experimental standard boundary-column reduction, NOT native Morse matching.
// H0/H1 only, fixed original labels, no quotient or G. Edge sorting is paid once;
// triangles are generated only when their last edge is admitted. Reduced pivot
// columns and live births persist across checkpoints; no old column is reduced
// again. Equal-weight edges precede triangles, both in colex order.
class StreamingF2
{
public:
    struct Pair { std::size_t dimension; double birth,death; bool operator==(const Pair&) const = default; };
    struct Counts { std::size_t new_edges=0,new_triangles=0,total_edges=0,total_triangles=0,pivots=0; };
    StreamingF2(const DistanceMatrix& distances,double cutoff);
    Counts advance(double epsilon);
    [[nodiscard]] std::vector<Pair> diagram(std::size_t dimension) const;
    [[nodiscard]] std::uint64_t stateHash() const;
private:
    using Id=std::size_t;
    struct Edge { Id a,b,bindex; double weight; };
    struct Triangle { std::array<Id,3> vertices; Id bindex; };
    struct Birth { std::size_t dimension; double weight; bool alive; };
    std::size_t n_,cursor_=0,triangle_count_=0;
    double cutoff_,epsilon_=-1.;
    std::vector<Edge> edges_;
    std::vector<std::vector<std::uint8_t>> adjacency_;
    std::vector<std::vector<Id>> edge_ids_;
    std::vector<Birth> births_;
    std::unordered_map<Id,std::vector<Id>> pivots_;
    std::vector<Pair> finite_;
    void reduce(std::vector<Id> column,std::size_t dimension,double weight,Id birth_id);
};

struct PipelineConfig;
void runStreamingF2(const DistanceMatrix& distances,const PipelineConfig& config);

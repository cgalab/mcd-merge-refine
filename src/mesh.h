#ifndef MESH_H
#define MESH_H

#include <set>
#include <boost/operators.hpp>
#include "gmobj.h"

class Vertex;
class Edge;
using VertexPtr = std::shared_ptr<Vertex>;
using EdgePtr = std::shared_ptr<Edge>;

class Vertex 
: private boost::equivalent<Vertex>
, boost::totally_ordered<Vertex> {
    class NeighborComp {
    public:
        NeighborComp(const GmPnt &center);

        bool operator()(const VertexPtr &lhs, const VertexPtr &rhs);

    private:
        GmPnt m_center;
    };

public:
    Vertex(int id = -1, const GmPnt &pnt = GmPnt(), bool bIsOnConvexHull = false);
    
    bool operator<(const Vertex &other) const;
    
    int id() const;
    const GmPnt & pnt() const;
    bool isOnConvexHull() const;
    const std::set<VertexPtr, NeighborComp> & neighbors() const;
    bool findNeighbor(VertexPtr &neighbor, const VertexPtr &vertex, bool bNext);
    
    void setIsOnConvexHull(bool bIsOnConvexHull);
    void addNeighbor(const VertexPtr &neighbor);
    void removeNeighbor(const VertexPtr &neighbor);
    
private:
    int m_id;
    GmPnt m_pnt;
    bool m_bIsOnConvexHull;
    std::set<VertexPtr, NeighborComp> m_neighbors;
};

class Edge
: private boost::equivalent<Edge>
, boost::totally_ordered<Edge> {
public:
    Edge();
    Edge(const VertexPtr &start, const VertexPtr &end, bool bIsOnConvexHull = false);
    
    bool operator<(const Edge &other) const;
    
    const VertexPtr & start() const;
    const VertexPtr & end() const;
    const GmSeg & seg() const;
    const std::pair<int, int> & id() const;
    bool isOnConvexHull() const;
    bool isConvex() const;
    std::vector<std::pair<int, int>> neighborIds() const;
    
    void setIsOnConvexHull(bool bIsOnConvexHull);
    
private:
    VertexPtr m_start;
    VertexPtr m_end;
    bool m_bIsOnConvexHull;
    GmSeg m_seg;
    std::pair<int, int> m_id;
};

class Face {
};

#endif /* MESH_H */


#include <iostream>
#include <assert.h>
#include <math.h>
#include <algorithm>
#include "mesh.h"

Vertex::NeighborComp::NeighborComp(const GmPnt &center) 
: m_center(center) {
}

bool Vertex::NeighborComp::operator()(const VertexPtr &lhs, const VertexPtr &rhs) {
    double alpha = m_center.angle(lhs->pnt()), beta = m_center.angle(rhs->pnt());
    
    if (!GmObj::epsEqual(alpha, beta)) {
        return alpha < beta;
    }
    
    return false;
}

Vertex::Vertex(int id, const GmPnt &pnt, bool bIsOnConvexHull)
: m_id(id)
, m_pnt(pnt)
, m_bIsOnConvexHull(bIsOnConvexHull)
, m_neighbors(std::set<VertexPtr, NeighborComp>(NeighborComp(pnt))) {
}

bool Vertex::operator<(const Vertex &other) const {
    return m_id < other.id();
}

int Vertex::id() const {
    return m_id;
}

const GmPnt & Vertex::pnt() const {
    return m_pnt;
}

bool Vertex::isOnConvexHull() const {
    return m_bIsOnConvexHull;
}

const std::set<VertexPtr, Vertex::NeighborComp> & Vertex::neighbors() const {
    return m_neighbors;
}

bool Vertex::findNeighbor(VertexPtr &neighbor, const VertexPtr &vertex, 
                          bool bNext) {
    const auto it = m_neighbors.find(vertex);

    if (it != m_neighbors.end()) {
        if (bNext) {
            if (std::next(it) != m_neighbors.end()) {
                neighbor = *std::next(it);
            } else {
                neighbor = *m_neighbors.begin();
            }
        } else {
            if (it != m_neighbors.begin()) {
                neighbor = *std::prev(it);
            } else {
                neighbor = *m_neighbors.rbegin();
            }
        }
        
        return true;
    }
    
    return false;
}

void Vertex::setIsOnConvexHull(bool bIsOnConvexHull) {
    m_bIsOnConvexHull = bIsOnConvexHull;
}

void Vertex::addNeighbor(const VertexPtr &neighbor) {
    m_neighbors.insert(neighbor);
}

void Vertex::removeNeighbor(const VertexPtr &neighbor) {
    if (m_neighbors.find(neighbor) != m_neighbors.end()) {
        m_neighbors.erase(m_neighbors.find(neighbor));
    }
}

Edge::Edge() {
}

Edge::Edge(const VertexPtr &start, const VertexPtr &end, bool bIsOnConvexHull)
: m_start(start)
, m_end(end)
, m_bIsOnConvexHull(bIsOnConvexHull) {
    assert(*m_start != *m_end);
    if (*m_start > *m_end) {
        std::swap(m_start, m_end);
    }
    
    m_seg = GmSeg(m_start->pnt(), m_end->pnt());
    m_id = std::make_pair(m_start->id(), m_end->id());
}

bool Edge::operator<(const Edge &other) const {
    if (*m_start != *other.start()) {
        return *m_start < *other.start();
    }
    
    return *m_end < *other.end();
}

const VertexPtr & Edge::start() const {
    return m_start;
}

const VertexPtr & Edge::end() const {
    return m_end;
}

const GmSeg & Edge::seg() const {
    return m_seg;
}

const std::pair<int, int> & Edge::id() const {
    return m_id;
}

bool Edge::isOnConvexHull() const {
    return m_bIsOnConvexHull;
}

bool Edge::isConvex() const {
    const auto fIsConvex
            = [&](const VertexPtr &v, const VertexPtr & w) {
                const auto p = v->pnt();
                VertexPtr n1, n2;
                bool bOk1 = v->findNeighbor(n1, w, false),
                        bOk2 = v->findNeighbor(n2, w, true);

                if (bOk1 && bOk2) {
                    const auto q = n1->pnt(), r = n2->pnt();
                    double alpha = p.angle(q), beta = p.angle(r),
                            diff = alpha < beta ? (beta - alpha) : (2. * M_PI - alpha + beta);
                    
                    return diff < M_PI;
                }

                return false;
            };

    return !m_bIsOnConvexHull && fIsConvex(m_start, m_end) && fIsConvex(m_end, m_start);
}

std::vector<std::pair<int, int>> Edge::neighborIds() const {
    std::vector<std::pair<int, int>> neighborIds;

    const auto fFindNeighbors
            = [&](bool bIsStart) {
                const auto v = bIsStart ? m_start : m_end,
                        w = bIsStart ? m_end : m_start;
                VertexPtr n1, n2;
                bool bOk1 = v->findNeighbor(n1, w, false),
                        bOk2 = v->findNeighbor(n2, w, true);

                if (bOk1) {
                    neighborIds.push_back(std::minmax({v->id(), n1->id()}));
                }

                if (bOk2) {
                    neighborIds.push_back(std::minmax({v->id(), n2->id()}));
                }
            };

    fFindNeighbors(false);
    fFindNeighbors(true);
            
    return neighborIds;
}

void Edge::setIsOnConvexHull(bool bIsOnConvexHull) {
    m_bIsOnConvexHull = bIsOnConvexHull;
}

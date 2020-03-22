#ifndef CONVEX_PART_H
#define CONVEX_PART_H

#include <set>
#include <stack>
#ifdef ENABLE_VIEW
#include <QWidget>
#include <QPainter>
#include <CGAL/Qt/GraphicsItem.h>
#include <unordered_map>
#endif
#include "gmobj.h"
#include "file_io.h"

enum class BaseStrategy {
    Delaunay,
};

enum class RefinementStrategy {
    EdgeRemoval,
};

class ConvexDecomp {
public:
    ConvexDecomp(const std::string &inFilePath,
                 const std::string &outFilePath = "",
                 int status_fd = -1,
                 BaseStrategy base = BaseStrategy::Delaunay,
                 RefinementStrategy ref = RefinementStrategy::EdgeRemoval);
    
    const std::unordered_map<int, VertexPtr> & verts() const {
        return m_verts;
    }

    const std::unordered_map<std::pair<int, int>, EdgePtr, HashPair> & edges() const {
        return m_edges;
    }
    
private:
    void compDelaunay();
    void removeEdges();
    
    std::unordered_map<int, VertexPtr> m_verts;
    std::unordered_map<int, VertexPtr> m_convexVerts;
    std::unordered_map<std::pair<int, int>, EdgePtr, HashPair> m_edges;
    std::unordered_map<std::pair<int, int>, EdgePtr, HashPair> m_convexEdges;
    int m_nFaces = 0;
    unsigned m_seed = 0;
};

#ifdef ENABLE_VIEW

class ConvexDecompGraphicsItem : public CGAL::Qt::GraphicsItem {
    Q_OBJECT

public:
    ConvexDecompGraphicsItem(const ConvexDecomp &cd);

    QRectF boundingRect() const override;

    void paint(QPainter *painter, const QStyleOptionGraphicsItem *option,
               QWidget *widget = nullptr) override;

public slots:
    void modelChanged() override;
    void onPrevEv();
    void onNextEv();
    
private:
    ConvexDecomp m_cd;
    int m_index = 0;
    double m_scale = 0.;
    QRectF m_boundingRect;
};

#endif

#endif /* CONVEX_PART_H */

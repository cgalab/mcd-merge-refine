#include <math.h>
#include <time.h>
#include <iostream>
#include <sstream>
#include <algorithm>
#include <assert.h>
#include <chrono>
#include <random>
#include <map>
#include "decomp.h"
extern "C" {
#include <triangle.h>
}


#include "../gitversion.h"

/* timing stuff */
#include <sys/time.h>
#include <sys/resource.h>
#include <errno.h>
#include <string.h>

static double get_current_rtime(void) {
  struct rusage usage;
  if (getrusage(RUSAGE_SELF, &usage) < 0) {
    fprintf(stderr, "getrusage() failed: %s\n", strerror(errno));
    exit(1);
  }
  return usage.ru_utime.tv_sec + (double)usage.ru_utime.tv_usec/1e6;
}

/* we are too fast for sane use of getrusage.  Also try a wallclock timer */
#include <chrono>
typedef std::chrono::high_resolution_clock Clock;

static long get_maxrss(void) {
  struct rusage usage;
  if (getrusage(RUSAGE_SELF, &usage) < 0) {
    fprintf(stderr, "getrusage() failed: %s\n", strerror(errno));
    exit(1);
  }
  return usage.ru_maxrss;
}


ConvexDecomp::ConvexDecomp(const std::string &inFilePath,
                           const std::string &outFilePath,
                           int status_fd,
                           BaseStrategy base, RefinementStrategy ref)
: m_seed(std::random_device("/dev/urandom")()) {
    std::cout << "random_seed:     " << m_seed << std::endl;
    
    const auto fr = FileReader(inFilePath);
    m_verts = fr.vertices();

    double start_rtime = get_current_rtime();
    auto start_hirestime = Clock::now();

    /*std::chrono::high_resolution_clock::time_point t1 =
            std::chrono::high_resolution_clock::now();*/
    
    if (m_verts.size() > 2) {
        switch (base) {
            case BaseStrategy::Delaunay:
            {
                compDelaunay();
                break;
            }
            default:
            {
                break;
            }
        }

        switch (ref) {
            case RefinementStrategy::EdgeRemoval:
            {
                removeEdges();
                break;
            }
            default:
            {
                break;
            }
        }
    }

    /*std::chrono::high_resolution_clock::time_point t2
            = std::chrono::high_resolution_clock::now();*/

    /*long long duration
            = std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1).count();*/
    
    size_t nFaces = (m_edges.size() - m_verts.size() + 1);
    std::cout << "num_cvx_areas:   " << nFaces << "\n";

    double end_rtime = get_current_rtime();
    auto end_hirestime = Clock::now();
    long rmem = get_maxrss();
    if (status_fd >= 0) {
       FILE *status = fdopen(status_fd, "a");
       if (!status) {
          fprintf(stderr, "Cannot open status FD %d: %s\n", status_fd, strerror(errno));
          exit(-1);
       }

       fprintf(status, "[STATUS] VERSION: %s\n", GITVERSION);
       fprintf(status, "[STATUS] GENERATOR: mcd-ref-merge\n");
       fprintf(status, "[STATUS] INPUT_SIZE: %ld\n", m_verts.size());
       fprintf(status, "[STATUS] CPUTIME: %.6lf\n", end_rtime - start_rtime);
       fprintf(status, "[STATUS] WALLTIME: %.9lf\n",
         double(std::chrono::duration_cast<std::chrono::nanoseconds>(end_hirestime - start_hirestime).count())/1e9);
       fprintf(status, "[STATUS] MAXRSS: %ld\n", rmem);
    }


    if (!outFilePath.empty()) {
        const auto name = fr.instanceName();
        FileWriter(outFilePath, name, m_verts, m_edges)();
    }
}

void ConvexDecomp::compDelaunay() {
    struct triangulateio in, out;

    /* Define input points. */
    in.numberofpoints = m_verts.size();
    in.numberofpointattributes = 0;
    in.pointlist = (REAL *) malloc(in.numberofpoints * 2 * sizeof(REAL));
    for (size_t i = 0, j = 0; i < m_verts.size(); i++, j += 2) {
        in.pointlist[j] = m_verts.at(i)->pnt().x();
        in.pointlist[j + 1] = m_verts.at(i)->pnt().y();
    }

    in.pointattributelist = (REAL *) malloc(in.numberofpoints *
                                            in.numberofpointattributes *
                                            sizeof(REAL));
    in.pointmarkerlist = (int *) malloc(in.numberofpoints * sizeof(int));
    in.numberofsegments = 0;
    in.numberofholes = 0;
    in.numberofregions = 1;
    in.regionlist = (REAL *) malloc(in.numberofregions * 4 * sizeof(REAL));

    out.pointlist = (REAL *) NULL; /* Not needed if -N switch used. */
    /* Not needed if -N switch used or number of point attributes is zero: */
    out.pointattributelist = (REAL *) NULL;
    out.pointmarkerlist = (int *) NULL; /* Not needed if -N or -B switch used. */
    out.trianglelist = (int *) NULL; /* Not needed if -E switch used. */
    /* Not needed if -E switch used or number of triangle attributes is zero: */
    out.triangleattributelist = (REAL *) NULL;
    out.neighborlist = (int *) NULL; /* Needed only if -n switch used. */
    /* Needed only if segments are output (-p or -c) and -P not used: */
    out.segmentlist = (int *) NULL;
    /* Needed only if segments are output (-p or -c) and -P and -B not used: */
    out.segmentmarkerlist = (int *) NULL;
    out.edgelist = (int *) NULL; /* Needed only if -e switch used. */
    out.edgemarkerlist = (int *) NULL; /* Needed if -e used and -B not used. */

    /* Triangulate the points.  Switches are chosen to read and write a  */
    /*   PSLG (p), preserve the convex hull (c), number everything from  */
    /*   zero (z), assign a regional attribute to each element (A), and  */
    /*   produce an edge list (e), a Voronoi diagram (v), and a triangle */
    /*   neighbor list (n).                                              */
    char options[] = "pczen";
    triangulate(options, &in, &out, (struct triangulateio *) NULL);
    
    trifree(in.pointlist);
    trifree(in.pointattributelist);
    trifree(in.pointmarkerlist);
    trifree(in.regionlist);
    
    for (int i = 0; i < 2 * out.numberofedges; i += 2) {
        int id1 = out.edgelist[i], id2 = out.edgelist[i + 1];
        const auto edge = std::make_shared<Edge>(m_verts.at(id1), m_verts.at(id2));
        
        m_edges[edge->id()] = edge;
        m_verts[id1]->addNeighbor(m_verts.at(id2));
        m_verts[id2]->addNeighbor(m_verts.at(id1));
    }

    std::map<std::pair<int, int>, int> counts;
    for (int i = 0; i < 3 * out.numberoftriangles; i += 3) {
        int id1 = out.trianglelist[i], id2 = out.trianglelist[i + 1],
                id3 = out.trianglelist[i + 2];

        counts[std::minmax(id1, id2)]++;
        counts[std::minmax(id1, id3)]++;
        counts[std::minmax(id2, id3)]++;
    }
    
    trifree(out.pointlist);
    trifree(out.pointattributelist);
    trifree(out.pointmarkerlist);
    trifree(out.trianglelist);
    trifree(out.triangleattributelist);
    trifree(out.neighborlist);
    trifree(out.segmentlist);
    trifree(out.segmentmarkerlist);
    trifree(out.edgelist);
    trifree(out.edgemarkerlist);
    
    for (const auto &count : counts) {
        auto edge = m_edges.at(count.first);
        if (count.second < 2) {
            edge->setIsOnConvexHull(true);
            edge->start()->setIsOnConvexHull(true);
            edge->end()->setIsOnConvexHull(true);
        }
    }
    
    for (const auto &entry : m_edges) {
        const auto edge = entry.second;
        
        if (edge->isConvex()) {
            m_convexEdges[edge->id()] = edge;
        }
    }
    
    /*for (const auto &entry : m_verts) {
        const auto v = entry.second;
        
        bool bIsConvex = true;
        for (const auto w : v->neighbors()) {
            const auto edge = m_edges.at(std::minmax(v->id(), w->id()));
            
            if (edge->isConvex()) {
                bIsConvex = false;
                break;
            }
        }
        
        if (bIsConvex) {
            m_convexVerts[v->id()] = v;
        }
    }*/
}

void ConvexDecomp::removeEdges() {
    std::cout << "Removing Edges ...\n";
    const auto fRemoveConvexEdge
            = [&](const EdgePtr & edge) {
                assert(edge->isConvex());
                const auto v = edge->start(), w = edge->end();

                const auto neighborIds = edge->neighborIds();
                v->removeNeighbor(w);
                w->removeNeighbor(v);

                if (m_convexEdges.find(edge->id()) != m_convexEdges.end()) {
                    m_convexEdges.erase(m_convexEdges.find(edge->id()));
                }

                if (m_edges.find(edge->id()) != m_edges.end()) {
                    m_edges.erase(m_edges.find(edge->id()));
                }

                for (const auto &neighborId : neighborIds) {
                    if (m_convexEdges.find(neighborId) != m_convexEdges.end()) {
                        const auto edge = m_convexEdges.at(neighborId);

                        if (!edge->isConvex()) {
                            m_convexEdges.erase(m_convexEdges.find(edge->id()));
                        }
                    }
                }
            };

    std::vector<std::pair<int, int>> edgeIds;
    for (const auto &entry : m_convexEdges) {
        const auto edgeId = entry.first;
        edgeIds.push_back(edgeId);
    }

    srand48(m_seed);
    std::random_shuffle(edgeIds.begin(), edgeIds.end(),
                        [](int i) {
                            return lrand48() % i;
                        });
    for (const auto &edgeId : edgeIds) {
        if (m_convexEdges.find(edgeId) != m_convexEdges.end()) {
            const auto edge = m_convexEdges.at(edgeId);
            fRemoveConvexEdge(edge);
        }
    }
            
    /*srand(time(NULL));
    while (!m_convexEdges.empty()) {
        size_t size = m_convexEdges.size();
        int steps = rand() % size;

        auto it = m_convexEdges.begin();
        std::advance(it, steps);
        const auto edge = it->second;
        fRemoveConvexEdge(edge);
    }*/
}

#ifdef ENABLE_VIEW

ConvexDecompGraphicsItem::ConvexDecompGraphicsItem(const ConvexDecomp &cd)
: m_cd(cd) {
    double width = 0., height = 0.;
    for (const auto &entry : m_cd.verts()) {
        const auto pnt = entry.second->pnt();
        
        if (pnt.x() > width) {
            width = pnt.x();
        }
        
        if (pnt.y() > height) {
            height = pnt.y();
        }
    }
    
    m_scale = 1000. / std::max(width, height);
    m_boundingRect = QRectF(0., 0., width, height);
}

QRectF ConvexDecompGraphicsItem::boundingRect() const {
    return m_boundingRect;
}

void ConvexDecompGraphicsItem::paint(QPainter *painter,
                                     const QStyleOptionGraphicsItem*,
                                     QWidget*) {
    painter->setPen(QPen(Qt::black, 1.));
    painter->setBrush(Qt::black);
    
    for (const auto &val : m_cd.edges()) {
        const auto edge = val.second;
        const auto from = edge->start(), to = edge->end();

        painter->drawLine(QPointF(from->pnt().x() * m_scale, from->pnt().y() * m_scale),
                          QPointF(to->pnt().x() * m_scale, to->pnt().y() * m_scale));
    }
    
    for (const auto &entry : m_cd.verts()) {
        const auto v = entry.second;
        double x = v->pnt().x(), y = v->pnt().y();
        painter->drawEllipse(QPointF(x, y) * m_scale, 1., 1.);
    }
}

void ConvexDecompGraphicsItem::modelChanged() {
    update();
}

void ConvexDecompGraphicsItem::onPrevEv() {
    m_index--;
    update();
}

void ConvexDecompGraphicsItem::onNextEv() {
    m_index++;
    update();
}

#include "moc_decomp.cpp"

#endif

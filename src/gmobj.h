#ifndef GMOBJ_H
#define GMOBJ_H

#include <vector>
#include <memory>
#include <limits>

class GmPnt;
class GmMatrix;
class GmRay;
class GmArc;

const double EPS_1 = 1.e-12;
const double EPS_2 = 1.e-9;
const double INFTY = 1.e4;
const int MATR_DIM = 4;

using Matr4 = std::array<std::array<double, MATR_DIM>, MATR_DIM>;

enum class GmType {
    None, Pnt, Line, Seg, Ray, Arc, Circ, Chain
};

enum class Side {
    None, Left, Right
};

class GmObj {
public:
    friend std::ostream & operator<<(std::ostream &os, GmObj *obj);
    
    GmObj() {}
    virtual ~GmObj() = default;
    
    GmObj(const GmObj &obj) = default;
    GmObj & operator=(const GmObj &obj) = default;
    GmObj(GmObj &&obj) = default;
    GmObj & operator=(GmObj &&obj) = default;

    virtual GmType type() const {
        return GmType::None;
    }
    
    virtual int isect(std::vector<GmPnt> &, const GmObj &) const {
        return -1;
    }
    
    virtual Side side(const GmPnt &) const {
        return Side::None;
    }
    
    static bool epsEqual(double lhs, double rhs, double eps = EPS_1);
    static bool epsLessEqual(double lhs, double rhs, double eps = EPS_1);
};

class DoubleComp {
public:
    bool operator()(double lhs, double rhs) const {
        if (!GmObj::epsEqual(lhs, rhs)) {
            return lhs < rhs;
        }
        
        return false;
    }
};

class GmPnt : public GmObj {
public:
    friend GmPnt operator*(double v, const GmPnt &pnt);
    friend GmPnt operator*(const GmPnt &pnt, double v);
    friend GmPnt operator*(const GmMatrix &matrix, const GmPnt &pnt);
    friend GmPnt operator*(const GmPnt &pnt, const GmMatrix &matrix);
    friend GmPnt operator/(const GmPnt &pnt, double v);

    GmPnt(double x = 0., double y = 0., double z = 0.);

    bool operator<(const GmPnt &other) const;
    bool operator<=(const GmPnt &other) const;
    bool operator>(const GmPnt &other) const;
    bool operator>=(const GmPnt &other) const;
    bool operator==(const GmPnt &other) const;
    bool operator!=(const GmPnt &other) const;
    
    GmPnt operator+(const GmPnt &pnt) const;
    GmPnt operator-(const GmPnt &pnt) const;

    double x() const {
        return m_x;
    }

    double y() const {
        return m_y;
    }

    double z() const {
        return m_z;
    }

    GmType type() const override {
        return GmType::Pnt;
    }

    double dist(const GmPnt &other) const;
    double cross(const GmPnt &other) const; // TODO: Do not use for 3D points!
    double dot(const GmPnt &other) const;
    double angle(const GmPnt &other) const;
    int isect(std::vector<GmPnt> &pnts, const GmObj &obj) const override;
    std::string toString() const;

private:
    double m_x;
    double m_y;
    double m_z;
};

class GmLine : public GmObj {
public:
    GmLine(double m, double  b);
    GmLine(const GmPnt &start = GmPnt(), const GmPnt &end = GmPnt(1., 0.));
    
    virtual ~GmLine() {
    }
    
    const GmPnt & start() const {
        return m_start; 
    }
    
    const GmPnt & end() const { 
        return m_end; 
    }
    
    double m() const { 
        return m_m; 
    }
    
    double b() const { 
        return m_b; 
    }
    
    bool isVertical() const {
        return m_bIsVertical;
    }
    
    double yAt(double x) const {
        return m_m * x + m_b;
    }
    
    double xAt(double y) const {
        return (y - m_b) / m_m;
    }
    
    GmType type() const override {
        return GmType::Line;
    }
    
    int isect(std::vector<GmPnt> &pnts, const GmObj &obj) const override;
    virtual bool includes(const GmPnt &pnt) const;
    
    GmLine normal(const GmPnt &pnt) const;
    GmPnt project(const GmPnt &pnt) const;
    GmPnt reflect(const GmPnt &pnt) const;
    Side side(const GmPnt &pnt) const override;
    bool isLeft(const GmPnt &pnt) const;
    std::pair<GmRay, GmRay> splitAt(const GmPnt &pnt) const;
    std::vector<GmPnt> filterPnts(const std::vector<GmPnt> &pnts) const;
    
private:
    GmPnt m_start;
    GmPnt m_end;
    double m_m;
    double m_b;
    bool m_bIsVertical = false;
};

class GmSeg : public GmLine {
public:
    GmSeg(const GmPnt &start = GmPnt(), const GmPnt &end = GmPnt(1., 0.));
    
    double length() const {
        return start().dist(end());
    }
    
    GmType type() const override {
        return GmType::Seg;
    }
    
    GmLine toGmLine() const {
        return GmLine(start(), end());
    }
    
    bool includes(const GmPnt &pnt) const override;
    int isect(std::vector<GmPnt> &pnts, const GmObj &obj) const override;
    GmPnt midPnt() const;
};

class GmRay : public GmLine {
public:
    GmRay(const GmPnt &start = GmPnt(), const GmPnt &direction = GmPnt(1., 0.));
    
    GmType type() const override {
        return GmType::Ray;
    }
    
    bool includes(const GmPnt &pnt) const override;
    int isect(std::vector<GmPnt> &pnts, const GmObj &obj) const override;
};

class GmCirc : public GmObj {
public:
    GmCirc(const GmPnt &center = GmPnt(), double radius = 0., bool ccw = true);
    GmCirc(const GmPnt &p1, const GmPnt &p2, bool ccw = true);
    GmCirc(const GmPnt &p1, const GmPnt &p2, const GmPnt &p3, bool ccw = true);
    
    GmType type() const override {
        return GmType::Circ;
    }
    
    const GmPnt & center() const {
        return m_center;
    }
    
    double radius() const {
        return m_radius;
    }
    bool ccw() const {
        return m_ccw;
    }
    
    bool contains(const GmCirc &other) const;
    bool contains(const GmPnt &pnt) const;
    virtual bool includes(const GmPnt &pnt) const;
    int isect(std::vector<GmPnt> &pnts, const GmObj &obj) const override;
    Side side(const GmPnt &pnt) const override;

    // FIRST contains the left and SECOND the right arc relative to the 
    // straight-line segment defined by START and END.
    std::pair<GmArc, GmArc> split(const GmPnt &start, const GmPnt &end);
    
private:
    int doIsect(std::vector<GmPnt> &pnts, const GmLine &line) const;
    int doIsect(std::vector<GmPnt> &pnts, const GmCirc &other) const;
    
    GmPnt m_center;
    double m_radius;
    bool m_ccw;
};

class GmArc : public GmCirc { // TODO: Arcs are always oriented counter-clockwise.
public:
    GmArc(const GmPnt &center = GmPnt(), double radius = 0.,
          double from = 0., double to = 0.); // from and two are in radians!
    GmArc(const GmPnt &center, const GmPnt &start, const GmPnt &end);

    GmType type() const override {
        return GmType::Arc;
    }
            
    const GmPnt & start() const {
        return m_start;
    }

    const GmPnt & end() const {
        return m_end;
    }
    
    double from() const {
        return m_from;
    }
    
    double to() const {
        return m_to;
    }
    
    bool includes(const GmPnt &pnt) const override;
    int isect(std::vector<GmPnt> &pnts, const GmObj &obj) const override;
            
private:
    GmPnt m_start;
    GmPnt m_end;
    double m_from;
    double m_to;
};

class GmChain : public GmObj {
public:
    GmChain(const std::vector<GmPnt> &pnts = std::vector<GmPnt>());
    
    const std::vector<GmPnt> & pnts() const {
        return m_pnts;
    }
    
    void addPnt(const GmPnt &pnt) {
        m_pnts.push_back(pnt);
    }

    GmType type() const override {
        return GmType::Chain;
    }
    
    int isect(std::vector<GmPnt> &pnts, const GmObj &obj) const override;
    
private:
    std::vector<GmPnt> m_pnts;
};

class GmMatrix {
public:
    GmMatrix(const Matr4 &matrix = {});
    
    virtual ~GmMatrix() {}
    
    GmMatrix operator*(const GmMatrix &other) const;
    
    double at(int i, int j) const;
    void set(int i, int j, double value);
    GmMatrix inv() const;
    double det() const;
    
private:
    Matr4 m_matrix;
};

class GmTransl : public GmMatrix {
public:
    GmTransl(const GmPnt &pnt = GmPnt());
};

class GmRot : public GmMatrix { // TODO: Rotates exclusively around the z-axis!
public:
    GmRot(double rad = 0., const GmPnt &pnt = GmPnt());
};

#endif /* GMOBJ_H */

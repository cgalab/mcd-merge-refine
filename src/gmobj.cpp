#include <complex>
#include <algorithm>
#include <assert.h>
#include <iostream>
#include "gmobj.h"

std::ostream & operator<<(std::ostream &os, GmObj *obj) {
    switch (obj->type()) {
        case GmType::Pnt:
        {
            const auto pnt = static_cast<const GmPnt *>(obj);
            os << pnt->toString();
            break;
        }
        case GmType::Line:
        case GmType::Ray:
        case GmType::Seg:
        {
            const auto line = static_cast<const GmLine *>(obj);
            os << "(" << static_cast<int>(obj->type()) << ")" 
                    << line->start().toString() << " -> " 
                    << line->end().toString();
            break;
        }
        case GmType::Circ:
        {
            const auto circ = static_cast<const GmCirc *> (obj);
            os << "c: " << circ->center().toString()
                    << " r: " << circ->radius();
            break;
        }
        case GmType::Arc:
        {
            break;
        }
        case GmType::None:
        case GmType::Chain:
        default:
        {
            break;
        }
    }
    
    return os;
}

bool GmObj::epsEqual(double lhs, double rhs, double eps) {
    return std::abs(lhs - rhs) < eps;
}

bool GmObj::epsLessEqual(double lhs, double rhs, double eps) {
    return lhs < rhs || GmObj::epsEqual(lhs, rhs, eps);
}

GmPnt::GmPnt(double x, double y, double z)
: m_x{x}
, m_y{y}
, m_z{z}
{
}

bool GmPnt::operator<(const GmPnt &other) const {
    if (!GmObj::epsEqual(x(), other.x(), EPS_2)) {
        return x() < other.x();
    } else if (!GmObj::epsEqual(y(), other.y(), EPS_2)) {
        return y() < other.y();
    } else if (!GmObj::epsEqual(z(), other.z(), EPS_2)) {
        return z() < other.z();
    }
    
    return false;
}

bool GmPnt::operator==(const GmPnt &other) const {
    return !(*this < other) && !(other < *this);
}

bool GmPnt::operator!=(const GmPnt &other) const {
    return !(*this == other);
}

bool GmPnt::operator>=(const GmPnt &other) const {
    return !(*this < other);
}

bool GmPnt::operator<=(const GmPnt &other) const {
    return (*this < other) || (*this == other);
}

bool GmPnt::operator>(const GmPnt &other) const {
    return !(*this <= other);
}

GmPnt GmPnt::operator+(const GmPnt &pnt) const {
    return GmPnt(x() + pnt.x(), y() + pnt.y(), z() + pnt.z());
}

GmPnt GmPnt::operator-(const GmPnt &pnt) const {
    return GmPnt(x() - pnt.x(), y() - pnt.y(), z() - pnt.z());
}

GmPnt operator*(double v, const GmPnt &pnt) {
    return GmPnt(v * pnt.x(), v * pnt.y(), v * pnt.z());
}

GmPnt operator*(const GmPnt &pnt, double v) {
    return operator*(v, pnt);
}

GmPnt operator*(const GmMatrix &matrix, const GmPnt &pnt) {
    return operator*(pnt, matrix);
}

GmPnt operator*(const GmPnt &pnt, const GmMatrix &matrix) {
    const auto a = std::array<double, MATR_DIM>{pnt.x(), pnt.y(), pnt.z(), 1.};
    std::array<double, MATR_DIM> c = {};

    for (int i = 0; i < MATR_DIM; i++) {
        for (int j = 0; j < MATR_DIM; j++) {
            c[i] += matrix.at(i, j) * a[j];
        }
    }

    return GmPnt(c[0], c[1], c[2]);
}

GmPnt operator/(const GmPnt &pnt, double v) {
    return GmPnt(pnt.x() * 1. / v, pnt.y() * 1. / v, pnt.z() * 1. / v);
}

double GmPnt::dist(const GmPnt &other) const {
    return std::sqrt(std::pow(x() - other.x(), 2.) + std::pow(y() - other.y(), 2.)
                     + std::pow(z() - other.z(), 2.));
}

double GmPnt::cross(const GmPnt &other) const {
    return x() * other.y() - y() * other.x();
}

double GmPnt::dot(const GmPnt &other) const {
    return x() * other.x() + y() * other.y() + z() * other.z();
}

double GmPnt::angle(const GmPnt &other) const {
    const auto p = other - *this;
    return atan2(p.y(), p.x()) + M_PI;
}

int GmPnt::isect(std::vector<GmPnt> &, const GmObj &) const {
    return -1;
}

std::string GmPnt::toString() const {
    std::stringstream strs;
    strs << "(" << m_x << ", " << m_y << ")";
    return strs.str();
}

GmLine::GmLine(double m, double b)
: m_m{m}
, m_b{b}
{
    double x1 = 0., x2 = 1.;
    m_start = GmPnt(x1, yAt(x1));
    m_end = GmPnt(x2, yAt(x2));
}

GmLine::GmLine(const GmPnt &start, const GmPnt &end)
: m_start{start}
, m_end{end}
, m_m{(m_end.y() - m_start.y()) / (m_end.x() - m_start.x())}
, m_b(m_start.y() - m() * m_start.x()) {
    assert(m_start != m_end);
    if (GmObj::epsEqual(m_end.x(), m_start.x())) {
        m_bIsVertical = true;
    }
}

int GmLine::isect(std::vector<GmPnt> &pnts, const GmObj &obj) const {
    switch (obj.type()) {
        case GmType::Ray:
        {
            const auto ray = static_cast<const GmRay*> (&obj);
            const auto other = GmLine(ray->start(), ray->end());
            
            std::vector<GmPnt> temps;
            isect(temps, other);
            pnts = ray->filterPnts(temps);
            
            return pnts.size();
        }
        case GmType::Seg:
        {
            const auto seg = static_cast<const GmSeg*> (&obj);
            const auto other = seg->toGmLine();
            
            std::vector<GmPnt> temps;
            isect(temps, other);
            pnts = seg->filterPnts(temps);
            
            return pnts.size();
        }
        case GmType::Line:
        {
            const auto other = static_cast<const GmLine *> (&obj);
            
            if (isVertical() && !other->isVertical()) {
                double x = m_start.x(),
                        y = other->yAt(m_start.x());
                pnts.push_back(GmPnt(x, y));
                return 1;
            } else if(!isVertical() && other->isVertical()) {
                double x = other->start().x(),
                        y = yAt(other->start().x());
                pnts.push_back(GmPnt(x, y));
                return 1;
            }

            double a = m(), c = b(), b = other->m(), d = other->b(),
                    x = (d - c) / (a - b), y = a * x + c;
            if (GmObj::epsEqual(a, b)) {
                return 0;
            } else {
                pnts.push_back(GmPnt(x, y));
                return 1;
            }
        }
        case GmType::Circ:
        case GmType::Arc:
        case GmType::Pnt:
        {
            return obj.isect(pnts, *this);
        }
        case GmType::Chain:
        case GmType::None:
        default:
        {
            break;
        }
    }

    return -1;
}

bool GmLine::includes(const GmPnt &pnt) const {
    if (m_bIsVertical) {
        return GmObj::epsEqual(pnt.x(), m_start.x());
    }
    
    return GmObj::epsEqual(pnt.y(), ((m_m * pnt.x())) + m_b);
}

GmLine GmLine::normal(const GmPnt &pnt) const {
    const auto vec1 = m_end - m_start,
            vec2 = pnt + GmPnt(vec1.y(), -vec1.x());

    return GmLine(pnt, vec2);
}

GmPnt GmLine::project(const GmPnt &pnt) const {
    if (!std::isinf(m_m)) {
        double x = (m_m * pnt.y() + pnt.x() - m_m * m_b) / (std::pow(m_m, 2.) + 1),
                y = (std::pow(m_m, 2.) * pnt.y() + m_m * pnt.x() + m_b) / (std::pow(m_m, 2.) + 1);
        return GmPnt(x, y);
    } else {
        return GmPnt(m_start.x(), pnt.y());
    }
}

GmPnt GmLine::reflect(const GmPnt &pnt) const {
    if (includes(pnt)) {
        return pnt;
    }
    
    return 2. * project(pnt) - pnt;
}

Side GmLine::side(const GmPnt &pnt) const {
    const auto result = (pnt.x() - m_start.x())*(m_end.y() - m_start.y())
            - (pnt.y() - m_start.y())*(m_end.x() - m_start.x());

    if (GmObj::epsEqual(result, 0.)) {
        return Side::None;
    } else if (result < 0.) {
        return Side::Left;
    }

    return Side::Right;
}

bool GmLine::isLeft(const GmPnt &pnt) const {
    return (side(pnt) == Side::Left);
}

std::pair<GmRay, GmRay> GmLine::splitAt(const GmPnt &pnt) const {
    double x = pnt.x(), y = pnt.y();
    return m_bIsVertical ? std::make_pair(GmRay(pnt, GmPnt(x, y + 1.)),
                                          GmRay(pnt, GmPnt(x, y - 1.)))
            : std::make_pair(GmRay(pnt, GmPnt(x + 1., yAt(x + 1.))),
                             GmRay(pnt, GmPnt(x - 1., yAt(x - 1.))));
}

std::vector<GmPnt> GmLine::filterPnts(const std::vector<GmPnt> &pnts) const {
    std::vector<GmPnt> temps;
    const auto checkPoint
            = [&](const GmPnt &pnt) {
                if (includes(pnt)) {
                    temps.push_back(pnt);
                }
            };

    std::for_each(pnts.begin(), pnts.end(), checkPoint);
    
    return temps;
}

GmSeg::GmSeg(const GmPnt &start, const GmPnt &end)
: GmLine{start, end}
{
}

bool GmSeg::includes(const GmPnt &pnt) const {
    const auto a = start(), b = end();
    double dotProduct = (b - a).dot(pnt - a);

    return GmObj::epsEqual((b - a).cross(pnt - a), 0.)
            && GmObj::epsLessEqual(0., dotProduct, 1e-10)
            && GmObj::epsLessEqual(dotProduct, std::pow(length(), 2.), 1e-10);
}

int GmSeg::isect(std::vector<GmPnt> &pnts, const GmObj &obj) const {
    if (obj.type() == GmType::Seg) {
        const auto other = static_cast<const GmSeg *> (&obj);
        const auto p = start(),
                r = end() - start(),
                q = other->start(),
                s = other->end() - other->start();
        double t = (q - p).cross(s) / (r.cross(s)),
                u = (q - p).cross(r) / (r.cross(s));

        if (GmObj::epsEqual(r.cross(s), 0.)) {
            if (GmObj::epsEqual((q - p).cross(r), 0.)) {
                std::cout << "Line segments overlap\n";
                return 0;
            } else {
                // Line segments are parallel and non-intersecting.
                return 0;
            }
        } else {
            if (0. <= t && t <= 1. && 0. <= u && u <= 1.) {
                pnts.push_back(p + t * r);
                return 1;
            } else {
                // Line segments are not parallel but do not intersect.
                return 0;
            }
        }
    } else {
        return obj.isect(pnts, *this);
    }
}

GmPnt GmSeg::midPnt() const {
    return start() + .5 * (end() - start());
}

GmRay::GmRay(const GmPnt &start, const GmPnt &direction)
: GmLine{start, direction}
{
}

bool GmRay::includes(const GmPnt &pnt) const {
    if (!GmLine::includes(pnt)) {
        return false;
    }
    
    double s = 0., e = 0., p = 0.;
    if (GmObj::epsEqual(start().x(), end().x())) {
        s = start().y(); e = end().y(); p = pnt.y();
    } else {
        s = start().x(); e = end().x(); p = pnt.x();
    }

    return ((s > e) && (s > p)) || ((s < e) && (s < p));
}

int GmRay::isect(std::vector<GmPnt> &pnts, const GmObj &obj) const {
    if (obj.type() == GmType::Ray || obj.type() == GmType::Seg) {
        const auto fGetPrim
                = [&]() -> const GmLine * {
                    if (obj.type() == GmType::Ray) {
                        return static_cast<const GmRay*> (&obj);
                    } else {
                        return static_cast<const GmSeg*> (&obj);
                    }
                };
        const GmLine *line = fGetPrim();
        const auto other = GmLine(line->start(), line->end());

        std::vector<GmPnt> temps;
        other.isect(temps, *this);
        pnts = line->filterPnts(temps);

        return pnts.size();
    } else {
        return obj.isect(pnts, *this);
    }

    return -1;
}

GmCirc::GmCirc(const GmPnt &center, double radius, bool ccw)
: m_center{center}
, m_radius{radius}
, m_ccw{ccw}
{
}

GmCirc::GmCirc(const GmPnt &p0, const GmPnt &p1, bool ccw)
: m_center{(p0 + p1) / 2.}
, m_radius{p0.dist(p1) / 2.}
, m_ccw{ccw}
{
}

GmCirc::GmCirc(const GmPnt &p1, const GmPnt &p2, const GmPnt &p3, bool ccw)
: m_ccw{ccw} {
    double ax = p1.x(), ay = p1.y(),
            bx = p2.x(), by = p2.y(),
            cx = p3.x(), cy = p3.y(),
            x1 = (bx + ax) / 2.,
            y11 = (by + ay) / 2.,
            dy1 = bx - ax, dx1 = -(by - ay),
            x2 = (cx + bx) / 2., y2 = (cy + by) / 2.,
            dy2 = cx - bx, dx2 = -(cy - by),
            ox = (y11 * dx1 * dx2 + x2 * dx1 * dy2 - x1 * dy1 * dx2 - y2
                  * dx1 * dx2) / (dx1 * dy2 - dy1 * dx2),
            oy = (ox - x1) * dy1 / dx1 + y11,
            dx = ox - ax, dy = oy - ay;
    
    m_center = GmPnt(ox, oy);
    m_radius = std::sqrt(dx * dx + dy * dy);
}

bool GmCirc::contains(const GmCirc &other) const {
    return m_center.dist(other.center()) + other.radius() < m_radius;
}

bool GmCirc::contains(const GmPnt &pnt) const {
    return m_center.dist(pnt) < m_radius;
}

bool GmCirc::includes(const GmPnt &pnt) const {
    return GmObj::epsEqual(m_center.dist(pnt), m_radius);
}

int GmCirc::isect(std::vector<GmPnt> &pnts, const GmObj &obj) const {
    switch (obj.type()) {
        case GmType::Seg:
        {
            const auto seg = static_cast<const GmSeg *> (&obj);
            
            std::vector<GmPnt> temp;
            doIsect(temp, seg->toGmLine());

            const auto checkPoint
                    = [&](const GmPnt & pnt) {
                        if (seg->includes(pnt)) {
                            pnts.push_back(pnt);
                        }
                    };

            std::for_each(temp.begin(), temp.end(), checkPoint);

            if (pnts.size() == 2) {
                const auto p1 = pnts.at(0), p2 = pnts.at(1);

                if (p1.dist(seg->start()) > p2.dist(seg->start())) {
                    std::reverse(pnts.begin(), pnts.end());
                }
            }
            
            return pnts.size();
        }
        case GmType::Ray:
        {
            const auto ray = static_cast<const GmRay *> (&obj);
            std::vector<GmPnt> temp;
            doIsect(temp, GmLine(ray->start(), ray->end()));

            const auto checkPoint
                    = [&](const GmPnt & pnt) {
                        if (ray->includes(pnt)) {
                            pnts.push_back(pnt);
                        }
                    };

            std::for_each(temp.begin(), temp.end(), checkPoint);
            
            return pnts.size();
        }
        case GmType::Line:
        {
            const auto line = static_cast<const GmLine *> (&obj);
            return doIsect(pnts, *line);
        }
        case GmType::Circ:
        {
            const auto circle = static_cast<const GmCirc *> (&obj);
            return doIsect(pnts, *circle);
        }
        case GmType::Arc:
        case GmType::Pnt:
        {
            return obj.isect(pnts, *this);
        }
        case GmType::Chain:
        case GmType::None:
        default:
        {
            break;
        }
    }

    return -1;
}

Side GmCirc::side(const GmPnt &pnt) const {
    double d = m_center.dist(pnt);
    if (!GmObj::epsEqual(d, m_radius)) {
        return (m_ccw ? d < m_radius : d > m_radius) ? Side::Left : Side::Right;
    }
    
    return Side::None;
}

// FIRST contains the left and SECOND the right arc relative to the 
// straight-line segment defined by START and END.

std::pair<GmArc, GmArc> GmCirc::split(const GmPnt &start, const GmPnt &end) {
    return std::make_pair(GmArc(m_center, end, start), GmArc(m_center, start, end));
}

int GmCirc::doIsect(std::vector<GmPnt> &pnts, const GmLine &line) const {
    double x1 = line.start().x() - m_center.x(),
            y1 = line.start().y() - m_center.y(),
            x2 = line.end().x() - m_center.x(),
            y2 = line.end().y() - m_center.y(),
            dx = x2 - x1, dy = y2 - y1,
            dr = std::sqrt(std::pow(dx, 2.) + std::pow(dy, 2.)),
            D = x1 * y2 - x2 * y1;
    const auto fSign
            = [](double x) {
                return x < 0. ? -1 : 1;
            };
    double delta = std::pow(m_radius, 2.) * std::pow(dr, 2.) - std::pow(D, 2.);

    if (GmObj::epsEqual(delta, 0.)) {
        pnts.push_back(GmPnt(D * dy / std::pow(dr, 2.) + m_center.x(),
            (-D * dx) / std::pow(dr, 2.) + m_center.y()));
    } else if (delta > 0.) {
        pnts.push_back(GmPnt((D * dy - fSign(dy) * dx * std::sqrt(delta)) / std::pow(dr, 2.) + m_center.x(),
            (-D * dx - std::fabs(dy) * std::sqrt(delta)) / std::pow(dr, 2.) + m_center.y()));
        pnts.push_back(GmPnt((D * dy + fSign(dy) * dx * sqrt(delta)) / std::pow(dr, 2.) + m_center.x(),
            (-D * dx + std::fabs(dy) * sqrt(delta)) / std::pow(dr, 2.) + m_center.y()));
    }
    
    return pnts.size();
}

int GmCirc::doIsect(std::vector<GmPnt> &pnts, const GmCirc &other) const {
    // TODO: Handle case when circles (or arcs) overlap!
    assert(m_center != other.center() || !GmObj::epsEqual(m_radius, other.radius()));

    const auto larger = m_radius > other.radius() ? * this : other,
            smaller = m_radius > other.radius() ? other : * this;

    const auto vec = larger.center() - smaller.center();
    double angle = vec.angle(GmPnt(1., 0.));

    if (vec.y() > 0.) {
        angle = M_PI - angle;
    }

    const auto matr = GmRot(angle) * GmTransl{-1 * larger.center()};
    const auto c2 = smaller.center() * matr;
    double d = c2.x(), r = smaller.radius(), R = larger.radius();

    if (GmObj::epsEqual(r + R, d)) {
        pnts.push_back(GmPnt(R, 0.) * matr.inv());
        return 1;
    }
    
    if (GmObj::epsEqual(d, 0.)) {
        return 0;
    }

    double x = (std::pow(d, 2.) - std::pow(r, 2.) + std::pow(R, 2.)) / (2. * d),
            a = (1. / d) * std::sqrt((-d + r - R) * (-d - r + R) * (-d + r + R) * (d + r + R));

    if (!std::isnan(a)) {
        pnts.push_back(GmPnt(x, a / 2.) * matr.inv());

        if (!GmObj::epsEqual(a, 0.)) {
            pnts.push_back(GmPnt(x, -(a / 2.)) * matr.inv());
        }
    }

    return pnts.size();
}

GmArc::GmArc(const GmPnt &center, double radius, double from, double to)
: GmCirc{center, radius}
, m_start{center + radius * GmPnt{std::cos(from), std::sin(from)}}
, m_end{center + radius * GmPnt{std::cos(to), std::sin(to)}}
, m_from{from}
, m_to{to}
{
    assert(from < 2 * M_PI && to < 2 * M_PI);
}

GmArc::GmArc(const GmPnt &center, const GmPnt &start, const GmPnt &end)
: GmCirc{center, center.dist(start)}
, m_start{start}
, m_end{end}
, m_from{std::atan2((start - center).y(), (start - center).x())} // TODO: Check if correct.
, m_to{std::atan2((end - center).y(), (end - center).x())}
{
    const auto fCheckAngle
            = [](double angle) {
                return angle < 0. ? 2. * M_PI + angle : angle;
            };

    m_from = fCheckAngle(m_from);
    m_to = fCheckAngle(m_to);

    assert(m_from < 2. * M_PI && m_to < 2. * M_PI);
    //assert(GmObj::epsEqual(center.dist(start), center.dist(end)));
}

bool GmArc::includes(const GmPnt &pnt) const {
    if (!GmObj::epsEqual(pnt.dist(center()), radius())) {
        return false;
    }

    const auto matr = GmTransl(-1 * center());
    const auto p = pnt*matr;
    double angle = std::atan2(p.y(), p.x());

    if (angle < 0.) {
        angle = angle + 2. * M_PI;
    }

    if (GmObj::epsEqual(angle, m_from) || GmObj::epsEqual(angle, m_to)) {
        return true;
    }

    return (m_from < m_to) ? (m_from < angle && angle < m_to) : (m_from < angle || angle < m_to);
}

int GmArc::isect(std::vector<GmPnt> &pnts, const GmObj &obj) const {
    std::vector<GmPnt> temp;

    if (obj.type() == GmType::Arc) { // Avoid infinite recursion! TODO: Simply this!
        const auto arc = static_cast<const GmArc *> (&obj);
        GmCirc::isect(temp, GmCirc(arc->center(), arc->radius()));

        for (const auto &pnt : temp) {
            if (includes(pnt) && arc->includes(pnt)) {
                pnts.push_back(pnt);
            }
        }
    } else {
        GmCirc::isect(temp, obj);

        for (const auto &pnt : temp) {
            if (includes(pnt)) {
                pnts.push_back(pnt);
            }
        }
    }

    return pnts.size();
}

GmChain::GmChain(const std::vector<GmPnt> &pnts)
: m_pnts{pnts}
{
}

int GmChain::isect(std::vector<GmPnt> &, const GmObj &) const {
    return -1;
}

GmMatrix::GmMatrix(const Matr4 &matrix)
: m_matrix{matrix}
{
}

GmMatrix GmMatrix::operator*(const GmMatrix &other) const {
    GmMatrix c;

    for (int i = 0; i < MATR_DIM; i++)
        for (int j = 0; j < MATR_DIM; j++) {
            double value = 0.;
            for (int k = 0; k < MATR_DIM; k++) {
                value += at(i, k) * other.at(k, j);
            }
            c.set(i, j, value);
        }

    return c;
}

double GmMatrix::at(int i, int j) const {
    assert(i < MATR_DIM && j < MATR_DIM);
    return m_matrix[i][j];
}

void GmMatrix::set(int i, int j, double value) {
    assert(i < MATR_DIM && j < MATR_DIM);
    m_matrix[i][j] = value;
}

// A^-1 = 1/det(A) * adj(a)
GmMatrix GmMatrix::inv() const {
    Matr4 inv;

    inv[0][0] = m_matrix[1][1] * m_matrix[2][2] * m_matrix[3][3] -
            m_matrix[1][1] * m_matrix[2][3] * m_matrix[3][2] -
            m_matrix[2][1] * m_matrix[1][2] * m_matrix[3][3] +
            m_matrix[2][1] * m_matrix[1][3] * m_matrix[3][2] +
            m_matrix[3][1] * m_matrix[1][2] * m_matrix[2][3] -
            m_matrix[3][1] * m_matrix[1][3] * m_matrix[2][2];

    inv[1][0] = -m_matrix[1][0] * m_matrix[2][2] * m_matrix[3][3] +
            m_matrix[1][0] * m_matrix[2][3] * m_matrix[3][2] +
            m_matrix[2][0] * m_matrix[1][2] * m_matrix[3][3] -
            m_matrix[2][0] * m_matrix[1][3] * m_matrix[3][2] -
            m_matrix[3][0] * m_matrix[1][2] * m_matrix[2][3] +
            m_matrix[3][0] * m_matrix[1][3] * m_matrix[2][2];

    inv[2][0] = m_matrix[1][0] * m_matrix[2][1] * m_matrix[3][3] -
            m_matrix[1][0] * m_matrix[2][3] * m_matrix[3][1] -
            m_matrix[2][0] * m_matrix[1][1] * m_matrix[3][3] +
            m_matrix[2][0] * m_matrix[1][3] * m_matrix[3][1] +
            m_matrix[3][0] * m_matrix[1][1] * m_matrix[2][3] -
            m_matrix[3][0] * m_matrix[1][3] * m_matrix[2][1];

    inv[3][0] = -m_matrix[1][0] * m_matrix[2][1] * m_matrix[3][2] +
            m_matrix[1][0] * m_matrix[2][2] * m_matrix[3][1] +
            m_matrix[2][0] * m_matrix[1][1] * m_matrix[3][2] -
            m_matrix[2][0] * m_matrix[1][2] * m_matrix[3][1] -
            m_matrix[3][0] * m_matrix[1][1] * m_matrix[2][2] +
            m_matrix[3][0] * m_matrix[1][2] * m_matrix[2][1];

    inv[0][1] = -m_matrix[0][1] * m_matrix[2][2] * m_matrix[3][3] +
            m_matrix[0][1] * m_matrix[2][3] * m_matrix[3][2] +
            m_matrix[2][1] * m_matrix[0][2] * m_matrix[3][3] -
            m_matrix[2][1] * m_matrix[0][3] * m_matrix[3][2] -
            m_matrix[3][1] * m_matrix[0][2] * m_matrix[2][3] +
            m_matrix[3][1] * m_matrix[0][3] * m_matrix[2][2];

    inv[1][1] = m_matrix[0][0] * m_matrix[2][2] * m_matrix[3][3] -
            m_matrix[0][0] * m_matrix[2][3] * m_matrix[3][2] -
            m_matrix[2][0] * m_matrix[0][2] * m_matrix[3][3] +
            m_matrix[2][0] * m_matrix[0][3] * m_matrix[3][2] +
            m_matrix[3][0] * m_matrix[0][2] * m_matrix[2][3] -
            m_matrix[3][0] * m_matrix[0][3] * m_matrix[2][2];

    inv[2][1] = -m_matrix[0][0] * m_matrix[2][1] * m_matrix[3][3] +
            m_matrix[0][0] * m_matrix[2][3] * m_matrix[3][1] +
            m_matrix[2][0] * m_matrix[0][1] * m_matrix[3][3] -
            m_matrix[2][0] * m_matrix[0][3] * m_matrix[3][1] -
            m_matrix[3][0] * m_matrix[0][1] * m_matrix[2][3] +
            m_matrix[3][0] * m_matrix[0][3] * m_matrix[2][1];

    inv[3][1] = m_matrix[0][0] * m_matrix[2][1] * m_matrix[3][2] -
            m_matrix[0][0] * m_matrix[2][2] * m_matrix[3][1] -
            m_matrix[2][0] * m_matrix[0][1] * m_matrix[3][2] +
            m_matrix[2][0] * m_matrix[0][2] * m_matrix[3][1] +
            m_matrix[3][0] * m_matrix[0][1] * m_matrix[2][2] -
            m_matrix[3][0] * m_matrix[0][2] * m_matrix[2][1];

    inv[0][2] = m_matrix[0][1] * m_matrix[1][2] * m_matrix[3][3] -
            m_matrix[0][1] * m_matrix[1][3] * m_matrix[3][2] -
            m_matrix[1][1] * m_matrix[0][2] * m_matrix[3][3] +
            m_matrix[1][1] * m_matrix[0][3] * m_matrix[3][2] +
            m_matrix[3][1] * m_matrix[0][2] * m_matrix[1][3] -
            m_matrix[3][1] * m_matrix[0][3] * m_matrix[1][2];

    inv[1][2] = -m_matrix[0][0] * m_matrix[1][2] * m_matrix[3][3] +
            m_matrix[0][0] * m_matrix[1][3] * m_matrix[3][2] +
            m_matrix[1][0] * m_matrix[0][2] * m_matrix[3][3] -
            m_matrix[1][0] * m_matrix[0][3] * m_matrix[3][2] -
            m_matrix[3][0] * m_matrix[0][2] * m_matrix[1][3] +
            m_matrix[3][0] * m_matrix[0][3] * m_matrix[1][2];

    inv[2][2] = m_matrix[0][0] * m_matrix[1][1] * m_matrix[3][3] -
            m_matrix[0][0] * m_matrix[1][3] * m_matrix[3][1] -
            m_matrix[1][0] * m_matrix[0][1] * m_matrix[3][3] +
            m_matrix[1][0] * m_matrix[0][3] * m_matrix[3][1] +
            m_matrix[3][0] * m_matrix[0][1] * m_matrix[1][3] -
            m_matrix[3][0] * m_matrix[0][3] * m_matrix[1][1];

    inv[3][2] = -m_matrix[0][0] * m_matrix[1][1] * m_matrix[3][2] +
            m_matrix[0][0] * m_matrix[1][2] * m_matrix[3][1] +
            m_matrix[1][0] * m_matrix[0][1] * m_matrix[3][2] -
            m_matrix[1][0] * m_matrix[0][2] * m_matrix[3][1] -
            m_matrix[3][0] * m_matrix[0][1] * m_matrix[1][2] +
            m_matrix[3][0] * m_matrix[0][2] * m_matrix[1][1];

    inv[0][3] = -m_matrix[0][1] * m_matrix[1][2] * m_matrix[2][3] +
            m_matrix[0][1] * m_matrix[1][3] * m_matrix[2][2] +
            m_matrix[1][1] * m_matrix[0][2] * m_matrix[2][3] -
            m_matrix[1][1] * m_matrix[0][3] * m_matrix[2][2] -
            m_matrix[2][1] * m_matrix[0][2] * m_matrix[1][3] +
            m_matrix[2][1] * m_matrix[0][3] * m_matrix[1][2];

    inv[1][3] = m_matrix[0][0] * m_matrix[1][2] * m_matrix[2][3] -
            m_matrix[0][0] * m_matrix[1][3] * m_matrix[2][2] -
            m_matrix[1][0] * m_matrix[0][2] * m_matrix[2][3] +
            m_matrix[1][0] * m_matrix[0][3] * m_matrix[2][2] +
            m_matrix[2][0] * m_matrix[0][2] * m_matrix[1][3] -
            m_matrix[2][0] * m_matrix[0][3] * m_matrix[1][2];

    inv[2][3] = -m_matrix[0][0] * m_matrix[1][1] * m_matrix[2][3] +
            m_matrix[0][0] * m_matrix[1][3] * m_matrix[2][1] +
            m_matrix[1][0] * m_matrix[0][1] * m_matrix[2][3] -
            m_matrix[1][0] * m_matrix[0][3] * m_matrix[2][1] -
            m_matrix[2][0] * m_matrix[0][1] * m_matrix[1][3] +
            m_matrix[2][0] * m_matrix[0][3] * m_matrix[1][1];

    inv[3][3] = m_matrix[0][0] * m_matrix[1][1] * m_matrix[2][2] -
            m_matrix[0][0] * m_matrix[1][2] * m_matrix[2][1] -
            m_matrix[1][0] * m_matrix[0][1] * m_matrix[2][2] +
            m_matrix[1][0] * m_matrix[0][2] * m_matrix[2][1] +
            m_matrix[2][0] * m_matrix[0][1] * m_matrix[1][2] -
            m_matrix[2][0] * m_matrix[0][2] * m_matrix[1][1];

    double det = m_matrix[0][0] * inv[0][0] + m_matrix[0][1] * inv[1][0]
            + m_matrix[0][2] * inv[2][0] + m_matrix[0][3] * inv[3][0];

    if (GmObj::epsEqual(det, 0.)) {
        // TODO: Handle case when det is 0.
        assert(false);
        return GmMatrix();
    }

    det = 1. / det;

    for (int i = 0; i < MATR_DIM; i++) {
        for (int j = 0; j < MATR_DIM; j++) {
            inv[i][j] = inv[i][j] * det;
        }
    }

    return GmMatrix(inv);
}

double GmMatrix::det() const {
    return m_matrix[0][3] * m_matrix[1][2] * m_matrix[2][1] * m_matrix[3][0] - m_matrix[0][2] * m_matrix[1][3] * m_matrix[2][1] * m_matrix[3][0] -
            m_matrix[0][3] * m_matrix[1][1] * m_matrix[2][2] * m_matrix[3][0] + m_matrix[0][1] * m_matrix[1][3] * m_matrix[2][2] * m_matrix[3][0] +
            m_matrix[0][2] * m_matrix[1][1] * m_matrix[2][3] * m_matrix[3][0] - m_matrix[0][1] * m_matrix[1][2] * m_matrix[2][3] * m_matrix[3][0] -
            m_matrix[0][3] * m_matrix[1][2] * m_matrix[2][0] * m_matrix[3][1] + m_matrix[0][2] * m_matrix[1][3] * m_matrix[2][0] * m_matrix[3][1] +
            m_matrix[0][3] * m_matrix[1][0] * m_matrix[2][2] * m_matrix[3][1] - m_matrix[0][0] * m_matrix[1][3] * m_matrix[2][2] * m_matrix[3][1] -
            m_matrix[0][2] * m_matrix[1][0] * m_matrix[2][3] * m_matrix[3][1] + m_matrix[0][0] * m_matrix[1][2] * m_matrix[2][3] * m_matrix[3][1] +
            m_matrix[0][3] * m_matrix[1][1] * m_matrix[2][0] * m_matrix[3][2] - m_matrix[0][1] * m_matrix[1][3] * m_matrix[2][0] * m_matrix[3][2] -
            m_matrix[0][3] * m_matrix[1][0] * m_matrix[2][1] * m_matrix[3][2] + m_matrix[0][0] * m_matrix[1][3] * m_matrix[2][1] * m_matrix[3][2] +
            m_matrix[0][1] * m_matrix[1][0] * m_matrix[2][3] * m_matrix[3][2] - m_matrix[0][0] * m_matrix[1][1] * m_matrix[2][3] * m_matrix[3][2] -
            m_matrix[0][2] * m_matrix[1][1] * m_matrix[2][0] * m_matrix[3][3] + m_matrix[0][1] * m_matrix[1][2] * m_matrix[2][0] * m_matrix[3][3] +
            m_matrix[0][2] * m_matrix[1][0] * m_matrix[2][1] * m_matrix[3][3] - m_matrix[0][0] * m_matrix[1][2] * m_matrix[2][1] * m_matrix[3][3] -
            m_matrix[0][1] * m_matrix[1][0] * m_matrix[2][2] * m_matrix[3][3] + m_matrix[0][0] * m_matrix[1][1] * m_matrix[2][2] * m_matrix[3][3];
}

GmTransl::GmTransl(const GmPnt &pnt)
: GmMatrix{Matr4{std::array<double, MATR_DIM>{1., 0., 0., pnt.x()},
                 std::array<double, MATR_DIM>{0., 1., 0., pnt.y()},
                 std::array<double, MATR_DIM>{0., 0., 1., pnt.z()},
                 std::array<double, MATR_DIM>{0., 0., 0., 1.}}} {
}

GmRot::GmRot(double rad, const GmPnt &pnt)
: GmMatrix{Matr4{std::array<double, MATR_DIM>{std::cos(rad), -std::sin(rad), -pnt.x() * std::cos(rad) + pnt.y() * std::sin(rad) + pnt.x(), 0.},
                 std::array<double, MATR_DIM>{std::sin(rad), std::cos(rad), -pnt.x() * std::sin(rad) - pnt.y() * std::cos(rad) + pnt.y(), 0.},
                 std::array<double, MATR_DIM>{0., 0., 1., 0.},
                 std::array<double, MATR_DIM>{0., 0., 0., 1.}}} {
}
                 
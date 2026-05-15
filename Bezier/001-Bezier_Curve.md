# Bezier 곡선
## 1. Bézier 곡선의 정의
- n차 Bézier 곡선:
```math
\mathbf{C}(u)=\sum _{i=0}^nB_{i,n}(u)\, \mathbf{P_{\mathnormal{i}}},\quad 0\leq u\leq 1
```
- 여기서 $\mathbf{P_{\mathnormal{i}}}$ 는 control point,
- $B_{i,n}(u)$ 는 n차 Bernstein 다항식:
```math
B_{i,n}(u)={n \choose i}u^i(1-u)^{n-i}=\frac{n!}{i!(n-i)!}u^i(1-u)^{n-i}
```
## 2. Bernstein 다항식의 성질 (네가 그림에 적어둔 것들 수식화)
### 2.1 비음수성 (nonnegativity)
```math
B_{i,n}(u)\geq 0,\quad \forall i,\  n,\  0\leq u\leq 1
```
### 2.2 partition of unity
```math
\sum _{i=0}^nB_{i,n}(u)=1,\quad \forall 0\leq u\leq 1
```
### 2.3 끝점 값
```math
B_{0,n}(0)=1,\quad B_{i>0,n}(0)=0
```

```math
B_{n,n}(1)=1,\quad B_{i\lt n,n}(1)=0
```

### 2.4 최대값 위치
```math
B_{i,n}(u)\  \mathrm{는}\  u=\frac{i}{n}\  \mathrm{에서\  정확히\  한\  번\  최대값을\  갖는다.}
```
- 정확히 쓰면:
```math
\arg \max _{u\in [0,1]}B_{i,n}(u)=\frac{i}{n}
```
### 2.5 대칭성 (symmetry)
```math
B_{i,n}(u)=B_{n-i,n}(1-u)
```
- 즉, $B_{i,n}(u)$ 는 에 대해 대칭.
### 2.6 재귀 정의 (de Casteljau와 연결되는 형태)
```math
B_{i,n}(u)=(1-u)\, B_{i,n-1}(u)+u\, B_{i-1,n-1}(u)
```
여기서 i<0 또는 i>n이면 $B_{i,n}(u)\equiv 0$ 로 정의.
### 2.7 도함수

- 경계에서:
```math
B_{-1,n-1}(u)\equiv 0,\quad B_{n,n-1}(u)\equiv 0
```
## 3. 저차 예시들 (n=1,2)
### 3.1 n=1
```math
B_{0,1}(u)=1-u,\quad B_{1,1}(u)=u
```
```math
\mathbf{C}(u)=(1-u)\mathbf{P_{\mathnormal{0}}}+u\mathbf{P_{\mathnormal{1}}}
```
### 3.2 n=2
```math
B_{0,2}(u)=(1-u)^2
```
```math
B_{1,2}(u)=2u(1-u)
```
```math
B_{2,2}(u)=u^2
```
---
## Bernstein (또는 Bézier) basis의 도함수 수식
```math
\frac{d}{du} B_{i,n}(u)
= B'_{i,n}(u)
= n\big( B_{i-1,n-1}(u) - B_{i,n-1}(u) \big)
```

```math
B_{-1,n-1}(u) \equiv 0,\quad B_{n,n-1}(u) \equiv 0
```
---
## 소스 코드
```cpp
#pragma once
#include "geom.h"
#include "tarray.hpp"

class BezierCurve {
public:
    int dim;
    int degree;
    TArray<Point4D> ctrl;

    BezierCurve();
    ~BezierCurve() = default;

    BezierCurve(const BezierCurve& other);
    BezierCurve& operator=(const BezierCurve& other);

    BezierCurve(BezierCurve&& other) noexcept;
    BezierCurve& operator=(BezierCurve&& other) noexcept;

    explicit BezierCurve(const TArray<Point4D>& controlPoints);
    explicit BezierCurve(const TArray<Point3D>& points);
    explicit BezierCurve(int degree, const TArray<Point4D>& controlPoints);
    bool CreateEmpty(int degree);
    [[nodiscard]] bool IsRational() const;
    [[nodiscard]] Point3D EvalPoint(Real u) const;
    Point4D EvalPoint4D(double u) const;

    bool Evaluate(double t, int der_count, int v_stride, double* v) const;

    void ElevateDegree(int t);

    bool Create(const TArray<Point4D>& controlPoints);
    bool Create(const std::vector<Point4D> &controlPoints);
    bool Create(const TArray<Point3D> &points);
    bool Create(const std::vector<Point3D> &points);


    bool Split(Real u, BezierCurve& left, BezierCurve& right) const;

    [[nodiscard]] bool IsValid() const;

    bool Trim(Real t0, Real t1, BezierCurve& curve) const;
};


struct CubicBezier {
    Point3D B0, B1, B2, B3;
};

bool on_build_conic_bezier(
    const Point3D& P0,
    const Point3D& P1,
    double w1,
    const Point3D& P2,
    BezierCurve& bez);

bool on_make_bezier_conic_arc(
    const Point3D& P0,
    const Vector3D& T0,
    const Point3D& P2,
    const Vector3D& T2,
    const Point3D& P,
    Point3D& P1_out,
    double& w1_out
);

Point3D on_get_bezier_euclid(
    const BezierCurve& bc,
    int idx
);

double on_calc_bezier_ctrl_polygon_length(
    const BezierCurve& bc
);

double on_bezier_chord_length(
    const BezierCurve& bc
);

bool on_bezier_arc_length_recursive(
    const BezierCurve& bc,
    double& eps,
    double& len);

void on_compute_bezier_degree3_chord_parameters(
    const Point3D Q[4],
    double u[4]);

bool on_fit_cubic_bezier_through4_points(
    const Point3D Q[4],
    Point3D Pw_out[4]);

bool on_bezier_cubic_interpolation(
    const Point3D Q[4],
    BezierCurve& out_bez);


bool on_fit_bezier_through_points(
    const TArray<Point3D>& Q,
    int n,
    std::vector<Point3D>& Pw_out);

bool on_bezier_global_interpolation(
    const TArray<Point3D>& Q,
    BezierCurve& out_bez);

Matrix on_bezier_to_power_matrix(int degree);

Point3D on_eval_bezier3(const std::array<Point3D, 4>& q, double t);

CubicBezier on_hermite_to_bezier(
    const Point3D& P0,
    const Point3D& P1,
    const Vector3D& T0,
    const Vector3D& T1,
    double h);

Vector3D on_bezier_first_derivative_at_start(
        const CubicBezier& s);

Vector3D on_bezier_first_derivative_at_end(
        const CubicBezier& s);

Point3D on_evaluate_cubic_bezier_point(
        const CubicBezier& s,
        double t);

Point3D on_evaluate_cubic_bezier_point_de_casteljau(
        const CubicBezier& s,
        double t);

bool on_bezier_split(
    const std::vector<Point4D>& pw,
    double u,
    std::vector<Point4D>& qw,
    std::vector<Point4D>& rw);

bool on_bezier_rec_len_internal(
    const std::vector<Point4D>& pw,
    double& eps,
    double& len);

bool on_bezier_arc_length_recursive(
    const std::vector<Point4D>& pw,
    double eps,
    double& out_len);

bool on_bezier_arc_length(
    const std::vector<Point4D>& pw,
    double tol,
    double& length);

std::tuple<SimpleArray<Point4D>,SimpleArray<Point4D>>
    on_split_bezier_1d_curve(const SimpleArray<Point4D>& a, Real t);

std::vector<std::vector<Real>>
on_bezier_curve_degree_elevation_matrix(
    int degree,
    int increment);

Point2D on_bezier_curve(
    const std::vector<Point2D>& control_points,
    Real u);

Point2D on_bezier_curve_derivative(
    const std::vector<Point2D>& control_points,
    Real u);

std::pair<std::vector<Point4D>, std::vector<Point4D>>
on_split_bezier_curve_lerp(
    std::vector<Point4D> a,
    Real t);

std::pair<Point3D, Point3D>
on_eval_bezier_curve_d1(
    const std::vector<Point4D>& ctrl,
    Real t);

std::pair<std::vector<Point4D>, std::vector<Point4D>>
on_bezier_curve_decasteljau_split(
    const std::vector<Point4D>& ctrl,
    Real t);

std::vector<Point4D>
on_reduce_bezier_curve_degree(
    const std::vector<Point4D>& ctrl,
    int target);

Real on_bezier_degree_reduction_error(
    const std::vector<Point4D>& original,
    const std::vector<Point4D>& reduced,
    int samples);

std::vector<Point4D>
on_elevate_bezier_curve_once(
    const std::vector<Point4D>& ctrl);

std::vector<Point4D>
on_elevate_bezier_curve(
    const std::vector<Point4D>& ctrl,
    int times);
```
```cpp
#include "headers.h"

BezierCurve::BezierCurve() {
    this->degree = 0;
    this->dim = 3;
}


BezierCurve::BezierCurve(const BezierCurve& other)
    : dim(other.dim),
      degree(other.degree),
      ctrl(other.ctrl)
{
}

BezierCurve& BezierCurve::operator=(const BezierCurve& other)
{
    if (this != &other)
    {
        dim = other.dim;
        degree = other.degree;
        ctrl = other.ctrl;
    }
    return *this;
}

BezierCurve::BezierCurve(BezierCurve&& other) noexcept
    : dim(other.dim),
      degree(other.degree),
      ctrl(std::move(other.ctrl))
{
    other.dim = 3;
    other.degree = 0;
}

BezierCurve& BezierCurve::operator=(BezierCurve&& other) noexcept
{
    if (this != &other)
    {
        dim = other.dim;
        degree = other.degree;
        ctrl = std::move(other.ctrl);

        other.dim = 3;
        other.degree = 0;
    }
    return *this;
}

BezierCurve::BezierCurve(const TArray<Point4D> &controlPoints) {
    if (controlPoints.Count() == 0) return;
    this->ctrl = controlPoints;
    this->dim = 3;
    this->degree = controlPoints.Count() - 1;
}

BezierCurve::BezierCurve(const TArray<Point3D> &points) {
    if (points.Count() == 0) return;
    this->ctrl.SetSize(points.Count());
    for (int i = 0; i < points.Count(); ++i) {
        this->ctrl[i] = Point4D(points[i], 0.0);
    }
    this->dim = 3;
    this->degree = points.Count() - 1;
}

BezierCurve::BezierCurve(const int degree, const TArray<Point4D> &controlPoints) {
    if (controlPoints.Count() != degree + 1) return;
    this->degree = degree;
    this->ctrl = controlPoints;
    this->dim = 3;
}

bool BezierCurve::CreateEmpty(const int degree) {
    if (degree < 1) return false;
    this->degree = degree;
    this->ctrl = TArray<Point4D>(degree + 1);
    this->dim = 3;
    return true;
}

bool BezierCurve::IsRational() const {
    bool isRational = false;
    for (int i = 0; i < this->ctrl.Count(); i++) {
        if (ctrl[i].w != 0.0f) {
            isRational = true;
            break;
        }
    }
    return isRational;
}

/// Compute a point on a Bezier curve function
Point3D BezierCurve::EvalPoint(const double u) const {
    const auto p = this->degree;
    auto pt = Point4D(0.0f, 0.0f, 0.0f, 0.0f);
    for (int i = 0; i<=p; i++) {
        const double b = basis::on_bernstein(p, i, u);
        pt.x += b * this->ctrl[i].x;
        pt.y += b * this->ctrl[i].y;
        pt.z += b * this->ctrl[i].z;
        pt.w += b * this->ctrl[i].w;
    }
    return pt.ToEuclid();
}

/// B_CFNEVN: Compute a point on a Bezier curve function
Point4D BezierCurve::EvalPoint4D(const double u) const {
    const auto p = this->degree;
    auto pt = Point4D(0.0f, 0.0f, 0.0f, 0.0f);
    for (int i = 0; i<=p; i++) {
        const double b = basis::on_bernstein(p, i, u);
        pt.x += b * this->ctrl[i].x;
        pt.y += b * this->ctrl[i].y;
        pt.z += b * this->ctrl[i].z;
        pt.w += b * this->ctrl[i].w;
    }
    return pt;
}

bool BezierCurve::Evaluate(
    const double t,
    const int der_count,
    const int v_stride,
    double *v) const
{
    if (v == nullptr)
        return false;

    if (dim <= 0)
        return false;

    if (degree < 0)
        return false;

    if (ctrl.Count() != degree + 1)
        return false;

    if (v_stride < dim)
        return false;

    const bool is_rat = IsRational();
    const int order = degree + 1;
    constexpr int cv_stride = 4; // Point4D = x,y,z,w

    std::vector cv(order * cv_stride, 0.0);
    for (int i = 0; i < order; ++i)
    {
        cv[i * cv_stride + 0] = ctrl[i].x;
        cv[i * cv_stride + 1] = ctrl[i].y;
        cv[i * cv_stride + 2] = ctrl[i].z;
        cv[i * cv_stride + 3] = ctrl[i].w;
    }

    return basis::on_evaluate_bezier(
        dim,
        is_rat,
        order,
        cv_stride,
        cv.data(),
        0.0, 1.0,
        der_count,
        t,
        v_stride,
        v);
}

/// B_CDEGEL: Elevate the degree of a Bezier curve
void BezierCurve::ElevateDegree(int t) {
    if (t <= 0)
        return;

    const int old_degree = this->degree;
    const int new_degree = old_degree + t;

    const auto mat = basis::on_degree_elevation_matrix(old_degree, t);
    auto new_controlPoints = TArray<Point4D>(new_degree+1);
    new_controlPoints.Fill(Point4D(0.0, 0.0, 0.0, 0.0));

    for (int i = 0; i <= new_degree; i++) {
        const int a = (i - t > 0) ? (i - t) : 0;
        const int b = (i < old_degree) ? i : old_degree;
        Point4D q(0.0, 0.0, 0.0, 0.0);

        for (int j = a; j <= b; ++j){
            const double value = mat.Get(i, j);
            q.x += value * this->ctrl[j].x;
            q.y += value * this->ctrl[j].y;
            q.z += value * this->ctrl[j].z;
            q.w += value * this->ctrl[j].w;
        }
        new_controlPoints[i] = q;
    }
    this->ctrl = new_controlPoints;
    this->degree = new_degree;
}

bool BezierCurve::Create(const TArray<Point4D> &controlPoints) {
    if (controlPoints.Count() == 0) return false;
    this->ctrl = controlPoints;
    this->dim = 3;
    this->degree = controlPoints.Count() - 1;
    return true;
}

bool BezierCurve::Create(const std::vector<Point4D> &controlPoints) {

    if (controlPoints.size() == 0) return false;
    this->ctrl.SetData(controlPoints.data(), controlPoints.size());
    this->dim = 3;
    this->degree = controlPoints.size() - 1;
    return true;
}

bool BezierCurve::Create(const TArray<Point3D> &points) {
    if (points.Count() == 0) return false;
    this->ctrl.SetSize(points.Count());
    for (int i = 0; i < points.Count(); ++i) {
        this->ctrl[i] = Point4D(points[i], 0.0);
    }
    this->dim = 3;
    this->degree = points.Count() - 1;
    return true;
}

bool BezierCurve::Create(const std::vector<Point3D> &points) {
    if (points.size() == 0) return false;
    this->ctrl.SetSize(points.size());
    for (int i = 0; i < points.size(); ++i) {
        this->ctrl[i] = Point4D(points[i], 0.0);
    }
    this->dim = 3;
    this->degree = points.size() - 1;
    return true;
}

bool BezierCurve::IsValid() const {
    if (dim <= 0)
        return false;
    if (ctrl.GetCount() == 0)
        return false;
    if (ctrl.GetCount() != degree + 1)
        return false;
    return true;
}

bool BezierCurve::Split(const Real u, BezierCurve& left, BezierCurve& right) const {
    if (!IsValid()) return false;
    auto p = this->degree;
    auto a = this->ctrl;
    auto left_ctrl = TArray<Point4D>(this->degree + 1);
    auto right_ctrl = TArray<Point4D>(this->degree + 1);

    left_ctrl[0] = a[0];
    right_ctrl[p] = a[p];
    for (int k = 1; k <= this->degree; k++) {
        for (int i = 0; i<= p - k; i++) {
            a[i] = a[i].lerp(a[i + 1], u);
        }
        left_ctrl[k] = a[0];
        right_ctrl[p - k] = a[p - k];
    }
    left.Create(left_ctrl);
    right.Create(right_ctrl);
    return true;
}

bool BezierCurve::Trim(Real t0,  Real t1, BezierCurve& curve) const {
    // Simple normalization and clamping
    if (t0 > t1) {
        on_swap(t0, t1);
    }
    if (t1 <= 0.0 || t0 >= 1.0) {
        return false;
    }
    t0 = maths::on_clamp(t0, 0.0, 1.0);
    t1 = maths::on_clamp(t1, 0.0, 1.0);
    if (fabs(t1 - t0) < 1e-15) {
        return false;
    }

    BezierCurve cl0, cr0;
    Split(t1, cl0, cr0);
    if (t0 <= 0.0) {
        curve = cl0;
        return true;
    }
    const double local = t0 / t1;
    BezierCurve cl1, cr1;
    cl0.Split(local, cl1, cr1);
    curve = cr1;
    return true;
}

bool on_build_conic_bezier(
    const Point3D& P0,
    const Point3D& P1,
    const double w1,
    const Point3D& P2,
    BezierCurve& bez
)
{
    const Point4D C0(P0.x, P0.y, P0.z, 1.0);
    const Point4D C1(P1.x, P1.y, P1.z, w1);
    const Point4D C2(P2.x, P2.y, P2.z, 1.0);

    TArray<Point4D> controlPoints(3);
    controlPoints[0] = C0;
    controlPoints[1] = C1;
    controlPoints[2] = C2;
    if (!bez.Create(controlPoints)) return false;

    return true;
}

// Build conic middle control and weight.
// Returns true on success.
bool on_make_bezier_conic_arc(
    const Point3D& P0,
    const Vector3D& T0,
    const Point3D& P2,
    const Vector3D& T2,
    const Point3D& P,
    Point3D& P1_out,
    double& w1_out
)
{
    // Build a plane frame using P0, P2, P (typical conic data are planar).
    Point3D O;
    Vector3D x_axis, y_axis, z_axis;
    if (!maths::on_make_plane_frame(P0, P2, P, O, x_axis, y_axis, z_axis))
        return false;

    // Project to 2D
    Point2D  p0 = maths::on_project_to_plane(P0, O, x_axis, y_axis);
    Point2D  p2 = maths::on_project_to_plane(P2, O, x_axis, y_axis);
    Point2D  pp = maths::on_project_to_plane(P, O, x_axis, y_axis);
    Vector2D t0 = maths::on_project_to_plane(T0, x_axis, y_axis);
    Vector2D t2 = maths::on_project_to_plane(T2, x_axis, y_axis);

    // 1) Try intersect L0 = (p0,t0) and L2 = (p2,t2)
    double tau0 = 0.0, tau2 = 0.0; Point2D p1_2d;
    bool ints = maths::on_intersect_lines_2d(p0, t0, p2, t2, tau0, tau2, p1_2d);

    if (ints)
    {
        // p1 is the intersection
        // To compute w1:
        //  - Intersect the segment p0-p2 with the line through p1 and p.
        Vector2D seg = (p2 - p0).ToVector();
        Vector2D dir = (pp - p1_2d).ToVector();
        double t_seg = 0.0, tl = 0.0; Point2D M;
        bool ints2 = maths::on_intersect_lines_2d(p0, seg, p1_2d, dir, t_seg, tl, M);
        if (!ints2) return false;

        // Require 0<=t_seg<=1 and denominators not near zero
        const double eps = 1e-15;
        if (t_seg < -1e-12 || t_seg > 1.0 + 1e-12) return false;
        if (std::fabs(1.0 - t_seg) <= eps) return false;

        // u from t_seg
        double a = std::sqrt(t_seg / (1.0 - t_seg));
        double u = a / (1.0 + a);

        // Vectors for dot products (all in 2D frame)
        Vector2D V0 = (pp - p0).ToVector();
        Vector2D V1 = (p1_2d - pp).ToVector();
        Vector2D V2 = (pp - p2).ToVector();

        double alf = V0.x * V1.x + V0.y * V1.y;
        double bet = V1.x * V2.x + V1.y * V2.y;
        double gam = V1.x * V1.x + V1.y * V1.y;

        double A = (1.0 - u) * (1.0 - u);
        double B = u * u;
        double C = 2.0 * u * (1.0 - u);

        double num = A * alf + B * bet;
        double den = C * gam;
        if (std::fabs(den) <= eps) return false;

        double w1 = num / den;
        w1_out = w1;

        // Lift p1 back to 3D
        Point3D P1 = O + p1_2d.x * x_axis + p1_2d.y * y_axis;
        P1_out = P1;
        return true;
    }
    // 2) Parallel tangents case (parabola): w1 = 0, P1 is a direction.
    // Intersect L0=(P,T0) with segment S=(P0->P2) to get parameter along T0.
    Point2D  A = pp;
    Vector2D U = t0;
    Point2D  B = p0;
    Vector2D V = (p2 - p0).ToVector();

    double tt = 0.0, ts = 0.0; Point2D X;
    bool ints3 = maths::on_intersect_lines_2d(A, U, B, V, tt, ts, X);
    if (!ints3) return false;

    const double eps = 1e-15;
    if (std::fabs(1.0 - ts) <= eps) return false;
    if (ts < -1e-12 || ts > 1.0 + 1e-12) return false;

    double a = std::sqrt(ts / (1.0 - ts));
    double u = a / (1.0 + a);
    double b = 2.0 * u * (1.0 - u);

    double num = -tt * (1.0 - b);
    if (std::fabs(b) <= eps) return false;

    double scale = num / b;
    w1_out = 0.0;

    // Return a 3D vector encoded in P1_out (no origin) along T0.
    Vector3D t0_u = T0;
    if (!t0_u.IsZero()) {
        // Do not normalize; preserve original scale to match param units
        Vector3D V3 = scale * t0_u;
        P1_out = Point3D(V3.x, V3.y, V3.z);
    }
    else {
        P1_out = Point3D::Origin;
    }
    return true;
}

Point3D on_get_bezier_euclid(
    const BezierCurve& bc,
    const int idx
)
{
    const Point4D cv = bc.ctrl[idx];
    return cv.ToEuclid();
}


double on_calc_bezier_ctrl_polygon_length(
    const BezierCurve& bc
)
{
    const int n = bc.ctrl.Count();
    if (n < 2) return 0.0;
    double len = 0.0;
    Point3D prev = on_get_bezier_euclid(bc, 0);
    for (int i = 1; i < n; ++i)
    {
        Point3D cur = on_get_bezier_euclid(bc, i);
        len += prev.DistanceTo(cur);
        prev = cur;
    }
    return len;
}

double on_bezier_chord_length(
  const BezierCurve& bc
)
{
    const int n = bc.ctrl.Count();
    if (n < 2) return 0.0;
    const Point3D a = on_get_bezier_euclid(bc, 0);
    const Point3D b = on_get_bezier_euclid(bc, n - 1);
    return a.DistanceTo(b);
}

/// Compute arc length of a Bezier curve using recursive subdivision.
///
/// Inputs:
/// - pw: control points Pw[0..p]
/// - eps: tolerance
///
/// Output:
/// - length (Real)
///
/// Notes:
/// - Uses standard recursive bound test:
///     UB = control polygon length
///     LB = straight-line distance P(0)–P(p)
/// - If UB - LB <= eps, approximate length = (UB + LB)/2
bool on_bezier_arc_length_recursive(
    const std::vector<Point4D>& pw,
    double eps,
    double& out_len)
{
    if (pw.size() < 2)
        return false;

    if (eps <= 0.0)
        return false;

    double eps_mut = eps;
    double total2 = 0.0;

    if (!on_bezier_rec_len_internal(pw, eps_mut, total2))
        return false;

    out_len = 0.5 * total2;
    return true;
}


/// Compute arc length of a Bezier curve using recursive subdivision.
bool on_bezier_arc_length_recursive(
    const BezierCurve& bc,
    double& eps,
    double& len
)
{
    if (bc.dim != 3) return false;

    // Bounds
    const double UB = 0.5 * on_calc_bezier_ctrl_polygon_length(bc);
    const double LB = 0.5 * on_bezier_chord_length(bc);

    if (const double del = std::fabs(UB - LB); del <= eps)
    {
        len += UB + LB;
        eps = del;
        return true;
    }

    // Split and recurse
    BezierCurve left, right;
    if (!bc.Split(0.5, left, right)) return false;

    double epl = eps * 0.5;
    if (!on_bezier_arc_length_recursive(left, epl, len)) return false;

    double epr = eps - epl;
    if (!on_bezier_arc_length_recursive(right, epr, len)) return false;

    eps = epl + epr;
    return true;
}

void on_compute_bezier_degree3_chord_parameters(const Point3D Q[4], double u[4])
{
    u[0] = 0.0;
    u[3] = 1.0;

    double d1 = Q[0].DistanceTo(Q[1]);
    double d2 = Q[1].DistanceTo(Q[2]);
    double d3 = Q[2].DistanceTo(Q[3]);
    double sum = d1 + d2 + d3;

    const double eps = ON_TOL12;
    if (!(sum > eps) || !(d1 / sum >= 0.0) || !(d2 / sum >= 0.0))
    {
        // fallback: uniform
        u[1] = 1.0 / 3.0;
        u[2] = 2.0 / 3.0;
    }
    else {
        u[1] = u[0] + d1 / sum;
        u[2] = u[1] + d2 / sum;
    }
}

bool on_fit_cubic_bezier_through4_points(
    const Point3D Q[4],
    Point3D Pw_out[4])
{
    double u[4];
    on_compute_bezier_degree3_chord_parameters(Q, u);

    // Build 4x4 system A * Pw = Q  (vector-valued; solve per component)
    // Rows: u=0, u=u1, u=u2, u=1
    double A[4][4] = { 0 };
    double B0[4], B1[4], B2[4];
    basis::on_bernstein_degree3(0.0, B0);
    basis::on_bernstein_degree3(u[1], B1);
    basis::on_bernstein_degree3(u[2], B2);
    // A rows
    for (int j = 0; j < 4; j++) { A[0][j] = B0[j]; }
    for (int j = 0; j < 4; j++) { A[1][j] = B1[j]; }
    for (int j = 0; j < 4; j++) { A[2][j] = B2[j]; }
    double B3[4]; basis::on_bernstein_degree3(1.0, B3);
    for (int j = 0; j < 4; j++) { A[3][j] = B3[j]; }

    // Solve for x, y, z
    double Ax[4][4], Ay[4][4], Az[4][4];
    for (int r = 0; r < 4; r++) for (int c = 0; c < 4; c++) {
        Ax[r][c] = A[r][c]; Ay[r][c] = A[r][c]; Az[r][c] = A[r][c];
    }
    double bx[4] = { Q[0].x, Q[1].x, Q[2].x, Q[3].x };
    double by[4] = { Q[0].y, Q[1].y, Q[2].y, Q[3].y };
    double bz[4] = { Q[0].z, Q[1].z, Q[2].z, Q[3].z };

    bool okx = maths::on_solve4x4(Ax, bx);
    bool oky = maths::on_solve4x4(Ay, by);
    bool okz = maths::on_solve4x4(Az, bz);
    if (!(okx && oky && okz)) return false;

    for (int i = 0; i < 4; i++) {
        Pw_out[i].x = bx[i];
        Pw_out[i].y = by[i];
        Pw_out[i].z = bz[i];
    }
    return true;
}

bool on_bezier_cubic_interpolation(
    const Point3D Q[4],
    BezierCurve& out_bez)
{
    TArray<Point3D> P(4);
    if (!on_fit_cubic_bezier_through4_points(Q, P.GetData())) return false;
    if (!out_bez.Create(P)) return false; // dim=3, non-rational, order=4
    return true;
}


// Compute Bezier control points Pw[0..n] interpolating Q[0..n] at chord params.
// Pw_out must have size n+1.
bool on_fit_bezier_through_points(
    const TArray<Point3D>& Q,
    int n,
    std::vector<Point3D>& Pw_out)
{
    if (n < 1) return false;
    Pw_out.assign(n + 1, Point3D::UnsetPoint);

    // Parameters
    std::vector<double> u = fitting::on_chord_length_params(Q);

    // Build (n+1)x(n+1) system: A * Pw = Q (vector-valued; solve per component)
    std::vector A(n + 1, std::vector<double>(n + 1, 0.0));
    std::vector<double> B;

    // first and last rows enforce end interpolation
    A[0][0] = 1.0;
    A[n][n] = 1.0;

    for (int i = 1; i < n; ++i)
    {
        basis::on_all_bernstein_degree(n, u[i], B);
        for (int j = 0; j <= n; ++j) A[i][j] = B[j];
    }

    // Solve for x, y, z components
    auto solve_component = [&](auto getter)->std::vector<double>
    {
        std::vector<std::vector<double>> Ac = A;
        std::vector<double> b(n + 1, 0.0);
        for (int i = 0; i <= n; ++i) b[i] = getter(Q[i]);
        if (!maths::on_solve_linear_system(Ac, b)) return std::vector<double>();
        return b;
    };

    auto getx = [](const Point3D& P) { return P.x; };
    auto gety = [](const Point3D& P) { return P.y; };
    auto getz = [](const Point3D& P) { return P.z; };

    const std::vector<double> cx = solve_component(getx);
    const std::vector<double> cy = solve_component(gety);
    const std::vector<double> cz = solve_component(getz);
    if (cx.empty() || cy.empty() || cz.empty()) return false;

    for (int i = 0; i <= n; ++i)
        Pw_out[i].Set(cx[i], cy[i], cz[i]);

    return true;
}


bool on_bezier_global_interpolation(
  const TArray<Point3D>& Q,
  BezierCurve& out_bez)
{
    const int n = Q.Count() - 1;
    std::vector<Point3D> Pw;
    if (!on_fit_bezier_through_points(Q, n, Pw)) return false;
    if (!out_bez.Create(Pw)) return false; // dim=3, non-rational, order=n+1
    return true;
}

/// Bezier(n) -> Power(n)
/// power = T · bezier
/// T[k][i] = ∑_{j} coeff
///   B_i^n(t) = ∑_{k=i..n} C(n,i) C(n-i, k-i) (-1)^{k-i} t^k
/// ⇒ T[k][i] = C(n,i) C(n-i, k-i) (-1)^{k-i}, if k<i 0
Matrix on_bezier_to_power_matrix(const int degree)
{
    const int n = degree;
    Matrix M(n + 1, n + 1);
    for (int m = 0; m <= n; ++m)
    {
        for (int i = 0; i <= n; ++i)
        {
            M[m][i] = 0.0;
        }
    }

    for (int m = 0; m <= n; ++m)
    {
        for (int i = 0; i <= n; ++i)
        {
            if (m >= i)
            {
                const int k = m - i;
                const double val = basis::on_binomial_real(n, i)
                    * basis::on_binomial_real(n - i, k) * (k % 2 ? -1.0 : 1.0);
                M[m][i] = val;
            }
        }
    }
    return M;
}


Point3D on_eval_bezier3(
    const std::array<Point3D, 4>& q,
    const double t)
{
    const double s = 1.0 - t;

    const Point3D a = q[0] * s + q[1] * t;
    const Point3D b = q[1] * s + q[2] * t;
    const Point3D c = q[2] * s + q[3] * t;
    const Point3D d = a * s + b * t;
    const Point3D e = b * s + c * t;
    return d * s + e * t;
}

// Hermite → Bezier (G1)
CubicBezier on_hermite_to_bezier(
  const Point3D& P0,
  const Point3D& P1,
  const Vector3D& T0,
  const Vector3D& T1,
  const double h)
{
    CubicBezier S;
    S.B0 = P0;
    S.B3 = P1;
    S.B1 = P0 + T0 * (h / 3.0);
    S.B2 = P1 - T1 * (h / 3.0);
    return S;
}

Vector3D on_bezier_first_derivative_at_start(
    const CubicBezier& s)
{
    return (s.B1 - s.B0).ToVector() * 3.0;
}

Vector3D on_bezier_first_derivative_at_end(
    const CubicBezier& s)
{
    return (s.B3 - s.B2).ToVector() * 3.0;
}


Point3D on_evaluate_cubic_bezier_point(
    const CubicBezier& s,
    const double t)
{
    const double u = 1.0 - t;

    const double b0 = u * u * u;
    const double b1 = 3.0 * u * u * t;
    const double b2 = 3.0 * u * t * t;
    const double b3 = t * t * t;

    Point3D p;
    p.x = b0 * s.B0.x + b1 * s.B1.x + b2 * s.B2.x + b3 * s.B3.x;
    p.y = b0 * s.B0.y + b1 * s.B1.y + b2 * s.B2.y + b3 * s.B3.y;
    p.z = b0 * s.B0.z + b1 * s.B1.z + b2 * s.B2.z + b3 * s.B3.z;
    return p;
}

Point3D on_evaluate_cubic_bezier_point_de_casteljau(
    const CubicBezier& s,
    const double t)
{
    const Point3D q0 = maths::on_lerp(s.B0, s.B1, t);
    const Point3D q1 = maths::on_lerp(s.B1, s.B2, t);
    const Point3D q2 = maths::on_lerp(s.B2, s.B3, t);

    const Point3D r0 = maths::on_lerp(q0, q1, t);
    const Point3D r1 = maths::on_lerp(q1, q2, t);

    return maths::on_lerp(r0, r1, t);
}


/// Split a Bezier curve at parameter `u` using de Casteljau.
///
/// Bezier curve:
///     C(u) = Σ_{i=0..p} B_i^p(u) * Pw[i]
///
/// This routine splits the Bezier segment into:
///     - left:  defined on [0, u]
///     - right: defined on [u, 1]
///
/// Inputs:
/// - `pw`: original control points Pw[0..p] (length = p+1)
/// - `u` : split parameter in [0,1]
///
/// Outputs:
/// - `(qw, rw)`:
///     qw[0..p] = left segment control points
///     rw[0..p] = right segment control points
///
/// Notes:
/// - This is pure Bezier (no knot vector).
/// - Uses de Casteljau triangle; diagonal gives left/right control nets.
bool on_bezier_split(
    const std::vector<Point4D>& pw,
    double u,
    std::vector<Point4D>& qw,
    std::vector<Point4D>& rw)
{
    const int p = static_cast<int>(pw.size()) - 1;
    if (p < 0) return false;
    if (u < 0.0 || u > 1.0) return false;

    qw.clear();
    rw.clear();

    // special case
    if (u == 0.0)
    {
        qw.assign(p + 1, pw.front());
        rw = pw;
        return true;
    }

    if (u == 1.0)
    {
        qw = pw;
        rw.assign(p + 1, pw.back());
        return true;
    }

    std::vector<Point4D> aw = pw;

    qw.reserve(p + 1);
    rw.resize(p + 1);

    qw.push_back(pw[0]);
    rw[p] = pw[p];

    for (int k = 1; k <= p; ++k)
    {
        for (int i = 0; i <= p - k; ++i)
        {
            aw[i] = aw[i].lerp(aw[i + 1], u);
        }

        qw.push_back(aw[0]);
        rw[p - k] = aw[p - k];
    }
    return true;
}

/// Compute the length of a control polygon
void on_comp_polygon_length(const std::vector<Point4D> &pw, const int p, double &ub) {
    for (int i = 0; i < p; ++i)
        ub += maths_extra::on_euclid_dist(pw[i], pw[i + 1]);
}

/// G_DISCPP: Compute the distance between two control points
double on_chord_length(const std::vector<Point4D> &pw, const int p) {
    const double lb = maths_extra::on_euclid_dist(pw[0], pw[p]);
    return lb;
}

/// Internal recursive routine
/// A method that approximates the length by subdividing the Bezier
/// until the difference between the control polygon length
/// and the chord length becomes sufficiently small,
/// rather than a method that calculates the length using exact integration.
bool on_bezier_rec_len_internal(
    const std::vector<Point4D>& pw,
    double& eps,
    double& len)
{
    const int p = static_cast<int>(pw.size()) - 1;

    // upper bound
    // the control polygon length
    double ub = 0.0;
    on_comp_polygon_length(pw, p, ub);

    // lower bound
    // the chord length
    const double lb = on_chord_length(pw, p);

    if (double del = std::abs(ub - lb); del <= eps)
    {
        len += ub + lb;
        eps = del;
        return true;
    }

    std::vector<Point4D> qw, rw;

    if (!on_bezier_split(pw, 0.5, qw, rw))
        return false;

    double eps_left = eps * 0.5;
    if (!on_bezier_rec_len_internal(qw, eps_left, len))
        return false;

    double eps_right = eps - eps_left;
    if (!on_bezier_rec_len_internal(rw, eps_right, len))
        return false;

    eps = eps_left + eps_right;
    return true;
}

/// Compute arc length of a Bezier curve.
///
/// This is a wrapper around the recursive routine:
///     bezier_reclen_internal()
///
/// Inputs:
/// - pw  : control points Pw[0..p]
/// - tol : tolerance
///
/// Output:
/// - Ok(length)
bool on_bezier_arc_length(
    const std::vector<Point4D>& pw,
    const double tol,
    double& length)
{
    if (pw.size() < 2)
        return false;

    double eps = 2.0 * tol;
    double total = 0.0;

    if (!on_bezier_rec_len_internal(pw, eps, total))
        return false;

    length = 0.5 * total;
    return true;
}

std::tuple<SimpleArray<Point4D>,SimpleArray<Point4D>>
    on_split_bezier_1d_curve(const SimpleArray<Point4D>& a, const Real t)
{
    auto p = a.Count() - 1;
    auto left = SimpleArray<Point4D>(p + 1);
    left.SetCount(p+1);
    auto right = SimpleArray<Point4D>(p + 1);
    right.SetCount(p+1);
    left[0] = a[0];
    right[p] = a[p];

    SimpleArray<Point4D> tmp_a = a;
    for (int k = 1; k<=p; k++)
    {
        for (int i = 0; i <=(p - k); i++) {
            tmp_a[i] = tmp_a[i].Lerp(tmp_a[i + 1], t);
        }
        left[k] = a[0];
        right[p - k] = a[p - k];
    }
    return {left, right};

}

std::vector<std::vector<Real>>
on_bezier_curve_degree_elevation_matrix(
    const int degree,
    const int increment)
{
    if (degree < 0 || increment < 0)
        return {};

    const int p = degree;
    const int t = increment;
    const int q = p + t;

    std::vector dm(
        static_cast<size_t>(q + 1),
        std::vector(static_cast<size_t>(p + 1), 0.0));

    const auto bin = basis::on_pascal_triangle_real(q);
    if (bin.empty())
        return {};

    dm[0][0] = 1.0;
    dm[q][p] = 1.0;

    if (q == 0)
        return dm;

    const int r = q / 2;

    // First half.
    for (int i = 1; i <= r; ++i)
    {
        const Real inv = 1.0 / bin[q][i];

        const int j0 = (std::max)(0, i - t);
        const int j1 = (std::min)(p, i);

        for (int j = j0; j <= j1; ++j)
        {
            dm[i][j] =
                inv *
                bin[p][j] *
                bin[t][i - j];
        }
    }
    // Second half by symmetry:
    //
    // E(i, j) = E(q - i, p - j)
    //
    // This mirrors the original optimized implementation.
    for (int i = r + 1; i < q; ++i)
    {
        const int j0 = (std::max)(0, i - t);
        const int j1 = (std::min)(p, i);

        for (int j = j0; j <= j1; ++j)
        {
            dm[i][j] = dm[q - i][p - j];
        }
    }
    return dm;
}

Point2D on_bezier_curve(
    const std::vector<Point2D>& control_points,
    Real u)
{
    if (control_points.empty())
        return Point2D{};

    const std::size_t n = control_points.size() - 1;

    Point2D result;

    for (std::size_t i = 0; i < control_points.size(); ++i)
    {
        const Real b = basis::on_bernstein(n, i, u);

        result.x += control_points[i].x * b;
        result.y += control_points[i].y * b;
    }
    return result;
}

Point2D on_bezier_curve_derivative(
    const std::vector<Point2D>& control_points,
    Real u)
{
    if (control_points.size() < 2)
        return Point2D{};

    const std::size_t n = control_points.size() - 1;

    Point2D result;

    for (std::size_t i = 0; i < control_points.size(); ++i)
    {
        const Real b_der = basis::on_bernstein_der(n, i, u);

        result.x += control_points[i].x * b_der;
        result.y += control_points[i].y * b_der;
    }

    return result;
}

std::pair<std::vector<Point4D>, std::vector<Point4D>>
on_split_bezier_curve_lerp(std::vector<Point4D> a, const Real t)
{
    const int p = static_cast<int>(a.size()) - 1;

    std::vector<Point4D> left(p + 1);
    std::vector<Point4D> right(p + 1);

    left[0] = a[0];
    right[p] = a[p];

    for (int k = 1; k <= p; ++k)
    {
        for (int i = 0; i <= p - k; ++i)
            a[i] = a[i].Lerp(a[i + 1], t);

        left[k] = a[0];
        right[p - k] = a[p - k];
    }

    return {left, right};
}

std::pair<Point3D, Point3D>
on_eval_bezier_curve_d1(
    const std::vector<Point4D>& ctrl,
    Real t)
{
    const int p = std::max(0, static_cast<int>(ctrl.size()) - 1);

    const auto ders = basis::on_all_ber_ders_1d(p, t, 1);

    const auto& B  = ders[0];
    const auto& dB = ders[1];

    Point3D pos{};
    Point3D du{};

    for (size_t i = 0; i < ctrl.size(); ++i)
    {
        const auto& c = ctrl[i];

        pos.x += B[i] * c.x;
        pos.y += B[i] * c.y;
        pos.z += B[i] * c.z;

        du.x += dB[i] * c.x;
        du.y += dB[i] * c.y;
        du.z += dB[i] * c.z;
    }

    return {pos, du};
}

std::pair<std::vector<Point4D>, std::vector<Point4D>>
on_bezier_curve_decasteljau_split(
    const std::vector<Point4D>& ctrl,
    const Real t)
{
    const int n = static_cast<int>(ctrl.size()) - 1;

    std::vector tri(
        n+1, std::vector<Point4D>(n+1));

    for (int i=0;i<=n;++i)
        tri[0][i] = ctrl[i];

    for (int r=1;r<=n;++r)
    {
        for (int i=0;i<=n-r;++i)
        {
            tri[r][i] = tri[r-1][i].Lerp(tri[r-1][i+1], t);
        }
    }

    std::vector<Point4D> left(n+1), right(n+1);

    for (int i=0;i<=n;++i)
    {
        left[i] = tri[i][0];
        right[i] = tri[n-i][i];
    }

    return {left, right};
}

std::vector<Point4D>
on_reduce_bezier_curve_degree(
    const std::vector<Point4D>& ctrl,
    const int target)
{
    int p = static_cast<int>(ctrl.size()) - 1;

    if (target >= p)
        return ctrl;

    std::vector<Point4D> out(target + 1);
    for (int i = 0; i <= target; ++i)
    {
        const Real t = static_cast<double>(i) / target;
        auto [pt, _] = on_eval_bezier_curve_d1(ctrl, t);
        out[i] = {pt.x, pt.y, pt.z, 1.0};
    }
    return out;
}

Real on_bezier_degree_reduction_error(
    const std::vector<Point4D>& original,
    const std::vector<Point4D>& reduced,
    const int samples)
{
    Real err = 0.0;

    for (int i=0;i<=samples;++i)
    {
        const Real t = static_cast<double>(i) / samples;
        auto [p1,_1] = on_eval_bezier_curve_d1(original,t);
        auto [p2,_2] = on_eval_bezier_curve_d1(reduced,t);
        const Real dx = p1.x - p2.x;
        const Real dy = p1.y - p2.y;
        const Real dz = p1.z - p2.z;
        err += dx*dx + dy*dy + dz*dz;
    }
    return std::sqrt(err / (samples+1));
}

std::vector<Point4D>
on_elevate_bezier_curve_once(const std::vector<Point4D>& ctrl)
{
    if (ctrl.empty())
        return {};

    const int n = static_cast<int>(ctrl.size()) - 1;

    std::vector<Point4D> q(n + 2);

    q[0] = ctrl[0];
    q[n + 1] = ctrl[n];

    for (int i = 1; i <= n; ++i)
    {
        const Real alpha = static_cast<Real>(i) / static_cast<Real>(n + 1);
        const Real beta = 1.0 - alpha;

        q[i].x = alpha * ctrl[i - 1].x + beta * ctrl[i].x;
        q[i].y = alpha * ctrl[i - 1].y + beta * ctrl[i].y;
        q[i].z = alpha * ctrl[i - 1].z + beta * ctrl[i].z;
        q[i].w = alpha * ctrl[i - 1].w + beta * ctrl[i].w;
    }

    return q;
}

std::vector<Point4D>
on_elevate_bezier_curve(
    const std::vector<Point4D>& ctrl,
    int times)
{
    if (times <= 0)
        return ctrl;

    std::vector<Point4D> result = ctrl;

    for (int k = 0; k < times; ++k)
        result = on_elevate_bezier_curve_once(result);

    return result;
}
```

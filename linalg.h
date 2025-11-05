#ifndef LINALG_H
#define LINALG_H

#include <cmath>
#include <Eigen/Dense>

struct Plane {
  double A, B, C, D; // A*x + B*y + C*z + D = 0
  double AA, BB, DD; // z = AA*x + BB*y + DD
};

struct Frame {
  Eigen::Matrix<double, 6, 1> frame;
  Eigen::Matrix4d transform;
};

struct Frene {
  Eigen::Vector3d t{Eigen::Vector3d::Zero()};
  Eigen::Vector3d b{Eigen::Vector3d::Zero()};
  Eigen::Vector3d n{Eigen::Vector3d::Zero()};
  Eigen::Vector3d p{Eigen::Vector3d::Zero()};
  Eigen::Matrix4d transf{Eigen::Matrix4d::Identity()};

  Frene(const Eigen::Vector3d& t_,
        const Eigen::Vector3d& b_,
        const Eigen::Vector3d& n_,
        const Eigen::Vector3d& p_)
      : t(t_), b(b_), n(n_), p(p_) {
    updateTransform();
  }

  void updateTransform() {
    transf.setIdentity();
    transf.block<3,1>(0,0) = t;
    transf.block<3,1>(0,1) = b;
    transf.block<3,1>(0,2) = n;
    transf.block<3,1>(0,3) = p;
  }
};

Eigen::Matrix4d trMatrix4x4(const Eigen::Vector3d& delta) {
  Eigen::Matrix4d T = Eigen::Matrix4d::Identity();
  T.block<3,1>(0,3) = delta;
  return T;
}

Eigen::Matrix4d rotMatrix4x4(double angle, char axis) {
  Eigen::Matrix4d R = Eigen::Matrix4d::Identity();

  const double c = std::cos(M_PI * angle / 180.0);
  const double s = std::sin(M_PI * angle / 180.0);

  switch (axis) {
    case 'x':
      R(1,1) =  c; R(1,2) = -s;
      R(2,1) =  s; R(2,2) =  c;
      break;
    case 'y':
      R(0,0) =  c; R(0,2) =  s;
      R(2,0) = -s; R(2,2) =  c;
      break;
    case 'z':
      R(0,0) =  c; R(0,1) = -s;
      R(1,0) =  s; R(1,1) =  c;
      break;
    default:
      throw std::invalid_argument("axis must be 'x', 'y', or 'z'");
  }

  // Clamp near-zeros like the Python version
  constexpr double EPS = 1e-4;
  R = R.unaryExpr([](double v) { return std::abs(v) <= EPS ? 0.0 : v; });

  return R;
}

Plane pointsToPlane(const Eigen::Ref<const Eigen::VectorXd>& x,
                    const Eigen::Ref<const Eigen::VectorXd>& y,
                    const Eigen::Ref<const Eigen::VectorXd>& z) {

  const long long n = x.size();

  // Build normal-equations matrix U and RHS V for [AA, BB, DD]^T
  Eigen::Matrix3d U;
  U << x.squaredNorm(), x.dot(y),        x.sum(),
       x.dot(y),        y.squaredNorm(), y.sum(),
       x.sum(),         y.sum(),         static_cast<double>(n);

  Eigen::Vector3d V;
  V << x.dot(z), y.dot(z), z.sum();

  // Solve U * P = V robustly
  const auto qr = U.colPivHouseholderQr();
  const Eigen::Vector3d P = qr.solve(V);

  Plane plane;
  plane.AA = P[0];
  plane.BB = P[1];
  plane.DD = P[2];

  // Convert to normalized implicit form with C > 0
  const double denom = plane.AA * plane.AA + plane.BB * plane.BB + 1.0;
  plane.C = std::sqrt(1.0 / denom);
  plane.A = -plane.AA * plane.C;
  plane.B = -plane.BB * plane.C;
  plane.D = -plane.DD * plane.C;

  return plane;
}

Eigen::Vector3d poly(double x0, double x1, double x2,
                     double y0, double y1, double y2)
{
  Eigen::Matrix3d A;
  A << x0 * x0, x0, 1.0,
      x1 * x1, x1, 1.0,
      x2 * x2, x2, 1.0;

  Eigen::Vector3d B(y0, y1, y2);
  return A.colPivHouseholderQr().solve(B); // [a, b, c]
}


Frame getBeltFrame(const Eigen::Vector3d& o,
                       const Eigen::Ref<const Eigen::VectorXd>& x,
                       const Eigen::Ref<const Eigen::VectorXd>& y,
                       const Eigen::Ref<const Eigen::VectorXd>& z)
{
  // 1) Fit plane: A x + B y + C z + D = 0, ||(A,B,C)|| = 1, C >= 0
  const Plane pl = pointsToPlane(x, y, z);

  // 2) Axes: n = belt z; t = projected (and flipped) X-axis; b = n × t (right-handed)
  Eigen::Vector3d n(pl.A, pl.B, pl.C);
  if (n.norm() == 0.0) throw std::runtime_error("Degenerate plane normal");
  n.normalize();

  // pick a helper axis not parallel to n
  const Eigen::Vector3d helper = (std::abs(n.x()) < 0.9) ? Eigen::Vector3d::UnitX()
                                                         : Eigen::Vector3d::UnitY();
  Eigen::Vector3d t = helper - helper.dot(n) * n; // project onto plane
  t.normalize();
  t = -t; // match Python sign convention

  Eigen::Vector3d b = n.cross(t);
  b.normalize();
  t = b.cross(n); // re-orthogonalize for a clean right-handed triad

  // 3) Homogeneous transform
  Eigen::Matrix4d T = Eigen::Matrix4d::Identity();
  T.block<3,1>(0,0) = t;
  T.block<3,1>(0,1) = b;
  T.block<3,1>(0,2) = n;
  T.block<3,1>(0,3) = o;

  // 4) ZYX Euler (yaw-pitch-roll), in degrees
  auto clamp = [](double v){ return std::max(-1.0, std::min(1.0, v)); };
  constexpr double RAD2DEG = 180.0 / 3.14159265358979323846;

  const double A = std::atan2(T(1,0), T(0,0)) * RAD2DEG; // yaw
  const double B = std::asin(clamp(-T(2,0))) * RAD2DEG;  // pitch
  const double C = std::atan2(T(2,1), T(2,2)) * RAD2DEG; // roll

  Frame out;
  out.frame << o.x(), o.y(), o.z(), A, B, C;
  out.transform = T;
  return out;
}

Frene getFreneByPoly(const Eigen::Vector3d& p0,
                     const Eigen::Vector3d& u1,
                     const Eigen::Vector3d& u2,
                     const Eigen::Vector3d& v1)
{
  // 1) Fit y = a0*t^2 + a1*t + a2 through (t,y) points with t=x
  const Eigen::Vector3d coef = poly(
      u1.x(), p0.x(), u2.x(),   // x0, x1, x2
      u1.y(), p0.y(), u2.y()    // y0, y1, y2
      );
  const double a0 = coef(0), a1 = coef(1);

  // 2) Tangent along the curve (x=t, y=..., z=const) at t = p0.x()
  const double tparam = p0.x();
  const double dy_dt  = 2.0 * a0 * tparam + a1;
  Eigen::Vector3d tanu(1.0, dy_dt, 0.0);
  double nt = tanu.norm();
  tanu /= nt;
  if (tanu.x() < 0.0) tanu = -tanu;

  // 3) Build a second direction from v1 - p0 (normalize, handle degeneracy)
  Eigen::Vector3d tanv = v1 - p0;
  tanv.normalize();

  // 4) Orthonormal triad: n = tanu × tanv; b = n × tanu
  Eigen::Vector3d n = tanu.cross(tanv);
  n.normalize();

  Eigen::Vector3d b = n.cross(tanu).normalized();

  return Frene(tanu, b, n, p0);
}


#endif // LINALG_H

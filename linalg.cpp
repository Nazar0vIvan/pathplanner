#include "linalg.h"

// ------------ Math ------------
Eigen::Matrix4d trMatrix4x4(const Eigen::Vector3d& delta)
{
  Eigen::Matrix4d T = Eigen::Matrix4d::Identity();
  T.block<3,1>(0,3) = delta;
  return T;
}

Eigen::Matrix4d rotMatrix4x4(double angleDeg, char axis)
{
  Eigen::Matrix4d R = Eigen::Matrix4d::Identity();
  const double ang = angleDeg * M_PI / 180.0;
  const double c = std::cos(ang), s = std::sin(ang);
  switch (axis) {
    case 'x': R(1,1)=c; R(1,2)=-s; R(2,1)= s; R(2,2)=c; break;
    case 'y': R(0,0)=c; R(0,2)= s; R(2,0)=-s; R(2,2)=c; break;
    case 'z': R(0,0)=c; R(0,1)=-s; R(1,0)= s; R(1,1)=c; break;
  }
  constexpr double EPS = 1e-4;
  R = R.unaryExpr([](double v){ return std::abs(v) <= EPS ? 0.0 : v; });
  return R;
}

Plane pointsToPlane(const Eigen::Ref<const Eigen::VectorXd>& x,
                    const Eigen::Ref<const Eigen::VectorXd>& y,
                    const Eigen::Ref<const Eigen::VectorXd>& z)
{
  const double n = static_cast<double>(x.size());
  Eigen::Matrix3d U;
  U << x.squaredNorm(), x.dot(y),        x.sum(),
       x.dot(y),        y.squaredNorm(), y.sum(),
       x.sum(),         y.sum(),         n;
  Eigen::Vector3d V(x.dot(z), y.dot(z), z.sum());

  const auto qr = U.colPivHouseholderQr();
  const Eigen::Vector3d P = qr.solve(V);

  Plane plane;
  plane.AA = P[0]; plane.BB = P[1]; plane.DD = P[2];
  const double denom = plane.AA*plane.AA + plane.BB*plane.BB + 1.0;
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
  A << x0*x0, x0, 1.0,
       x1*x1, x1, 1.0,
       x2*x2, x2, 1.0;
  Eigen::Vector3d B(y0,y1,y2);
  return A.colPivHouseholderQr().solve(B);
}

EulerSolution rot2euler(const Eigen::Matrix3d &R, bool is_deg)
{
  const double B1 = -std::asin(R(2,0));
  const double B2 = M_PI + std::asin(R(2,0));

  const double cB1 = std::cos(B1);
  const double cB2 = std::cos(B2);

  const double C1 = std::atan2(R(2,1) / cB1, R(2,2) / cB1);
  const double C2 = std::atan2(R(2,1) / cB2, R(2,2) / cB2);

  const double A1 = std::atan2(R(1,0) / cB1, R(0,0) / cB1);
  const double A2 = std::atan2(R(1,0) / cB2, R(0,0) / cB2);

  const double k = is_deg ? (180.0 / M_PI) : 1.0;

  return { k*A1, k*A2, k*B1, k*B2, k*C1, k*C2 };
}

Eigen::Matrix3d euler2rot(double A, double B, double C, bool is_deg)
{
  const double k = is_deg ? (180.0 / M_PI) : 1.0;
  Eigen::Matrix3d R =
      (Eigen::AngleAxisd(A*k, Eigen::Vector3d::UnitZ()) *
       Eigen::AngleAxisd(B*k, Eigen::Vector3d::UnitY()) *
       Eigen::AngleAxisd(C*k, Eigen::Vector3d::UnitX())).toRotationMatrix();
  double eps = 1e-4;
  R = R.unaryExpr([&](double v) { return std::abs(v) <= eps ? 0.0 : v; });
  return R;
}

// ------------ Frene ------------
Frene::Frene(const Eigen::Vector3d& t_,
             const Eigen::Vector3d& b_,
             const Eigen::Vector3d& n_,
             const Eigen::Vector3d& p_)
    : t(t_), b(b_), n(n_), p(p_)
{
  transf.setIdentity();
  transf.block<3,1>(0,0) = t;
  transf.block<3,1>(0,1) = b;
  transf.block<3,1>(0,2) = n;
  transf.block<3,1>(0,3) = p;
}

Frene getFreneByPoly(const Eigen::Vector3d& p0,
                     const Eigen::Vector3d& u1,
                     const Eigen::Vector3d& u2,
                     const Eigen::Vector3d& v1)
{
  const Eigen::Vector3d coef = poly(u1.x(), p0.x(), u2.x(), u1.y(), p0.y(), u2.y());
  const double a0 = coef[0], a1 = coef[1];
  const double dy_dt = 2.0*a0*p0.x() + a1;
  Eigen::Vector3d tanu(1.0, dy_dt, 0.0);
  tanu.normalize();
  if (tanu.x() < 0.0) tanu = -tanu;

  Eigen::Vector3d tanv = (v1 - p0).normalized();
  Eigen::Vector3d n = tanu.cross(tanv).normalized();
  Eigen::Vector3d b = n.cross(tanu).normalized();
  return Frene(tanu, b, n, p0);
}

Frene getFreneByCirc(const Eigen::Vector3d &pt0, const Eigen::Vector3d &ptc)
{
  Eigen::Vector3d v = pt0 - ptc;
  Eigen::Vector3d n = v.normalized();    // unit normal (radial)
  Eigen::Vector3d t(-n.y(), n.x(), 0.0); // in-plane tangent
  if (t.x() < 0.0) t = -t;
  Eigen::Vector3d b = n.cross(t);        // binormal

  return Frene(t, b, n, pt0);
}

// ------------ Blade ------------
static MatN3 jsonArrayToMat3(const QJsonArray& arr)
{
  MatN3 M(static_cast<int>(arr.size()), 3);
  for (int i=0; i<arr.size(); ++i) {
    const QJsonArray triple = arr[i].toArray();
    M(i,0) = triple[0].toDouble();
    M(i,1) = triple[1].toDouble();
    M(i,2) = triple[2].toDouble();
  }
  return M;
}

Airfoil loadBladeJson(const QString& filePath)
{
  QFile f(filePath);
  if (!f.open(QIODevice::ReadOnly))
    throw std::runtime_error(("Cannot open file: " + filePath).toStdString());

  QJsonParseError perr{};
  const QJsonDocument doc = QJsonDocument::fromJson(f.readAll(), &perr);
  if (perr.error)
    throw std::runtime_error(("JSON parse error: " + perr.errorString()).toStdString());
  if (!doc.isArray())
    throw std::runtime_error("Top-level JSON must be an array");

  const QJsonArray top = doc.array();
  Airfoil airfoil;
  airfoil.reserve(top.size());

  for (int i=0; i<top.size(); ++i) {
    const QJsonObject o = top[i].toObject();
    Profile p;
    p.cx = jsonArrayToMat3(o.value("cx").toArray());
    p.cv = jsonArrayToMat3(o.value("cv").toArray());
    p.le = jsonArrayToMat3(o.value("le").toArray());
    p.re = jsonArrayToMat3(o.value("re").toArray());
    airfoil.push_back(std::move(p));
  }
  return airfoil;
}


// ------------ Base ------------
Frame getBeltFrame(const Eigen::Vector3d& o,
                   const Eigen::Ref<const Eigen::VectorXd>& x,
                   const Eigen::Ref<const Eigen::VectorXd>& y,
                   const Eigen::Ref<const Eigen::VectorXd>& z)
{
  const Plane pl = pointsToPlane(x,y,z);
  Eigen::Vector3d n(pl.A,pl.B,pl.C);
  n.normalize();

  const Eigen::Vector3d helper = (std::abs(n.x()) < 0.9) ? Eigen::Vector3d::UnitX() : Eigen::Vector3d::UnitY();
  Eigen::Vector3d t = helper - helper.dot(n)*n;
  t.normalize();
  t = -t;

  Eigen::Vector3d b = n.cross(t).normalized();
  t = b.cross(n);

  Eigen::Matrix4d T = Eigen::Matrix4d::Identity();
  T.block<3,1>(0,0)=t;
  T.block<3,1>(0,1)=b;
  T.block<3,1>(0,2)=n;
  T.block<3,1>(0,3)=o;

  const double A = std::atan2(T(1,0), T(0,0)) * 180.0/M_PI;
  const double B = std::asin(-T(2,0)) * 180.0/M_PI;
  const double C = std::atan2(T(2,1), T(2,2)) * 180.0/M_PI;

  return { o.x(), o.y(), o.z(), A, B, C };
}


// ------------ Cylinder ------------
Cylinder Cylinder::fromAxis(const Eigen::Vector3d& c1,
                            const Eigen::Vector3d& c2,
                            const Eigen::Vector3d& o,
                            double R)
{
  // unit vector along cyl axis in WCS
  Eigen::Vector3d y = c1 - c2;
  y /= y.norm();
  // unit vector along cyl radius in WCS
  Eigen::Vector3d z = Eigen::Vector3d::UnitX() - Eigen::Vector3d::UnitX().dot(y) * y;
  z /= z.norm();
  // third unit vector form right system in WCS
  Eigen::Vector3d x = y.cross(z).normalized();
  // Transform matrix from CSCS to WCS
  Eigen::Matrix4d transform; transform.setIdentity();
  transform.block<3,3>(0,0) << x, y, z;
  transform.block<3,1>(0,3) = o;
  // get Frame
  EulerSolution angles = rot2euler(transform.topLeftCorner<3,3>(), true);
  Frame frame;
  frame << o, angles.A1, angles.B1, angles.C1;

  Pose pose = { frame, transform };

  return { R, pose };
}

Pose Cylinder::surfacePose(double deltaY) const
{
  Pose scs; // (surface coordinate system)

  const Eigen::Matrix3d Rw = pose.R();  // CCS rotation in WCS
  const Eigen::Vector3d o  = pose.t();  // CCS origin in WCS
  const Eigen::Vector3d y  = Rw.col(1); // CCS y-axis
  const Eigen::Vector3d z  = Rw.col(2); // CCS z-axis (radial)

  const Eigen::Vector3d os = o + deltaY * y + R * z;  // SCS origin in WCS

  scs.T = pose.T;
  scs.T.block<3,1>(0,3) = os;
  scs.frame << os, pose.frame(3), pose.frame(4), pose.frame(5); // same angles

  return scs;
}




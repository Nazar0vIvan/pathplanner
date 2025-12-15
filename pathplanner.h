#ifndef PATHPLANNER_H
#define PATHPLANNER_H

#include <Eigen/Core>
#include <Eigen/Geometry>
#include <QString>
#include <QFile>
#include <QJsonDocument>
#include <QJsonArray>
#include <QJsonObject>
#include <cmath>

// ------------ Math ------------

using Vec6d = Eigen::Matrix<double, 6, 1>; // row 6x1

struct Plane {
  double A, B, C, D;   // A*x + B*y + C*z + D = 0
  double AA, BB, DD;   // z = AA*x + BB*y + DD
};

struct Pose {
  Vec6d frame;       // [x,y,z,A,B,C] (deg)
  Eigen::Matrix4d T; // homogeneous transform

  Eigen::Vector3d t() const { return T.block<3,1>(0,3); }
  Eigen::Matrix3d R() const { return T.block<3,3>(0,0); }
};

struct Cylinder {
  static Cylinder fromPoints(const Eigen::Vector3d& c1,
                             const Eigen::Vector3d& c2,
                             const Eigen::Vector3d& pc,
                             double R, double L,
                             char axis = 'y');

  static Cylinder fromAxis(const Eigen::Vector3d& u,
                           const Eigen::Vector3d& pc,
                           double R, double L,
                           char axis = 'y');

  Pose surfacePose(char axis1, double val1,
                   char axisRot, double angleDeg,
                   char axis2, double val2) const;

  double R;
  double L;
  Pose pose;
};

struct EulerSolution {
  double A1, A2,
         B1, B2,
         C1, C2;
};

Eigen::Matrix4d trMatrix4x4(const Eigen::Vector3d& delta);
Eigen::Matrix4d rotMatrix4x4(double angleDeg, char axis);
Plane pointsToPlane(const Eigen::Ref<const Eigen::VectorXd>& x,
                    const Eigen::Ref<const Eigen::VectorXd>& y,
                    const Eigen::Ref<const Eigen::VectorXd>& z);
Eigen::Vector3d poly(double x0, double x1, double x2,
                     double y0, double y1, double y2);

EulerSolution rot2euler(const Eigen::Matrix3d& R, bool is_deg = false);
Eigen::Matrix3d euler2rot(double A, double B, double C, bool is_deg = false);
Eigen::Vector3d axisVec(char axis, double value);
Eigen::Vector3d prjPointToLine(const Eigen::Vector3d& l0,
                               const Eigen::Vector3d& v,
                               const Eigen::Vector3d& p);

// project vector vec onto plane ⟂ n (n assumed unit-ish)
Eigen::Vector3d prjToPerpPlane(const Eigen::Vector3d& vec, const Eigen::Vector3d& n);

// ------------ Frene ------------
struct Frene {
  Frene(const Eigen::Vector3d& t_, const Eigen::Vector3d& b_,
        const Eigen::Vector3d& n_, const Eigen::Vector3d& p_);

  Eigen::Vector3d t{Eigen::Vector3d::Zero()};
  Eigen::Vector3d b{Eigen::Vector3d::Zero()};
  Eigen::Vector3d n{Eigen::Vector3d::Zero()};
  Eigen::Vector3d p{Eigen::Vector3d::Zero()};
  Eigen::Matrix4d transf{Eigen::Matrix4d::Identity()};
};

Frene getFreneByPoly(const Eigen::Vector3d& p0,
                     const Eigen::Vector3d& u1,
                     const Eigen::Vector3d& u2,
                     const Eigen::Vector3d& v1);

Frene getFreneByCirc(const Eigen::Vector3d& pt0,
                     const Eigen::Vector3d& ptc);

// ------------ Blade ------------
using MatN3 = Eigen::Matrix<double, Eigen::Dynamic, 3, Eigen::RowMajor>;
struct Profile { MatN3 cx, cv, le, re; };
using Airfoil = std::vector<Profile>;

Vec6d getBeltFrame(const Eigen::Vector3d& o,
                   const Eigen::Ref<const Eigen::VectorXd>& x,
                   const Eigen::Ref<const Eigen::VectorXd>& y,
                   const Eigen::Ref<const Eigen::VectorXd>& z);

Airfoil loadBladeJson(const QString& filePath);

// ------------ Rsi Trajectory ------------
struct MotionParams {
  double v; // max contour velocity
  double a; // contour acceleration
};

namespace rsi {
  QVector<Vec6d> spline(const QVector<Vec6d> &ref_points, const MotionParams& mp, int decimals = 3);
  QVector<Vec6d> lin(const Vec6d& P1, const Vec6d& P2, const MotionParams& mp, int decimals = 3);
};

void writeOffsetsToJson(const QVector<Vec6d>& offsets, const QString& filePath, int decimals = 3);

#endif // PATHPLANNER_H

#include <QCoreApplication>
#include <QFileInfo>
#include <iostream>

#include "pathplanner.h"

int main(int argc, char *argv[])
{
  QCoreApplication a(argc, argv);

  /*
  // BELT LOCATING: B -> 0
  Eigen::Vector3d oT(1009.15, -16.49, 623.81);
  Eigen::VectorXd xT(9), yT(9), zT(9);
  xT << 996.14, 1010.89, 1010.89, 1023.99, 1014.15, 1014.15, 1004.89, 1004.89, 1009.15;
  yT << -16.14,  -29.24,    0.92,  -16.14,  -10.54,  -22.95,  -22.21,  -10.51,  -16.49;
  zT << 625.57,  623.52,  623.48,  622.35,  623.61,  622.86,  624.73,  624.40,  623.81;

  Vec6d beltFrame = getBeltFrame(oT, xT, yT, zT);

  // BLADE LOCATING: B -> F
  Eigen::Matrix4d ABF_T = trMatrix4x4({0.011, 0.047, 153.319});
  Eigen::Matrix4d ABF_O = rotMatrix4x4(-49.0, 'z');
  Eigen::Matrix4d ABF = ABF_T * ABF_O;
  */

  // ROLLER
  const Eigen::Vector3d ur(-0.0237168939, 0.9997013354, -0.0058948179);
  const Eigen::Vector3d Cr(854.512911, -16.511844, 623.196742); // A point on the axis (near the data “middle”)
  const double Rr = 19.991300;
  const double Lr = 20.0;

  Cylinder rl = Cylinder::fromAxis(ur, Cr, Rr, Lr, 'y');
  Pose rl_s = rl.surfacePose('y', -10.0, 'y', -45.0, 'z', rl.R, false);

  // std::cout << rl_s.frame << std::endl;

  // WORKPIECE
  Eigen::Vector3d Pc = { -0.113702, -0.012406, 111.290488 };

  Eigen::Vector3d Pc11 = { -0.151981, -0.002515, 120.0},
                  Pc12 = { -0.153901, -0.003125, 120.0 },
                  Pc21 = { -0.422887,  0.061220, 180.0 },
                  Pc22 = { -0.423638,  0.065223, 180.0};

  Eigen::Vector3d Pc1 = 0.5 * (Pc11 + Pc12),
                  Pc2 = 0.5 * (Pc21 + Pc22);

  double Rs[] = {
    12.991316,
    12.990244,
    12.998138,
    12.999339,
    13.008986,
    13.009134,
    13.019839,
    13.019753
  };

  const double Rc = Eigen::Map<const Eigen::Matrix<double,8,1>>(Rs).mean();
  const double Lc = 74.0;

  Cylinder cyl = Cylinder::fromTwoPoints(Pc1, Pc2, Pc, Rc, Lc, 'z');

  Eigen::Matrix4d AiC = cyl.surfacePose('z', 54.0+15.0, 'z', 0.0, 'y', -Rc, true).T;

  Eigen::Matrix4d AiT;
  AiT <<  -1.0, 0.0, 0.0, 0.0,
           0.0, 0.0, 1.0, 0.0,
           0.0, 1.0, 0.0, 0.0,
           0.0, 0.0, 0.0, 1.0; // ROLLER SURFACE

  Eigen::Matrix4d ACT = AiT * AiC.inverse();

  EulerSolution angles = rot2euler(ACT.topLeftCorner<3,3>(), true);
  const Eigen::Vector3d o = ACT.block<3,1>(0,3);
  Vec6d frame;
  frame << o.x(), o.y(), o.z(), angles.A1, angles.B1, angles.C1;

  std::cout << "P: " << std::endl << frame << std::endl;


  /*
  Eigen::Matrix4d AiT;
  AiT << -1.0, 0.0,  0.0, 0.0,
          0.0, 1.0,  0.0, 0.0,
          0.0, 0.0, -1.0, 0.0,
          0.0, 0.0,  0.0, 1.0; // BELT

  */

  /*
  AiT <<  0.0,-1.0, 0.0, 0.0,
         -1.0, 0.0, 0.0, 0.0,
          0.0, 0.0,-1.0, 0.0,
          0.0, 0.0, 0.0, 1.0; // WHEEL
  */

  /*
  const std::string menu = "Enter command:\n"
                           "1. read - read blade geo from 99.01.25.242.json\n"
                           "2. spline - generate rsi spline and write it into offsets.json\n"
                           "3. exit - exit program\n:";

  const QString jsonPath = QFileInfo(QString::fromUtf8(__FILE__)).absolutePath() + "/99.01.25.242.json";

  std::string cmd;
  std::cout << menu << std::flush;
  while (std::getline(std::cin, cmd)) {
    if (cmd == "read") {
      // reading blade geo
      auto blade = loadBladeJson(jsonPath);
      std::cout << "Loaded profiles: " << blade.size() << "\n";

      // rsi spline
    } else if (cmd == "spline") {
      Vec6d p0 = { 810.0, 0.0, 940.0, 180.0, 0.0, 180.0 };
      Vec6d p1 = { 1010.0, 0.0, 940.0, 180.0, 0.0, 180.0 };
      QVector<Vec6d> ref_points = { p0, p1 };
      QVector<Vec6d> offsets = rsi::spline(ref_points, { 10.0, 2.0 }, 4);
      const QString path = QFileInfo(QString::fromUtf8(__FILE__)).absolutePath() + "/offsets.json";
      writeOffsetsToJson(offsets, path, 4);

    } else if (cmd == "exit") {
      break;
    } else if (!cmd.empty()) {
      std::cout << "Unknown command: " << cmd << "\n";
    }
    std::cout << menu << std::flush;
  }
*/

  return 0;
}

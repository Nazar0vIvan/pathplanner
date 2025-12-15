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

  // std::cout << "R: " << rl.R << std::endl;
  // std::cout << "L: " << rl.L << std::endl;
  // std::cout << "frame: " << std::endl << rl.pose.frame << std::endl;
  // std::cout << "T: " << std::endl << rl.pose.T << std::endl;

  Pose scs = rl.surfacePose('y',0.0, 'y', -45.0, 'z', rl.R+10.0);
  // std::cout << std::endl << "frame: " << std::endl << scs.frame << std::endl;
  // std::cout << "T: " << std::endl << scs.T << std::endl;

  // WORKPIECE
  Eigen::Vector3d Cwp = { -0.113702, -0.012406, 111.290488 };

  Eigen::Vector3d C11 = { -0.151981, -0.002515, 120.0},
                  C12 = { -0.153901, -0.003125, 120.0 },
                  C21 = { -0.422887,  0.061220, 180.0 },
                  C22 = { -0.423638,  0.065223, 180.0};

  Eigen::Vector3d C1wp = 0.5 * (C11 + C12),
                  C2wp = 0.5 * (C21 + C22);

  double Rswp[] = {
    12.991316,
    12.990244,
    12.998138,
    12.999339,
    13.008986,
    13.009134,
    13.019839,
    13.019753
  };

  const double Rwp = Eigen::Map<const Eigen::Matrix<double,8,1>>(Rswp).mean();
  const double Lwp = 74.0;

  Cylinder wp = Cylinder::fromPoints(C1wp, C2wp, Cwp, Rwp, Lwp, 'z');

  Eigen::Matrix4d AiB = wp.surfacePose('z', Lwp/2, 'z', 0.0, 'x', Rwp).T;

  Eigen::Matrix4d AiT;
  AiT <<   0.0, -1.0, 0.0, 0.0,
           0.0,  0.0, 1.0, 0.0,
          -1.0,  0.0, 0.0, 0.0,
           0.0,  0.0, 0.0, 1.0; // ROLLER SURFACE

  Eigen::Matrix4d ABT = AiT * AiB.inverse();

  EulerSolution angles = rot2euler(ABT.topLeftCorner<3,3>(), true);
  const Eigen::Vector3d o = ABT.block<3,1>(0,3);
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

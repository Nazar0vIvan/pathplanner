#include <QCoreApplication>
#include <QFileInfo>
#include <iostream>

#include "pathplanner.h"

int main(int argc, char *argv[])
{
  QCoreApplication a(argc, argv);

  // BELT LOCATING: B -> 0
  Eigen::Vector3d oT(1009.15, -16.49, 623.81);
  Eigen::VectorXd xT(9), yT(9), zT(9);
  xT << 996.14, 1010.89, 1010.89, 1023.99, 1014.15, 1014.15, 1004.89, 1004.89, 1009.15;
  yT << -16.14,  -29.24,    0.92,  -16.14,  -10.54,  -22.95,  -22.21,  -10.51,  -16.49;
  zT << 625.57,  623.52,  623.48,  622.35,  623.61,  622.86,  624.73,  624.40,  623.81;

  Frame beltFrame = getBeltFrame(oT, xT, yT, zT);

  // BLADE LOCATING: B -> F
  Eigen::Matrix4d ABF_T = trMatrix4x4({0.011, 0.047, 153.319});
  Eigen::Matrix4d ABF_O = rotMatrix4x4(-49.0, 'z');
  Eigen::Matrix4d ABF = ABF_T * ABF_O;

  // SAMPLE LOCATING: S -> F
  // C = [0.012406, 111.290488, 0.113702]
  Eigen::Vector3d C11 = {  0.002515, 120.0, 0.151981 },
                  C12 = {  0.003125, 120.0, 0.153901 },
                  C21 = { -0.061220, 180.0, 0.422887 },
                  C22 = { -0.065223, 180.0, 0.423638 };

  Eigen::Vector3d C1 = 0.5 * (C11 + C12),
                  C2 = 0.5 * (C21 + C22);

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

  double R = Eigen::Map<const Eigen::Matrix<double,8,1>>(Rs).mean();

  Eigen::Matrix4d AiT;
  AiT << -1.0, 0.0,  0.0, 0.0,
          0.0, 1.0,  0.0, 0.0,
          0.0, 0.0, -1.0, 0.0,
          0.0, 0.0,  0.0, 1.0; // BELT
  /*
  AiT <<  0.0,-1.0, 0.0, 0.0,
         -1.0, 0.0, 0.0, 0.0,
          0.0, 0.0,-1.0, 0.0,
          0.0, 0.0, 0.0, 1.0; // WHEEL

  */

  const QString jsonPath = QFileInfo(QString::fromUtf8(__FILE__)).absolutePath() + "/99.01.25.242.json";

  std::string cmd;
  std::cout << "Enter command (read/exit): " << std::flush;
  while (std::getline(std::cin, cmd)) {
    if (cmd == "read") {
      // READING BLADE GEO
      auto blade = loadBladeJson(jsonPath);
      std::cout << "Loaded profiles: " << blade.size() << "\n";


    } else if (cmd == "exit") {
      break;
    } else if (!cmd.empty()) {
      std::cout << "Unknown command: " << cmd << "\n";
    }
    std::cout << "Enter command (read/exit): " << std::flush;
  }

  return 0;
}

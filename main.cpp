#include <QCoreApplication>
#include <Eigen/Dense>
#include <iostream>

#include "linalg.h"

int main(int argc, char *argv[])
{
  QCoreApplication a(argc, argv);

  Eigen::Vector3d oT(1009.15, -16.49, 623.81);
  Eigen::VectorXd xT(9), yT(9), zT(9);
  xT << 996.14, 1010.89, 1010.89, 1023.99, 1014.15, 1014.15, 1004.89, 1004.89, 1009.15;
  yT <<  -16.14,  -29.24,    0.92,  -16.14,  -10.54,  -22.95,  -22.21,  -10.51,  -16.49;
  zT <<  625.57,  623.52,  623.48,  622.35,  623.61,  622.86,  624.73,  624.40,  623.81;

  Frame beltFrame = getBeltFrame(oT, xT, yT, zT);

  Eigen::Matrix4d ABF_T = trMatrix4x4({0.011, 0.047, 153.319});
  Eigen::Matrix4d ABF_O = rotMatrix4x4(-49.0, 'z');
  Eigen::Matrix4d ABF = ABF_T * ABF_O;



  return a.exec();
}

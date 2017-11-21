#include <iostream>
#include "kalman_filter.h"

using Eigen::MatrixXd;
using Eigen::VectorXd;
using namespace std;

// Please note that the Eigen library does not initialize 
// VectorXd or MatrixXd objects with zeros upon creation.

KalmanFilter::KalmanFilter() {}

KalmanFilter::~KalmanFilter() {}

void KalmanFilter::Init(VectorXd &x_in, MatrixXd &P_in, MatrixXd &F_in,
                        MatrixXd &H_in, MatrixXd &R_in, MatrixXd &Q_in) {
  x_ = x_in;
  P_ = P_in;
  F_ = F_in;
  H_ = H_in;
  R_ = R_in;
  Q_ = Q_in;
}

void KalmanFilter::Predict() {
  /**
  TODO:
    * predict the state
  */
  x_ = F_ * x_;
  MatrixXd Ft = F_.transpose();
  P_ = F_ * P_ * Ft + Q_;
}

void KalmanFilter::Update(const VectorXd &z) {
  /**
  TODO:
    * update the state by using Kalman Filter equations
  */
  VectorXd y = z - H_ * x_;

  UpdateKF(y);
}

void KalmanFilter::UpdateEKF(const VectorXd &z) {
  /**
  TODO:
    * update the state by using Extended Kalman Filter equations
  */
  float rho = sqrt(x_(0) * x_(0) + x_(1) * x_(1));
  if (rho < 0.0000001f)
    rho = 0.0000001f;
  float phi = atan2(x_(1),x_(0));
  /*if (phi < -M_PI) {
    while (phi < -M_PI)
      phi += 2*M_PI;
  } else if (phi > M_PI) {
    while (phi > M_PI)
      phi -= 2*M_PI;
  }*/
  float rhodot = (x_(0) * x_(2) + x_(1) * x_(3)) / rho;
  VectorXd hx = VectorXd(3);
  hx << rho, phi, rhodot;
  VectorXd y = z - hx;
  
  float theta = y(1);
  if (theta < -M_PI) {
    while (theta < -M_PI)
      theta += 2*M_PI;
  } else if (theta > M_PI) {
    while (theta > M_PI)
      theta -= 2*M_PI;
  }
  y(1) = theta;
  UpdateKF(y);
}

void KalmanFilter::UpdateKF(const VectorXd &y) {
	MatrixXd S = H_ * P_ * H_.transpose() + R_;
	MatrixXd K = P_ * H_.transpose() * S.inverse();

	MatrixXd I = MatrixXd::Identity(x_.size(), x_.size() );

  x_ = x_ + K * y;
  P_ = (I - K * H_) * P_;
}

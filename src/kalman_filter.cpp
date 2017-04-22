#include "kalman_filter.h"

using Eigen::MatrixXd;
using Eigen::VectorXd;

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
    * predict the state (L5:13)
	For this project, however, we do not need to use F_​j. 
	If we had been using a non-linear model in the prediction step, we would need to replace the F matrix with its Jacobian, F_j.
	However, we are using a linear model for the prediction step. 
	So, for the prediction step, we can still use the regular Kalman filter equations and the F matrix rather than the extended Kalman filter equations.
  */
	x_ = F_ * x_;
	MatrixXd Ft = F_.transpose();
	P_ = F_ * P_ * Ft + Q_;
}

void KalmanFilter::Update(const VectorXd &z) {
  /**
  TODO:
    * update the state by using Kalman Filter equations (L5:13)
  */
	VectorXd z_pred = H_ * x_;
	VectorXd y = z - z_pred;
	MatrixXd Ht = H_.transpose();
	MatrixXd S = H_ * P_ * Ht + R_;
	MatrixXd Si = S.inverse();
	MatrixXd PHt = P_ * Ht;
	MatrixXd K = PHt * Si;

	//new estimate
	x_ = x_ + (K * y);
	long x_size = x_.size();
	MatrixXd I = MatrixXd::Identity(x_size, x_size);
	P_ = (I - K * H_) * P_;
}

void KalmanFilter::UpdateEKF(const VectorXd &z) {
  /**
  TODO:
    * update the state by using Extended Kalman Filter equations
  */

	// To get z_pred for RADAR
	double px = x_(0);
	double py = x_(1);
	double vx = x_(2);
	double vy = x_(3);

	double rho_pred = sqrt(pow(px, 2) + pow(py, 2));

	double phi_pred = 0.0;
	if (fabs(px) > 0.001) { //float or double absolute value
		phi_pred = atan2(py, px);
	}

	double rhodot_pred = 0.0;
	if (fabs(rho_pred) > 0.001) {
		rhodot_pred = (px*vx + py*vy) / rho_pred;
	}

	VectorXd z_pred(3);
	z_pred << rho_pred, phi_pred, rhodot_pred; // h(x')

	VectorXd y = z - z_pred;
	// take care that y[1] i.e. phi is not >pi and <-pi.
	if (y[1] > M_PI) {
		y[1] = y[1] - 2*M_PI;
	}
	else if (y[1] < -M_PI) {
		y[1] = y[1] + 2 * M_PI;
		}
	else {
		//do nothing keep it as it is
	}
	MatrixXd Hj = H_; // here H_ is Jacobian matrix Hj
	MatrixXd Ht = Hj.transpose();
	MatrixXd S = Hj * P_ * Ht + R_;
	MatrixXd Si = S.inverse();
	MatrixXd PHt = P_ * Ht;
	MatrixXd K = PHt * Si;

	//new estimate
	x_ = x_ + (K * y);
	long x_size = x_.size();
	MatrixXd I = MatrixXd::Identity(x_size, x_size);
	P_ = (I - K * Hj) * P_;

}

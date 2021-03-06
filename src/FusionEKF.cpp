#include "FusionEKF.h"
#include "tools.h"
#include "Eigen/Dense"
#include <iostream>

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::vector;

/*
 * Constructor.
 */
FusionEKF::FusionEKF() {
  is_initialized_ = false;

  previous_timestamp_ = 0;

  // initializing matrices
  R_laser_ = MatrixXd(2, 2);
  R_radar_ = MatrixXd(3, 3);
  H_laser_ = MatrixXd(2, 4);
  Hj_ = MatrixXd(3, 4);

  //measurement covariance matrix - laser
  R_laser_ << 0.0225, 0,
        0, 0.0225;

  //measurement covariance matrix - radar
  R_radar_ << 0.09, 0, 0,
        0, 0.0009, 0,
        0, 0, 0.09;

  /**
  TODO:
    * Finish initializing the FusionEKF.
    * Set the process and measurement noises
  */
  // set the process noise components
  H_laser_ << 1, 0, 0, 0,
	  0, 1, 0, 0;

  Hj_ << 1, 0, 0, 0,
	  0, 1, 0, 0,
	  0, 0, 1, 0;
  
  // set the acceleration measurement noise components
  noise_ax = 9;
  noise_ay = 9;

}

/**
* Destructor.
*/
FusionEKF::~FusionEKF() {}

void FusionEKF::ProcessMeasurement(const MeasurementPackage &measurement_pack) {


  /*****************************************************************************
   *  Initialization
   ****************************************************************************/
  if (!is_initialized_) {
    /**
    TODO:
      * Initialize the state ekf_.x_ with the first measurement.
      * Create the covariance matrix.
      * Remember: you'll need to convert radar from polar to cartesian coordinates.
    */
	float rho;
	float phi;
	float rhodot;
	float px;
	float py;
	float vx;
	float vy;

    // first measurement
    //cout << "EKF: " << endl;
    ekf_.x_ = VectorXd(4);
    ekf_.x_ << 1, 1, 1, 1;

	// Initialize transition matrix F
	ekf_.F_ = MatrixXd(4, 4);
	ekf_.F_ << 1, 0, 1, 0,
		0, 1, 0, 1,
		0, 0, 1, 0,
		0, 0, 0, 1;

	// Covariance matrix 
	ekf_.P_ = MatrixXd(4, 4);
	ekf_.P_ << 1, 0, 0, 0,
		0, 1, 0, 0,
		0, 0, 1000, 0,
		0, 0, 0, 1000;

	// Q matrix
	ekf_.Q_ = MatrixXd(4, 4);
	ekf_.Q_ << 1, 1, 1, 1,
		1, 1, 1, 1,
		1, 1, 1, 1,
		1, 1, 1, 1;

    if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
      /**
      Convert radar from polar to cartesian coordinates and initialize state.
      */
		rho = measurement_pack.raw_measurements_(0);
		phi = measurement_pack.raw_measurements_(1);
		rhodot = measurement_pack.raw_measurements_(2);

		px = rho * cos(phi);
		py = rho * sin(phi);
		vx = rhodot * cos(phi);
		vy = rhodot * sin(phi);

		ekf_.x_ << px, py, vx, vy;
    }
    else if (measurement_pack.sensor_type_ == MeasurementPackage::LASER) {
      /**
      Initialize state.
      */
	  ekf_.x_ << measurement_pack.raw_measurements_(0), measurement_pack.raw_measurements_(1), 0.0, 0.0;
    }

    // done initializing, no need to predict or update
    is_initialized_ = true;
	previous_timestamp_ = measurement_pack.timestamp_; // save the current timestamp for next calculation
    return;
  }

  /*****************************************************************************
   *  Prediction
   ****************************************************************************/

  /**
   TODO:
     * Update the state transition matrix F according to the new elapsed time.
      - Time is measured in seconds.
     * Update the process noise covariance matrix.
     * Use noise_ax = 9 and noise_ay = 9 for your Q matrix.
   */
   //elapsed time
  float dt = (measurement_pack.timestamp_ - previous_timestamp_) / 1000000.0;
  previous_timestamp_ = measurement_pack.timestamp_;

  //initialize some variables to get a clean code
  float dt_2 = dt * dt;
  float dt_3 = dt_2 * dt;
  float dt_4 = dt_3 * dt;

  //Modify the F matrix so that the time is integrated
  ekf_.F_(0, 2) = dt;
  ekf_.F_(1, 3) = dt;

  //set the process covariance matrix Q
  ekf_.Q_ << dt_4 / 4 * noise_ax, 0, dt_3 / 2 * noise_ax, 0,
	  0, dt_4 / 4 * noise_ay, 0, dt_3 / 2 * noise_ay,
	  dt_3 / 2 * noise_ax, 0, dt_2*noise_ax, 0,
	  0, dt_3 / 2 * noise_ay, 0, dt_2*noise_ay;

  ekf_.Predict();

  /*****************************************************************************
   *  Update
   ****************************************************************************/

  /**
   TODO:
     * Use the sensor type to perform the update step.
     * Update the state and covariance matrices.
   */
  float rho;
  float phi;
  float rhodot;

  if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
    // Radar updates
	  rho = measurement_pack.raw_measurements_(0);
	  phi = measurement_pack.raw_measurements_(1);
	  rhodot = measurement_pack.raw_measurements_(2);

	  VectorXd z_Radar(3);
	  z_Radar << rho, phi, rhodot;

	  // For calculating JacobianMatrix, create an instance of Tools package
	  Tools tools;
	  // Define the Hj Matrix
	  ekf_.H_ = tools.CalculateJacobian(ekf_.x_);
	  //Get the measurement covariance matrix
	  ekf_.R_ = R_radar_;
	  //Update the sensor measurements
	  ekf_.UpdateEKF(z_Radar);
  } else {
    // Laser updates
	  VectorXd z_Laser(2);
	  z_Laser << measurement_pack.raw_measurements_(0), measurement_pack.raw_measurements_(1);
	  ekf_.H_ = H_laser_; // H matrix for the case of laser
	  ekf_.R_ = R_laser_;// measurement covariance matrix
	  ekf_.Update(z_Laser);
  }

  // print the output
  //cout << "x_ = " << ekf_.x_ << endl;
  //cout << "P_ = " << ekf_.P_ << endl;
}

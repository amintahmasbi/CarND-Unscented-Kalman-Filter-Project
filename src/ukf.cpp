#include "ukf.h"
#include "tools.h"
#include "Eigen/Dense"
#include <iostream>

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::vector;

/**
 * Initializes Unscented Kalman filter
 */
UKF::UKF() {

  // set to false initially, set to true in first call of ProcessMeasurement
  is_initialized_ = false;

  //set state dimension
  n_x_ = 5;

  //set augmented dimension
  n_aug_ = 7;

  // set radar meas. dimensions
  n_z_radar_ = 3;

  // set lidar meas. dimensions
  n_z_laser_ = 2;

  //define spreading parameter
  lambda_ = 3 - n_aug_;

  // Initial time in us
  time_us_ = 0;

  // Initial NIS for radar
  NIS_radar_ = 0.0;

  // Initial NIS for lidar
  NIS_laser_ = 0.0;

  // if this is false, laser measurements will be ignored (except during init)
  use_laser_ = true;

  // if this is false, radar measurements will be ignored (except during init)
  use_radar_ = true;

  // initial state vector
  x_ = VectorXd::Zero(n_x_);

  // initial covariance matrix
  P_ = MatrixXd::Identity(n_x_,n_x_);

  // Process noise standard deviation longitudinal acceleration in m/s^2
  std_a_ = 0.25;

  // Process noise standard deviation yaw acceleration in rad/s^2
  std_yawdd_ = 20.0/100.0*M_PI;

  // Laser measurement noise standard deviation position1 in m
  std_laspx_ = 0.15;

  // Laser measurement noise standard deviation position2 in m
  std_laspy_ = 0.15;

  // Radar measurement noise standard deviation radius in m
  std_radr_ = 0.3;

  // Radar measurement noise standard deviation angle in rad
  std_radphi_ = 0.03;

  // Radar measurement noise standard deviation radius change in m/s
  std_radrd_ = 0.3;

}

UKF::~UKF() {}

/**
 * @param {MeasurementPackage} meas_package The latest measurement data of
 * either radar or laser.
 */
void UKF::ProcessMeasurement(MeasurementPackage meas_package) {


  /**
    * Initialize the state x_ with the first measurement.
    * Create the covariance matrix.
    * Remember: you'll need to convert radar from polar to Cartesian coordinates.
  */

  /*****************************************************************************
   *  Initialization
   ****************************************************************************/
  if (!is_initialized_) {

    if (use_radar_ && meas_package.sensor_type_ == MeasurementPackage::RADAR) {
      /**
      Convert radar from polar to Cartesian coordinates and initialize state.
      */
      float ro = meas_package.raw_measurements_(0);
      float phi = meas_package.raw_measurements_(1);
      // NOTE: ro_dot is not the actual speed (magnitude or direction), it is the speed in the direction of ro (range) vector
      float ro_dot = meas_package.raw_measurements_(2);
      x_ << ro * cos(phi), ro * sin(phi), ro_dot, 0, 0; //estimate the initial speed, too.
      is_initialized_ = true;

    }
    else if (use_laser_ && meas_package.sensor_type_ == MeasurementPackage::LASER) {
      /**
      Initialize state.
      */
      //set the state with the initial location and zero velocity
      x_ << meas_package.raw_measurements_(0), meas_package.raw_measurements_(1), 0, 0, 0;
      is_initialized_ = true;

    }

    time_us_ = meas_package.timestamp_;

    // done initializing, no need to predict or update

    //cout << "UKF Initialization " << endl;

    return;
  }

  /*****************************************************************************
   *  Prediction
   ****************************************************************************/

  /**
     * Update the state transition matrix F according to the new elapsed time.
      - Time is measured in seconds.
   */

  //compute the time elapsed between the current and previous measurements
  const float delta_t = (meas_package.timestamp_ - time_us_) / 1000000.0; //dt - expressed in seconds
  time_us_ = meas_package.timestamp_;


  Prediction(delta_t);

  /*****************************************************************************
   *  Update
   ****************************************************************************/

  /**
     * Use the sensor type to perform the update step.
     * Update the state and covariance matrices.
   */

  if (use_radar_ && meas_package.sensor_type_ == MeasurementPackage::RADAR) {
    // Radar measurement updates
    UpdateRadar(meas_package);

  } else if (use_laser_ && meas_package.sensor_type_ == MeasurementPackage::LASER){
    // Laser measurement updates
    UpdateLidar(meas_package);

  }

  // print the output
//  cout << "x_ = " << x_ << endl;
//  cout << "P_ = " << P_ << endl;
}

/**
 * Predicts sigma points, the state, and the state covariance matrix.
 * @param {double} delta_t the change in time (in seconds) between the last
 * measurement and this one.
 */
void UKF::Prediction(double delta_t) {
  /**
  Estimate the object's location. Modify the state
  vector, x_. Predict sigma points, the state, and the state covariance matrix.
  */

  //create sigma point matrix
  MatrixXd Xsig_aug = MatrixXd(n_aug_, 2 * n_aug_ + 1);
  AugmentedSigmaPoints(&Xsig_aug);
  SigmaPointPrediction(&Xsig_pred_, Xsig_aug, delta_t);
  PredictMeanAndCovariance(&x_, &P_, &weights_);
}

/**
* Creates augmented mean state, remember mean of noise is zero
* Creates augmented covariance matrix
* Creates square root matrix
* Creates and returns augmented sigma points
* @param {MatrixXd} Xsig_out
*/

void UKF::AugmentedSigmaPoints(MatrixXd* Xsig_out) {

  //create augmented mean vector
  VectorXd x_aug = VectorXd(n_aug_);

  //create augmented state covariance
  MatrixXd P_aug = MatrixXd(n_aug_, n_aug_);

  //create sigma point matrix
  MatrixXd Xsig_aug = MatrixXd(n_aug_, 2 * n_aug_ + 1);

  //create augmented mean state
  x_aug.head(5) = x_;
  x_aug(5) = 0;
  x_aug(6) = 0;

  //create augmented covariance matrix
  P_aug.fill(0.0);
  P_aug.topLeftCorner(5,5) = P_;
  P_aug(5,5) = std_a_*std_a_;
  P_aug(6,6) = std_yawdd_*std_yawdd_;

  //create square root matrix
  MatrixXd L = P_aug.llt().matrixL();

  //create augmented sigma points
  Xsig_aug.col(0)  = x_aug;
  for (int i = 0; i< n_aug_; i++)
  {
    Xsig_aug.col(i+1)       = x_aug + sqrt(lambda_+n_aug_) * L.col(i);
    Xsig_aug.col(i+1+n_aug_) = x_aug - sqrt(lambda_+n_aug_) * L.col(i);
  }

  //print result
//  std::cout << "Xsig_aug = " << std::endl << Xsig_aug << std::endl;

  //write result
  *Xsig_out = Xsig_aug;

}

/**
* Predict sigma points while avoiding division by zero
* @param {MatrixXd} Xsig_out
*/

void UKF::SigmaPointPrediction(MatrixXd* Xsig_out,const MatrixXd& Xsig_aug, const double delta_t) {

  //create matrix with predicted sigma points as columns
  MatrixXd Xsig_pred = MatrixXd(n_x_, 2 * n_aug_ + 1);

  //predict sigma points
  for (int i = 0; i< 2*n_aug_+1; i++)
  {
    //extract values for better readability
    double p_x = Xsig_aug(0,i);
    double p_y = Xsig_aug(1,i);
    double v = Xsig_aug(2,i);
    double yaw = Xsig_aug(3,i);
    double yawd = Xsig_aug(4,i);
    double nu_a = Xsig_aug(5,i);
    double nu_yawdd = Xsig_aug(6,i);

    //predicted state values
    double px_p, py_p;

    //avoid division by zero
    if (fabs(yawd) > 0.001) {
        px_p = p_x + v/yawd * ( sin (yaw + yawd*delta_t) - sin(yaw));
        py_p = p_y + v/yawd * ( cos(yaw) - cos(yaw+yawd*delta_t) );
    }
    else {
        px_p = p_x + v*delta_t*cos(yaw);
        py_p = p_y + v*delta_t*sin(yaw);
    }

    double v_p = v;
    double yaw_p = yaw + yawd*delta_t;
    double yawd_p = yawd;

    //add noise
    px_p = px_p + 0.5*nu_a*delta_t*delta_t * cos(yaw);
    py_p = py_p + 0.5*nu_a*delta_t*delta_t * sin(yaw);
    v_p = v_p + nu_a*delta_t;

    yaw_p = yaw_p + 0.5*nu_yawdd*delta_t*delta_t;
    yawd_p = yawd_p + nu_yawdd*delta_t;

    //write predicted sigma point into right column
    Xsig_pred(0,i) = px_p;
    Xsig_pred(1,i) = py_p;
    Xsig_pred(2,i) = v_p;
    Xsig_pred(3,i) = yaw_p;
    Xsig_pred(4,i) = yawd_p;
  }

  //print result
//  std::cout << "Xsig_pred = " << std::endl << Xsig_pred << std::endl;

  //write result
  *Xsig_out = Xsig_pred;

}

/**
* Set weights
* Predict state mean and covariance
* @param {VectorXd} x_out
* @param {MatrixXd} P_out
*/

void UKF::PredictMeanAndCovariance(VectorXd* x_out, MatrixXd* P_out, VectorXd* weights_out) {

  //create vector for weights
  VectorXd weights = VectorXd(2*n_aug_+1);

  //create vector for predicted state
  VectorXd x = VectorXd(n_x_);

  //create covariance matrix for prediction
  MatrixXd P = MatrixXd(n_x_, n_x_);

  // set weights
  double weight_0 = lambda_/(lambda_+n_aug_);
  weights(0) = weight_0;
  for (int i=1; i<2*n_aug_+1; i++) {  //2n+1 weights
    double weight = 0.5/(n_aug_+lambda_);
    weights(i) = weight;
  }

  //predicted state mean
  x.fill(0.0);
  for (int i = 0; i < 2 * n_aug_ + 1; i++) {  //iterate over sigma points
    x = x+ weights(i) * Xsig_pred_.col(i);
  }

  //predicted state covariance matrix
  P.fill(0.0);
  for (int i = 0; i < 2 * n_aug_ + 1; i++) {  //iterate over sigma points

    // state difference
    VectorXd x_diff = Xsig_pred_.col(i) - x;
    //angle normalization
    while (x_diff(3)> M_PI) x_diff(3)-=2.*M_PI;
    while (x_diff(3)<-M_PI) x_diff(3)+=2.*M_PI;

    P = P + weights(i) * x_diff * x_diff.transpose() ;
  }

  //print result
//  std::cout << "Predicted state" << std::endl;
//  std::cout << x << std::endl;
//  std::cout << "Predicted covariance matrix" << std::endl;
//  std::cout << P << std::endl;

  //write result
  *x_out = x;
  *P_out = P;
  *weights_out = weights;
}

/**
 * Updates the state and the state covariance matrix using a laser measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateLidar(MeasurementPackage meas_package) {
  /**
  Use lidar data to update the belief about the object's
  position. Modify the state vector, x_, and covariance, P_.

  You'll also need to calculate the lidar NIS.
  */

  VectorXd z = meas_package.raw_measurements_;

  /*****************************************************************************
   *  Predict Lidar Measurement
   ****************************************************************************/

   //create matrix for sigma points in measurement space
   MatrixXd Zsig = MatrixXd(n_z_laser_, 2 * n_aug_ + 1);

   //transform sigma points into measurement space
   for (int i = 0; i < 2 * n_aug_ + 1; i++) {  //2n+1 simga points

     // extract values for better readibility
     double p_x = Xsig_pred_(0,i);
     double p_y = Xsig_pred_(1,i);


     // measurement model
     Zsig(0,i) = p_x;                        //px
     Zsig(1,i) = p_y;                                 //py

   }

   //mean predicted measurement
   VectorXd z_pred = VectorXd(n_z_laser_);

   z_pred.fill(0.0);
   for (int i=0; i < 2*n_aug_+1; i++) {
       z_pred = z_pred + weights_(i) * Zsig.col(i);
   }

   //measurement covariance matrix S
   MatrixXd S = MatrixXd(n_z_laser_,n_z_laser_);
   S.fill(0.0);
   for (int i = 0; i < 2 * n_aug_ + 1; i++) {  //2n+1 simga points
     //residual
     VectorXd z_diff = Zsig.col(i) - z_pred;

     S = S + weights_(i) * z_diff * z_diff.transpose();
   }

   //add measurement noise covariance matrix
   MatrixXd R = MatrixXd(n_z_laser_,n_z_laser_);
   R <<    std_laspx_*std_laspx_, 0,
           0, std_laspy_*std_laspy_;
   S = S + R;

   //print result
//   std::cout << "z_pred: " << std::endl << z_pred << std::endl;
//   std::cout << "S: " << std::endl << S << std::endl;


   /*****************************************************************************
    *  Update State based on Lidar Measurement
    ****************************************************************************/

   //create matrix for cross correlation Tc
   MatrixXd Tc = MatrixXd(n_x_, n_z_laser_);

   //calculate cross correlation matrix
   Tc.fill(0.0);
   for (int i = 0; i < 2 * n_aug_ + 1; i++) {  //2n+1 simga points

     //residual
     VectorXd z_diff = Zsig.col(i) - z_pred;

     // state difference
     VectorXd x_diff = Xsig_pred_.col(i) - x_;
     //angle normalization
     while (x_diff(3)> M_PI) x_diff(3)-=2.*M_PI;
     while (x_diff(3)<-M_PI) x_diff(3)+=2.*M_PI;

     Tc = Tc + weights_(i) * x_diff * z_diff.transpose();
   }

   //Kalman gain K;
   MatrixXd K = Tc * S.inverse();

   //residual
   VectorXd z_diff = z - z_pred;

   //update state mean and covariance matrix
   x_ = x_ + K * z_diff;
   P_ = P_ - K*S*K.transpose();

   //print result
//   std::cout << "Updated state x: " << std::endl << x_ << std::endl;
//   std::cout << "Updated state covariance P: " << std::endl << P_ << std::endl;

   /*****************************************************************************
    *  NIS of Lidar Measurement
    ****************************************************************************/
   // Chi-Square 95-percentile  Probability for Lidar with 2 degrees of freedom is 5.991

   //calculate NIS value
   NIS_laser_ = z_diff.transpose()*S.inverse()*z_diff;

   //print result
//   std::cout << "NIS_laser: " << std::endl << NIS_laser_ << std::endl;
}

/**
 * Updates the state and the state covariance matrix using a radar measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateRadar(MeasurementPackage meas_package) {
  /**
  Use radar data to update the belief about the object's
  position. Modify the state vector, x_, and covariance, P_.

  You'll also need to calculate the radar NIS.
  */
  VectorXd z = meas_package.raw_measurements_;
  /*****************************************************************************
   *  Predict Radar Measurement
   ****************************************************************************/

   //create matrix for sigma points in measurement space
   MatrixXd Zsig = MatrixXd(n_z_radar_, 2 * n_aug_ + 1);

   //transform sigma points into measurement space
   for (int i = 0; i < 2 * n_aug_ + 1; i++) {  //2n+1 simga points

     // extract values for better readibility
     double p_x = Xsig_pred_(0,i);
     double p_y = Xsig_pred_(1,i);
     double v  = Xsig_pred_(2,i);
     double yaw = Xsig_pred_(3,i);

     double v1 = cos(yaw)*v;
     double v2 = sin(yaw)*v;

     // measurement model
     Zsig(0,i) = sqrt(p_x*p_x + p_y*p_y);                        //r
     Zsig(1,i) = atan2(p_y,p_x);                                 //phi
     Zsig(2,i) = (p_x*v1 + p_y*v2 ) / sqrt(p_x*p_x + p_y*p_y);   //r_dot
   }

   //mean predicted measurement
   VectorXd z_pred = VectorXd(n_z_radar_);

   z_pred.fill(0.0);
   for (int i=0; i < 2*n_aug_+1; i++) {
       z_pred = z_pred + weights_(i) * Zsig.col(i);
   }

   //measurement covariance matrix S
   MatrixXd S = MatrixXd(n_z_radar_,n_z_radar_);
   S.fill(0.0);
   for (int i = 0; i < 2 * n_aug_ + 1; i++) {  //2n+1 simga points
     //residual
     VectorXd z_diff = Zsig.col(i) - z_pred;

     //angle normalization
     while (z_diff(1)> M_PI) z_diff(1)-=2.*M_PI;
     while (z_diff(1)<-M_PI) z_diff(1)+=2.*M_PI;

     S = S + weights_(i) * z_diff * z_diff.transpose();
   }

   //add measurement noise covariance matrix
   MatrixXd R = MatrixXd(n_z_radar_,n_z_radar_);
   R <<    std_radr_*std_radr_, 0, 0,
           0, std_radphi_*std_radphi_, 0,
           0, 0,std_radrd_*std_radrd_;
   S = S + R;

   //print result
//   std::cout << "z_pred: " << std::endl << z_pred << std::endl;
//   std::cout << "S: " << std::endl << S << std::endl;


   /*****************************************************************************
    *  Update State based on Radar Measurement
    ****************************************************************************/

   //create matrix for cross correlation Tc
   MatrixXd Tc = MatrixXd(n_x_, n_z_radar_);

   //calculate cross correlation matrix
   Tc.fill(0.0);
   for (int i = 0; i < 2 * n_aug_ + 1; i++) {  //2n+1 simga points

     //residual
     VectorXd z_diff = Zsig.col(i) - z_pred;
     //angle normalization
     while (z_diff(1)> M_PI) z_diff(1)-=2.*M_PI;
     while (z_diff(1)<-M_PI) z_diff(1)+=2.*M_PI;

     // state difference
     VectorXd x_diff = Xsig_pred_.col(i) - x_;
     //angle normalization
     while (x_diff(3)> M_PI) x_diff(3)-=2.*M_PI;
     while (x_diff(3)<-M_PI) x_diff(3)+=2.*M_PI;

     Tc = Tc + weights_(i) * x_diff * z_diff.transpose();
   }

   //Kalman gain K;
   MatrixXd K = Tc * S.inverse();

   //residual
   VectorXd z_diff = z - z_pred;

   //angle normalization
   while (z_diff(1)> M_PI) z_diff(1)-=2.*M_PI;
   while (z_diff(1)<-M_PI) z_diff(1)+=2.*M_PI;

   //update state mean and covariance matrix
   x_ = x_ + K * z_diff;
   P_ = P_ - K*S*K.transpose();

   //print result
//   std::cout << "Updated state x: " << std::endl << x_ << std::endl;
//   std::cout << "Updated state covariance P: " << std::endl << P_ << std::endl;

   /*****************************************************************************
    *  NIS of Radar Measurement
    ****************************************************************************/
   // Chi-Square 95-percentile  Probability for Radar with 3 degrees of freedom is 7.815
   //calculate NIS value
   NIS_radar_ = z_diff.transpose()*S.inverse()*z_diff;


   //print result
//   std::cout << "NIS_radar: " << std::endl << NIS_radar_ << std::endl;

}








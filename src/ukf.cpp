#include "ukf.h"
#include "Eigen/Dense"
#include <iostream>
#include "tools.h"

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::vector;

/**
 * Initializes Unscented Kalman filter
 */
UKF::UKF() {
  // if this is false, laser measurements will be ignored (except during init)
  use_laser_ = true;

  // if this is false, radar measurements will be ignored (except during init)
  use_radar_ = true;

  // initial state vector
  x_ = VectorXd(5);

  // initial covariance matrix
  P_ = MatrixXd(5, 5);

  // Process noise standard deviation longitudinal acceleration in m/s^2
  std_a_ = 2;

  // Process noise standard deviation yaw acceleration in rad/s^2
  std_yawdd_ = 0.5;

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

  /**
  TODO:

  Complete the initialization. See ukf.h for other member properties.

  Hint: one or more values initialized above might be wildly off...
  */
    is_initialized_ = false;
    time_us_ = 0;
    
    ///* State dimension
    n_x_ = x_.size();
    
    ///* Augmented state dimension
    n_aug_ = n_x_ + 2;
    
    ///* Sigma point spreading parameter
    lambda_ = 3 - n_aug_;
    
    ///* Weights of sigma points
    weights_ = VectorXd(2*n_aug_+1);
    
    Xsig_pred_ = MatrixXd(n_x_, 2*n_aug_+1);;
    
}

UKF::~UKF() {}

/**
 * @param {MeasurementPackage} meas_package The latest measurement data of
 * either radar or laser.
 */
void UKF::ProcessMeasurement(MeasurementPackage measurement_pack) {
  /**
  TODO:

  Complete this function! Make sure you switch between lidar and radar
  measurements.
  */
    /*****************************************************************************
     *  Initialization
     ****************************************************************************/
    if (!is_initialized_) {
        // first measurement
        cout << "UKF: " << endl;
        x_ << 1, 1, 1, 1, 1;
        
        if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
            cout << "Initialize with radar" << endl;
            /**
             Convert radar from polar to cartesian coordinates and initialize state.
             */
            float rho = measurement_pack.raw_measurements_[0];
            float phi = measurement_pack.raw_measurements_[1];
            x_ << rho*cos(phi), rho*sin(phi), 0, 0, 0;
            
        }
        else if (measurement_pack.sensor_type_ == MeasurementPackage::LASER) {
            cout << "initialize with lidar" << endl;
            /**
             Initialize state.
             */
            x_ << measurement_pack.raw_measurements_[0], measurement_pack.raw_measurements_[1], 0, 0, 0;
            
        }
        
        P_ << 1, 0, 0, 0, 0,
        0, 1, 0, 0, 0,
        0, 0, 1, 0, 0,
        0, 0, 0, 1, 0,
        0, 0, 0, 0, 1;

        
        // done initializing, no need to predict or update
        time_us_ = measurement_pack.timestamp_;
        is_initialized_ = true;
        return;
    }
    
    float delta_t = (measurement_pack.timestamp_ - time_us_) / 1000000.0;	//dt - expressed in seconds
    time_us_ = measurement_pack.timestamp_;
    cout << "delta_t: " << delta_t << endl;
    
    cout << "starting prediction!" << endl;
    
    // a fix for when delta_t is too large -- break it down
    while (delta_t > 0.1)
    {
        const double dt = 0.05;
        // prediction
        
        Prediction(dt);
        
        // update value
        delta_t -= dt;
    }
    Prediction(delta_t);
    
    if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
        UpdateRadar(measurement_pack);
    }
    else if (measurement_pack.sensor_type_ == MeasurementPackage::LASER) {
        UpdateLidar(measurement_pack);
    }
}

/**
 * Predicts sigma points, the state, and the state covariance matrix.
 * @param {double} delta_t the change in time (in seconds) between the last
 * measurement and this one.
 */
void UKF::Prediction(double delta_t) {
  /**
  TODO:

  Complete this function! Estimate the object's location. Modify the state
  vector, x_. Predict sigma points, the state, and the state covariance matrix.
  */
    
    // generate sigma points
    //create augmented state covariance
    MatrixXd P_aug = MatrixXd(n_aug_, n_aug_);
    P_aug.fill(0.0);
    VectorXd x_aug = VectorXd(n_aug_);
    
    //create sigma point matrix
    MatrixXd Xsig_aug = MatrixXd(n_aug_, 2 * n_aug_ + 1);
    Xsig_aug.fill(0);
    
    //create augmented mean state
    x_aug.head(5) = x_;
    cout << "x_aug: " << x_aug << endl;
    Xsig_aug.col(0) = x_aug;
    //create augmented covariance matrix
    P_aug.topLeftCorner(5, 5) = P_;
    //create square root matrix
    MatrixXd Q = MatrixXd(2, 2);
    Q << std_a_*std_a_, 0,
    0, std_yawdd_*std_yawdd_;
    P_aug.bottomRightCorner(2, 2) = Q;
    cout << "P_aug: " << P_aug << endl;

    MatrixXd A = P_aug.llt().matrixL();
    cout << "A: " << A << endl;
    for (int i=0; i<n_aug_; i++) {
        Xsig_aug.col(i+1) = x_aug + sqrt(lambda_+n_aug_)*A.col(i);
        Xsig_aug.col(i+1+n_aug_) = x_aug - sqrt(lambda_ + n_aug_)*A.col(i);
    }
    cout << "Xsig_aug: " << Xsig_aug << endl;

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
        Xsig_pred_(0,i) = px_p;
        Xsig_pred_(1,i) = py_p;
        Xsig_pred_(2,i) = v_p;
        Xsig_pred_(3,i) = yaw_p;
        Xsig_pred_(4,i) = yawd_p;
    }
    cout << "Xsig_pred_: " << Xsig_pred_ << endl;
    
    // final predictions
    //set weights
    weights_(0) = lambda_ / (lambda_ + n_aug_);
    for (int i=1; i<2*n_aug_+1; i++) {
        weights_(i) = 0.5/(lambda_ + n_aug_);
    }
    cout << "weights_" << weights_ << endl;
    
    //predict state mean
    x_.fill(0.0);
    for (int i=0; i<2*n_aug_+1; i++) {
        x_ = x_ + weights_(i)*Xsig_pred_.col(i);
    }
    
    P_.fill(0.0);
    //predict state covariance matrix
    for (int i=0; i<2*n_aug_+1; i++) {
        VectorXd diff = Xsig_pred_.col(i) - x_;
        //angle normalization
        diff(3) = tools.NormalizeAngle(diff(3));
        P_ = P_ + weights_(i)*diff*diff.transpose();
    }
    
    cout << "x_" << x_ << endl;
    cout << "P_: " << P_ << endl;
}

/**
 * Updates the state and the state covariance matrix using a laser measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateLidar(MeasurementPackage meas_package) {
  /**
  TODO:

  Complete this function! Use lidar data to update the belief about the object's
  position. Modify the state vector, x_, and covariance, P_.

  You'll also need to calculate the lidar NIS.
  */
    cout << "Update with lidar" << endl;
    VectorXd z = VectorXd(2);
    z << meas_package.raw_measurements_(0),
        meas_package.raw_measurements_(1);
    
    //measurement covariance matrix S
    MatrixXd H_ = MatrixXd(2, 5);
    H_.fill(0.0);
    H_(0, 0) = 1.0;
    H_(1, 1) = 1.0;
    
    MatrixXd R_ = MatrixXd(2, 2);
    R_.fill(0.0);
    R_(0, 0) = std_laspx_*std_laspx_;
    R_(1, 1) = std_laspy_*std_laspy_;
    
    //cout << "H_*x_: " << H_*x_ << endl;
    //cout << "z: " << z << endl;
    VectorXd y = z - H_*x_;
    //cout << "y: " << y << endl;
    MatrixXd Ht = H_.transpose();
    MatrixXd S = H_ * P_ * Ht + R_;
    //cout << "S: " << S << endl;
    
    MatrixXd Si = S.inverse();
    //cout << "Si: " << Si << endl;
    MatrixXd PHt = P_ * Ht;
    //cout << "PHt: " << PHt << endl;
    MatrixXd K = PHt * (S.inverse());
    //cout << "K: " << K << endl;
    
    //new estimate
    x_ = x_ + (K * y);
    MatrixXd I = MatrixXd::Identity(5, 5);
    P_ = (I - K * H_) * P_;

    cout << x_ << endl;
    cout << P_ << endl;
}

/**
 * Updates the state and the state covariance matrix using a radar measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateRadar(MeasurementPackage meas_package) {
    cout << "Update with Radar" << endl;
  /**
  TODO:

  Complete this function! Use radar data to update the belief about the object's
  position. Modify the state vector, x_, and covariance, P_.

  You'll also need to calculate the radar NIS.
  */
    int n_z = 3;
    //create matrix for sigma points in measurement space
    MatrixXd Zsig = MatrixXd(n_z, 2 * n_aug_ + 1);
    
    //mean predicted measurement
    VectorXd z_pred = VectorXd(n_z);
    
    //measurement covariance matrix S
    MatrixXd S = MatrixXd(n_z,n_z);
    
    //transform sigma points into measurement space
    for (int i=0; i<2*n_aug_+1; i++) {
        VectorXd sig_pt = Xsig_pred_.col(i);
        double px = sig_pt(0);
        double py = sig_pt(1);
        double v = sig_pt(2);
        double yaw = sig_pt(3);
        double yawd = sig_pt(4);
        if (px == 0 && py == 0) {
            // check for zero division
            Zsig.col(i) << 0, 0, 0;
        } else {
            Zsig.col(i) << sqrt(px*px+py*py),
            atan2(py, px),
            (px*cos(yaw)*v + py*sin(yaw)*v)/ sqrt(px*px+py*py);
        }
    }
    //calculate mean predicted measurement
    z_pred.fill(0.0);
    for (int i=0; i<2*n_aug_+1; i++) {
        z_pred = z_pred + weights_(i)*Zsig.col(i);
    }
    
    //calculate measurement covariance matrix S
    S.fill(0.0);
    for (int i=0; i<2*n_aug_+1; i++) {
        VectorXd diff = Zsig.col(i) - z_pred;
        //angle normalization
        diff(1) = tools.NormalizeAngle(diff(1));
        S = S + weights_(i)*diff*diff.transpose();
    }
    MatrixXd R = MatrixXd(3, 3);
    R << std_radr_*std_radr_, 0, 0,
    0, std_radphi_*std_radphi_, 0,
    0, 0, std_radrd_*std_radrd_;
    S = S + R;

    
    // UKF Update
    VectorXd z = meas_package.raw_measurements_;
    
    //create matrix for cross correlation Tc
    MatrixXd Tc = MatrixXd(n_x_, n_z);
    Tc.fill(0.0);
    
    //calculate cross correlation matrix
    for (int i=0; i<2*n_aug_+1; i++) {
        
        VectorXd z_diff = Zsig.col(i) - z_pred;
        z_diff(1) = tools.NormalizeAngle(z_diff(1));
        
        VectorXd x_diff = Xsig_pred_.col(i) - x_;
        x_diff(3) = tools.NormalizeAngle(x_diff(3));
        Tc = Tc + weights_(i)*x_diff*z_diff.transpose();
    }
    
    //calculate Kalman gain K;
    MatrixXd K = Tc*S.inverse();
    
    //update state mean and covariance matrix
    VectorXd z_diff = z - z_pred;
    z_diff(1) = tools.NormalizeAngle(z_diff(1));
    x_ = x_ + K*z_diff;
    P_ = P_ - K*S*K.transpose();
    
    cout << x_ << endl;
    cout << P_ << endl;
}

#include "ukf.h"
#include "Eigen/Dense"
#include <iostream>
#include "tools.h"
using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::vector;

#define EPS 0.001 //small error

/**
 * Initializes Unscented Kalman filter
 */
UKF::UKF() {

    /*************************************************************
     * Parameters initialization
     ************************************************************/
  // if this is false, laser measurements will be ignored (except during init)
  use_laser_ = false;

  // if this is false, radar measurements will be ignored (except during init)
  use_radar_ = true;

  // initial state vector
  x_ = VectorXd(5);

  // initial covariance matrix
  P_ = MatrixXd(5, 5);

  // Process noise standard deviation longitudinal acceleration in m/s^2
  std_a_ = 1.5;

  // Process noise standard deviation yaw acceleration in rad/s^2
  std_yawdd_ = 0.57;

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


  // Filter initialization tag
  is_initialized_=false;

  // state dimension
  n_x_=x_.size();

  // Augmented state dimension
  n_aug_=n_x_+2;

  //sigma point spreading parameters
  lambda_=3-n_aug_;

  // Measurement dimensions
    n_z_radar=3;
    n_z_laser=2;

  // Weights of sigma points
  weights_=VectorXd(2*n_aug_+1);
  /*weights_.fill(0.0);
  weights_(0)=lambda_/(lambda_+n_aug_);
  weights_.tail(2*n_aug_).fill(0.5/(lambda_+n_aug_));
  cout<<weights_<<endl;*/
  // Initialize NIS for both sensors
  NIS_LASER=0;
  NIS_RADAR=0;

  // Initialize augmented vector for prediction
  Xsig_pred_=MatrixXd(n_x_,2*n_aug_+1);
  Xsig_pred_.fill(0.0);

  // Measurement covaraince Matrix
  R_laser_=MatrixXd(n_z_laser,n_z_laser);
  R_laser_<<std_laspx_*std_laspx_,0,
          0,std_laspy_*std_laspy_;
  R_radar_=MatrixXd(n_z_radar,n_z_radar);
  R_radar_<<std_radr_*std_radr_,0,0,
          0,std_radphi_*std_radphi_,0,
          0,0,std_radrd_*std_radrd_;

}

UKF::~UKF() {}

/**
 * @param {MeasurementPackage} meas_package The latest measurement data of
 * either radar or laser.
 */
void UKF::ProcessMeasurement(MeasurementPackage meas_package) {

  /*****************************************************************
   * Check for initialization
   ****************************************************************/
  if (!is_initialized_){
    is_initialized_= true;

    //check for the sensor in use
    if (meas_package.sensor_type_==MeasurementPackage::LASER){
      // initialize state vector using LIDAR
      x_<< meas_package.raw_measurements_[0],
           meas_package.raw_measurements_[1],
            0.0,
            0.0,
            0.0;
      // for small values
      if(x_(0)<EPS and x_(1)<EPS ){
        x_(0)=EPS;
        x_(1)=EPS;
      }

      }
    else{

      // initialize the state vector using RADAR measurements

      double rho=meas_package.raw_measurements_[0];
      double angel=meas_package.raw_measurements_[1];
      double rho_dot=meas_package.raw_measurements_[2];
      double p_x=rho*cos(angel);
      double p_y=rho*sin(angel);
      double v_x=rho_dot*cos(angel);
      double v_y=rho_dot*sin(angel);
      double v=sqrt(v_x*v_x + v_y*v_y);

      x_<< p_x,p_y,v,0,0;
    }

    // initialize covariance matrix
    P_.fill(0.0);
    P_(0,0)=1;
    P_(1,1)=1;
    P_(2,2)=1;
    P_(3,3)=1;
    P_(4,4)=1;


    // initialize the weights

    weights_.fill(0.0);
    weights_(0)=lambda_/(lambda_+n_aug_);
    weights_.tail(2*n_aug_).fill(0.5/(lambda_+n_aug_));

    // initialize time stamp
    time_us_=meas_package.timestamp_;
    //cout<<x_<<endl;

    return;
  }

  /****************************************************************
   * Prediction
   **************************************************************/

  double dt=(meas_package.timestamp_ - time_us_)/1000000.0;
  time_us_=meas_package.timestamp_;
  //cout<<is_initialized_<<endl;
  Prediction(dt);


  /*************************************************************
   * Update
   ************************************************************/

  if(meas_package.sensor_type_==MeasurementPackage::LASER && use_laser_) {
    //LASER update
    //cout<<"The sensor is LIDAR"<<endl;
    UpdateLidar(meas_package);
    //cout<<"The state is updated using a LIDAR measurements"<<endl;

  }
  else if (meas_package.sensor_type_==MeasurementPackage::RADAR && use_radar_) {
    //RADAR update
    UpdateRadar(meas_package);
  }//cout<<x_;

  }

/**
 * Get sigma points
 */
void UKF::GetSigmaPoints(double dt) {
  // State augmentation for process noise components
  //cout<<"Getting sigma points started"<<endl;
  VectorXd x_aug=VectorXd(n_aug_);
  x_aug.fill(0.0);
  x_aug.head(n_x_)=x_;

  // Covariance Matrix augmentation for process noise covariances
  MatrixXd Q=MatrixXd(2,2);
  Q<<std_a_*std_a_,0,
     0,std_yawdd_*std_yawdd_;
  MatrixXd P_aug=MatrixXd(n_aug_,n_aug_);
  P_aug.fill(0.0);
  P_aug.topLeftCorner(n_x_,n_x_)=P_;
  P_aug(5, 5) = std_a_ * std_a_;
  P_aug(6, 6) = std_yawdd_ * std_yawdd_;

  // Generate sigma point matrix
  MatrixXd Xsig_aug=MatrixXd(n_aug_,2*n_aug_+1);
    //Xsig_aug.fill(0.0);
  MatrixXd A=P_aug.llt().matrixL(); //Chelosky decomposition
  Xsig_aug.col(0)=x_aug;
  double sqrt_lambda=sqrt(lambda_+n_aug_);
  for(int i=0;i<n_aug_;i++){
      Xsig_aug.col(i+1)=x_aug+sqrt_lambda*A.col(i);
      Xsig_aug.col(i+1+n_aug_)=x_aug-sqrt_lambda*A.col(i);
  }

  // Predict sigma points
  for(int i=0;i<2*n_aug_+1;i++){

    // temporary variables to be used
    double p_x=Xsig_aug(0,i);
    double p_y=Xsig_aug(1,i);
    double v=Xsig_aug(2,i);
    double yaw=Xsig_aug(3,i);
    double yaw_dot=Xsig_aug(4,i);
    double nu_a=Xsig_aug(5,i);
    double nu_yaw_a=Xsig_aug(6,i);

    if (fabs(p_x) < EPS && fabs(p_y) < EPS) {
      p_x = EPS;
      p_y = EPS;
    }

    double px_pred,py_pred;
    // check for straight line motion
    if(fabs(yaw_dot)>EPS){
      px_pred=p_x+v/yaw_dot * (sin(yaw+yaw_dot*dt)-sin(yaw));
      py_pred=p_y+v/yaw_dot * (cos(yaw)-cos(yaw+yaw_dot*dt));
    }
    else{
      px_pred=p_x+v*dt*cos(yaw);
      py_pred=p_y+v*dt*sin(yaw);
    }


    double v_p = v;
    double yaw_p = yaw+yaw_dot*dt;
    double yawd_p = yaw_dot;


    //process noise Handling
    px_pred+=0.5*nu_a*dt*dt*cos(yaw);
    py_pred+=0.5*nu_a*dt*dt*sin(yaw);
    v_p+=nu_a*dt;
    yaw_p+=0.5*nu_yaw_a*dt*dt;
    yawd_p+=nu_yaw_a*dt;

    //predicted sigma points
    Xsig_pred_(0,i)=px_pred;
    Xsig_pred_(1,i)=py_pred;
    Xsig_pred_(2,i)=v_p;
    Xsig_pred_(3,i)=yaw_p;
    Xsig_pred_(4,i)=yawd_p;
  }
  //cout<<Xsig_pred_<<endl;
  //cout<<"Getting sigma points ended"<<endl;

}
/**
 * Predicts sigma points, the state, and the state covariance matrix.
 * @param {double} delta_t the change in time (in seconds) between the last
 * measurement and this one.
 */
void UKF::Prediction(double delta_t) {

    //cout<<"Prediction Testing begining"<<endl;
  // Get predicted sigma points
    GetSigmaPoints(delta_t);
  //cout<<Xsig_pred_<<endl;
  // predict state mean
  VectorXd x=VectorXd(n_x_);
  x.fill(0.0);
  x=Xsig_pred_*weights_;
  /*for(int i=0;i<2*n_aug_+1;i++){

    x=x+weights_(i)*Xsig_pred_.col(i);

  }*/

  //predict state covariance matrix
  MatrixXd P=MatrixXd(n_x_,n_x_);
  P.fill(0.0);

  for(int i=0;i<2*n_aug_+1;i++){
      VectorXd x_diff=Xsig_pred_.col(i)-x;

      //normalize angel
      Normalize_angel(x_diff(3));

      P=P+weights_(i)*x_diff*x_diff.transpose();
  }

  //debuging
  //cout<<"temp x is "<<x<<endl;
  x_=x;
  P_=P;

  //debuging
    //cout<<"P is"<<P_<<endl;

    //cout<<"prediction Testing Ending"<<endl;
}

/**
 * Updates the state and the state covariance matrix using a laser measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateLidar(MeasurementPackage meas_package) {

    //cout<<"LIDAR Testing begining"<<endl;
  // Use only the first two rows from sigma points
  MatrixXd z_sigma=Xsig_pred_.block(0,0,n_z_laser,2*n_aug_+1);
  VectorXd z_pred=z_sigma*weights_;

  MatrixXd S = MatrixXd(n_z_laser,n_z_laser);
  S.fill(0.0);

  for(int i=0;i<2*n_aug_+1;i++){
    VectorXd innov=(z_sigma.col(i)-z_pred);
    Normalize_angel(innov(1));
    S=S+weights_(i)*innov*innov.transpose();
  }






  S=S+R_laser_;
  //cout<<"Testing"<<endl;
  MatrixXd Tc = MatrixXd(n_x_, n_z_laser);
  Tc.fill(0.0);
  for(int i=0;i<2*n_aug_+1;i++){
    VectorXd x_diff= Xsig_pred_.col(i)-x_;
    Normalize_angel(x_diff(3));
    VectorXd innov=(z_sigma.col(i)-z_pred);

    Tc=Tc+weights_(i)*x_diff*innov.transpose();
  }

  //calculate Kalman gain K;
  MatrixXd K=Tc*S.inverse();
  Tools tools_;
  //cout<<"LASER NIS : "<<tools_.Nis(meas_package.raw_measurements_,z_pred,S)<<endl;
  //update state mean and covariance matrix
  x_=x_+K*(meas_package.raw_measurements_-z_pred);
    //cout<<"mean predicted"<<endl;
  P_=P_-K*S*K.transpose();

  cout<<"X = "<<x_<<endl;
  cout<<"P = "<<P_<<endl;

}

/**
 * Updates the state and the state covariance matrix using a radar measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateRadar(MeasurementPackage meas_package) {
    //cout<<"RADAR Testing begining"<<endl;

  MatrixXd z_sigma=MatrixXd(n_z_radar,2*n_aug_+1);
    //z_sigma.fill(0.0);
  for(int i=0;i<2*n_aug_+1;i++){
      double p_x=Xsig_pred_(0,i);
      double p_y=Xsig_pred_(1,i);
      double v=Xsig_pred_(2,i);
      double epsi=Xsig_pred_(3,i);
      z_sigma(0,i)=sqrt(p_x*p_x + p_y*p_y);  //rho
      z_sigma(1,i)=atan2(p_y,p_x);           // phi
      z_sigma(2,i)=(p_x*cos(epsi)*v + p_y*sin(epsi)*v)/z_sigma(0,i);    //rho_dot
  }

  // mean predicted measurements
    VectorXd z_pred=VectorXd(n_z_radar);
    //z_pred.fill(0.0);
  z_pred=z_sigma*weights_;

    //measurement covariance matrix S
    MatrixXd S = MatrixXd(n_z_radar,n_z_radar);
    S.fill(0.0);
    //VectorXd innov=VectorXd(z_pred.size());
    //innov.fill(0.0);
    for(int i=0;i<2*n_aug_+1;i++){
      VectorXd innov=(z_sigma.col(i)-z_pred);

      //Normalization
      Normalize_angel(innov(1));
        S=S+weights_(i)*innov*innov.transpose();
    }
    S=S+R_radar_;

  MatrixXd Tc = MatrixXd(n_x_, n_z_radar);
  Tc.fill(0.0);
    for(int i=0;i<2*n_aug_+1;i++){
      VectorXd x_diff=Xsig_pred_.col(i)-x_;
      VectorXd innov=(z_sigma.col(i)-z_pred);
        Normalize_angel(x_diff(3));
        Normalize_angel(innov(1));
        Tc=Tc+weights_(i)*x_diff*innov.transpose();
    }

    //calculate Kalman gain K;
    MatrixXd K=Tc*S.inverse();
    VectorXd innov=meas_package.raw_measurements_-z_pred;
    Normalize_angel(innov(1));

    Tools tools_;
   // cout<<"RADAR NIS : "<<tools_.Nis(meas_package.raw_measurements_,z_pred,S)<<endl;
    //update state mean and covariance matrix
    x_=x_+K*innov;
    P_=P_-K*S*K.transpose();
    cout<<"X = "<<x_<<endl;
    cout<<"P = "<<P_<<endl;

}

void UKF::Normalize_angel(double &x) {
    while(x>M_PI)
        x-=2.*M_PI;
    while(x<-M_PI)
        x+=2.*M_PI;
}
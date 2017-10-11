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
	// if this is false, laser measurements will be ignored (except during init)
	use_laser_ = true;

	// if this is false, radar measurements will be ignored (except during init)
	use_radar_ = true;

	///* State dimension
	n_x_ = 5;

	///* Augmented state dimension
	n_aug_ = 7;

	// initial state vector
	x_ = VectorXd(n_x_);

	// state vector: [pos1 pos2 vel_abs yaw_angle yaw_rate nu psi] in SI units and rad
	x_aug_ = VectorXd(n_aug_);

	// initial covariance matrix
	P_ = MatrixXd(n_x_, n_x_);

	// aug state covariance matrix
	P_aug_ = MatrixXd(n_aug_, n_aug_);

	// Process noise standard deviation longitudinal acceleration in m/s^2
	std_a_ = 1.0;

	// Process noise standard deviation yaw acceleration in rad/s^2
	std_yawdd_ = 1.0;

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

	///* var to capture previous timestamp
	previous_timestamp = 0;

	///* predicted sigma points matrix
	Xsig_pred_ = MatrixXd(n_x_, 2 * n_aug_ + 1);

	///* Sigma point spreading parameter
	lambda_ = 3 - n_aug_;

	// Measurement noise covariance matrices initialization
	R_radar_ = MatrixXd(3, 3);
	R_radar_ << std_radr_*std_radr_, 0, 0,
			0, std_radphi_*std_radphi_, 0,
			0, 0,std_radrd_*std_radrd_;
	R_lidar_ = MatrixXd(2, 2);
	R_lidar_ << std_laspx_*std_laspx_,0,
			0,std_laspy_*std_laspy_;

	//weights initialization
	weights_ = VectorXd(2*n_aug_+1);
	const float lambda_plus_n_aug_ = lambda_+n_aug_;
	weights_(0) = lambda_ / lambda_plus_n_aug_;
	float val = float(1/(2*lambda_plus_n_aug_));
	for(int i=1; i<2*n_aug_+1; i++){
		weights_(i) = val;
	}
}

UKF::~UKF() {}

/**
 * @param {MeasurementPackage} meas_package The latest measurement data of
 * either radar or laser.
 */
void UKF::ProcessMeasurement(MeasurementPackage meas_package) {
	cout << "ProcessMeasurement: " << endl;
	if(!is_initialized_){
		// first measurement
		cout << "Enter init: " << endl;
		x_.fill(0.0);

		//init state co-variance matrix
		//We can initialize this P based on measurement std_deviations, for now start with identity matrix
		P_ <<  1, 0, 0, 0, 0,
				0, 1, 0, 0, 0,
				0, 0, 1, 0, 0,
				0, 0, 0, 1, 0,
				0, 0, 0, 0, 1;

		if (meas_package.sensor_type_ == MeasurementPackage::RADAR) {
			//Convert radar from polar to cartesian coordinates and initialize state.
			float rho = meas_package.raw_measurements_(0);
			float theata = meas_package.raw_measurements_(1);
			float rho_dot = meas_package.raw_measurements_(2);
			double vx = rho_dot * cos(theata);
			double vy = rho_dot * sin(theata);
			x_(0) = rho*cos(theata);
			x_(1) = rho*sin(theata);
			x_(2) = sqrt(vx*vx + vy*vy);
		}
		else if (meas_package.sensor_type_ == MeasurementPackage::LASER) {
	      //Initialize state.
			x_(0) = meas_package.raw_measurements_(0);
			x_(1) = meas_package.raw_measurements_(1);
		}

		// Check for zeros in the initial values of px and py
		if (fabs(x_(0)) < 0.001 and fabs(x_(1)) < 0.001) {
			x_(0) = 0.001;
			x_(1) = 0.001;
		}

		previous_timestamp = meas_package.timestamp_;
		// done initializing, no need to predict or update
		is_initialized_ = true;
		std::cout<<"INIT Complete"<<std::endl;
		return;
	}

	float dt = (meas_package.timestamp_ - previous_timestamp) / 1000000.0;	//dt - expressed in seconds
	previous_timestamp = meas_package.timestamp_;

	//Predict the state for dt
	Prediction(dt);

	if (meas_package.sensor_type_ == MeasurementPackage::RADAR) {
		// Radar updates
		UpdateRadar(meas_package);
	} else {
		// Laser updates
		UpdateLidar(meas_package);
	}
}

/**
 * Predicts sigma points, the state, and the state covariance matrix.
 * @param {double} delta_t the change in time (in seconds) between the last
 * measurement and this one.
 */
void UKF::Prediction(double delta_t) {
	std::cout<<"Prediction start, deltaT = "<<delta_t<<std::endl;

	//create sigma point matrix
	MatrixXd Xsig_aug = MatrixXd(n_aug_, 2 * n_aug_ + 1);

	//update augmented mean state
	x_aug_.fill(0.0);
	x_aug_.head(n_x_) = x_;

	//update augmented covariance matrix
	P_aug_.fill(0.0);
	P_aug_.topLeftCorner(n_x_, n_x_) = P_;
	P_aug_.bottomRightCorner(2, 2) << std_a_*std_a_, 0,
			0, std_yawdd_*std_yawdd_;

	//create square root matrix
	MatrixXd B = P_aug_.llt().matrixL();

	//create augmented sigma points
	//First column is same as x_aug_
	Xsig_aug.col(0)  = x_aug_;

	//set remaining sigma points
	for (int i = 0; i < n_aug_; i++)
	{
		Xsig_aug.col(i+1)     = x_aug_ + sqrt(lambda_+n_aug_) * B.col(i);
		Xsig_aug.col(i+1+n_aug_) = x_aug_ - sqrt(lambda_+n_aug_) * B.col(i);
	}

	//predict sigma points
	for(int i=0; i<2*n_aug_+1; i++){
		MatrixXd tmp = Xsig_aug.col(i);
		VectorXd x = VectorXd(n_x_);
		x << tmp(0), tmp(1), tmp(2), tmp(3), tmp(4);
		const float v_k = tmp(2);
		const float psi = tmp(3);
		const float psi_dot = tmp(4);
		const float nu_a_k = tmp(5);
		const float nu_psi_dot_dot_k = tmp(6);

		VectorXd state_noise_vector = VectorXd(n_x_);
		VectorXd process_noise_vector = VectorXd(n_x_);
		process_noise_vector<<  (0.5*delta_t*delta_t*cos(psi)*nu_a_k),
				(0.5*delta_t*delta_t*sin(psi)*nu_a_k),
				(delta_t*nu_a_k),
				(0.5*delta_t*delta_t*nu_psi_dot_dot_k),
				(delta_t*nu_psi_dot_dot_k);
		//avoid division by zero when psi_dot is zero
		if(psi_dot != 0){
			state_noise_vector<<((v_k/psi_dot)*(sin(psi+psi_dot*delta_t)-sin(psi))),
					((v_k/psi_dot)*(-cos(psi+psi_dot*delta_t)+cos(psi))),
					0,
					psi_dot*delta_t,
					0;
		}
		else{
			state_noise_vector<<(v_k*cos(psi)*delta_t),
					(v_k*sin(psi)*delta_t),
					0,
					(psi_dot*delta_t),
					0;
		}
		//write predicted sigma points into correct column
		Xsig_pred_.col(i) = x + state_noise_vector + process_noise_vector;
	}

	//Calculate predicted state mean
	x_.fill(0.0);//reset x_
	for(int i=0; i<2*n_aug_+1; i++){
		x_ = x_+ weights_(i)*Xsig_pred_.col(i);
	}

	//Calculate predicted state covariance matrix
	P_.fill(0.0);//reset P_
	for(int i=0; i<2*n_aug_+1; i++){
		VectorXd tmp = Xsig_pred_.col(i) - x_;

		//angle normalization
		while (tmp(3)> M_PI) tmp(3)-=2.*M_PI;
		while (tmp(3)<-M_PI) tmp(3)+=2.*M_PI;

		P_ += (weights_(i)*tmp*tmp.transpose());
	}
}

/**
 * Updates the state and the state covariance matrix using a laser measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateLidar(MeasurementPackage meas_package) {
	std::cout<<"UpdateLidar start"<<std::endl;
	const int n_z = 2;

	//create matrix for sigma points in measurement space
	MatrixXd Zsig = MatrixXd(n_z, 2 * n_aug_ + 1);
	Zsig.fill(0.0);

	//mean predicted measurement
	VectorXd z_pred = VectorXd(n_z);
	z_pred.fill(0.0);

	//measurement covariance matrix S
	MatrixXd S = MatrixXd(n_z,n_z);
	S.fill(0.0);

	for(int i=0; i< 2 * n_aug_ + 1; ++i){
		VectorXd tmp = Xsig_pred_.col(i);
		float px = tmp(0);
		float py = tmp(1);

		Zsig.col(i) << px, py;
	}

	//calculate mean predicted measurement
	for(int i=0; i< 2 * n_aug_ + 1; ++i){
		z_pred += weights_(i) * Zsig.col(i);
	}

	//calculate measurement covariance matrix S
	for(int i=0; i< 2 * n_aug_ + 1; ++i){
		VectorXd tmp = Zsig.col(i) - z_pred;
		S = S + weights_(i)*tmp*tmp.transpose();
	}

	S = S + R_lidar_;

	//create matrix for cross correlation Tc
	MatrixXd Tc = MatrixXd(n_x_, n_z);
	Tc.fill(0.0);


	//calculate cross correlation matrix
	for(int i=0; i<(2 * n_aug_ + 1); ++i){
		VectorXd x_diff = Xsig_pred_.col(i) - x_;
		VectorXd z_diff = Zsig.col(i) - z_pred;
		Tc += weights_(i)*x_diff*z_diff.transpose();
	}

	//calculate Kalman gain K;
	MatrixXd K = Tc*S.inverse();

	//update state mean and covariance matrix
	VectorXd z_diff = meas_package.raw_measurements_ - z_pred;

	x_ = x_ + K * z_diff;

	P_ = P_ - K * S * K.transpose();

	//NIS Update
	NIS_laser_ = z_diff.transpose() * S.inverse() * z_diff;

	std::cout<<"UpdateLidar End"<<std::endl;
}

/**
 * Updates the state and the state covariance matrix using a radar measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateRadar(MeasurementPackage meas_package) {
	std::cout<<"UpdateRadar start"<<std::endl;
	const int n_z = 3;
	//create matrix for sigma points in measurement space
	MatrixXd Zsig = MatrixXd(n_z, 2 * n_aug_ + 1);
	Zsig.fill(0.0);

	//mean predicted measurement
	VectorXd z_pred = VectorXd(n_z);
	z_pred.fill(0.0);

	//measurement covariance matrix S
	MatrixXd S = MatrixXd(n_z,n_z);
	S.fill(0.0);

	//transform sigma points into measurement space
	for(int i=0; i< 2 * n_aug_ + 1; ++i){
		VectorXd tmp = Xsig_pred_.col(i);
		float px = tmp(0);
		float py = tmp(1);
		float v = tmp(2);
		float psi = tmp(3);
		float psi_dot = tmp(4);

		//check for zeros
		if (fabs(px) < 0.001) {
			px = 0.001;
		}
		if (fabs(py) < 0.001) {
			py = 0.001;
		}

		Zsig.col(i) << sqrt(px*px+py*py),
				atan2(py, px),
				((px*cos(psi)*v+py*sin(psi)*v)/(sqrt(px*px+py*py)));
	}

	//calculate mean predicted measurement
	for(int i=0; i< 2 * n_aug_ + 1; ++i){
		z_pred += weights_(i) * Zsig.col(i);
	}

	//std::cout<<"Calculating S"<<std::endl;
	//calculate measurement covariance matrix S
	for(int i=0; i< 2 * n_aug_ + 1; ++i){
		VectorXd tmp = Zsig.col(i) - z_pred;

		//angle normalization
		while (tmp(1)> M_PI) tmp(1)-=2.*M_PI;
		while (tmp(1)<-M_PI) tmp(1)+=2.*M_PI;

		S = S + weights_(i)*tmp*tmp.transpose();
	}

	//add measurement noise covariance matrix
	S = S+R_radar_;


	//create matrix for cross correlation Tc
	MatrixXd Tc = MatrixXd(n_x_, n_z);
	Tc.fill(0.0);

	//calculate cross correlation matrix
	for(int i=0; i<(2 * n_aug_ + 1); ++i){
		VectorXd x_diff = Xsig_pred_.col(i) - x_;
		while(x_diff(3) > M_PI) x_diff(3) -= 2.*M_PI;
		while(x_diff(3) < -M_PI) x_diff(3) += 2.*M_PI;

		VectorXd z_diff = Zsig.col(i) - z_pred;
		while(z_diff(1) > M_PI) z_diff(1) -= 2.*M_PI;
		while(z_diff(1) < -M_PI) z_diff(1) += 2.*M_PI;
		Tc += weights_(i)*x_diff*z_diff.transpose();
	}

	//calculate Kalman gain K;
	MatrixXd K = Tc*S.inverse();

	//update state mean and covariance matrix
	VectorXd z_diff = meas_package.raw_measurements_ - z_pred;
	while(z_diff(1) > M_PI) z_diff(1)-=2.*M_PI;
	while(z_diff(1) < -M_PI) z_diff(1)+=2.*M_PI;
	x_ = x_ + K * z_diff;

	P_ = P_ - K * S * K.transpose();

	//NIS Update
	NIS_radar_ = z_diff.transpose() * S.inverse() * z_diff;

	std::cout<<"UpdateRadar End"<<std::endl;
}

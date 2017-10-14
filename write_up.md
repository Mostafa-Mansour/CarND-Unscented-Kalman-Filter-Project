# Code organization
---

[image1]: ./images/data.png "data"

[image1]: ./data.png "data"

The data file can be checked [here](obj_pose-laser-radar-synthetic-input.txt).

* ![alt text][image1]
* Each row represents a sensor measurement where the first column tells you if the measurement comes from radar (R) or lidar (L).

* For a row containing radar data, the columns are: sensor_type, rho_measured, phi_measured, rhodot_measured, timestamp, x_groundtruth, y_groundtruth, vx_groundtruth, vy_groundtruth, yaw_groundtruth, yawrate_groundtruth.

* For a row containing lidar data, the columns are: sensor_type, x_measured, y_measured, timestamp, x_groundtruth, y_groundtruth, vx_groundtruth, vy_groundtruth, yaw_groundtruth, yawrate_groundtruth.

* Whereas radar has three measurements (rho, phi, rhodot), lidar has two measurements (x, y).

* I will use the measurement values and timestamp in Kalman filter algorithm. Groundtruth, which represents the actual path the bicycle took, is for calculating root mean square error.

---
The project consists of Three main classes:
* UKF class
* Tools class
*Measurement
---
Measurement class
* contains three private members, timestamp, enumeration of sensor type and raw measurements.
---
UKF class
* contains the parameters and the functions of Unscented Kalman filter.
---
Tools class

* contains auxiliary function to calculate RMSE.

___
* By using the data set in the text file (data set 1), UKF converges with RMSE (0.069,0.0605,0.17,0.175) after 500 measurements from both a Lidar and a RADAR.

* With the same data set, if Lidar measuements only are used, UKF will also converges but with RMSE higher than the previous case, RMSE is about (0.08,0.084,0.29,0.174).

* With the same data set, if Radar measuements only are used, UKF will also converges but with RMSE higher than in the Lidar case, RMSE is about (0.11,0.11,0.18,0.22).

* One can claim that using Radar only can give worse RMSE w.r.t. x and y position that Lidar. It is true because using a Lidar allow us to measure x and y directly which is not the case in Radar. However, Radar gives RMSE w.r.t. velocity better than Lidar and this is because velocity and bearing angle is directly measured using a Radar.





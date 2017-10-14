#include <iostream>
#include <fstream>
#include <vector>
#include "measurement_package.h"
#include "ukf.h"
#include "tools.h"
int main() {
    std::cout << "UKF implementation using RADAR & LIDAR signals for pedestrian detection" << std::endl;

    std::vector<MeasurementPackage> measurementList;
    std::vector<VectorXd> ground_truthList;
    std::vector<VectorXd> estimations;
    MeasurementPackage measurement;
    std::string filePath="../src/data/obj_pose-laser-radar-synthetic-input.txt";
    std::ifstream fileIn(filePath,std::ios::app);
    std::string line;
    std::string sensor;
    double x,y,x_true,y_true,vx_true,vy_true;
    double rho,phi,rhodot;
    long timeStamp;

    VectorXd ground_truth(4);

    if(fileIn.is_open()){
        std::cout<<"file is opened"<<std::endl;
        std::istringstream iss;
        while(std::getline(fileIn,line)){
            iss.str(line);
            iss>>sensor;

            // Check for sensor type
            //std::cout<<"Sensor compare";
            if(sensor.compare("L")==0){

                measurement.sensor_type_=MeasurementPackage::LASER;
                measurement.raw_measurements_=VectorXd(2);
                iss>>x;
                iss>>y;
                measurement.raw_measurements_<<x,y;
                std::cout<<measurement.raw_measurements_<<std::endl;
                iss>>timeStamp;
                measurement.timestamp_=timeStamp;
                iss>>x_true;
                iss>>y_true;
                iss>>vx_true;
                iss>>vy_true;

                ground_truth<<x_true,y_true,vx_true,vy_true;
            }
            else if(sensor.compare("R")==0){
                measurement.sensor_type_=MeasurementPackage::RADAR;
                measurement.raw_measurements_=VectorXd(3);
                iss>>rho;
                iss>>phi;
                iss>>rhodot;
                measurement.raw_measurements_<<rho,phi,rhodot;
                iss>>timeStamp;
                measurement.timestamp_=timeStamp;
                iss>>x_true;
                iss>>y_true;
                iss>>vx_true;
                iss>>vy_true;

                ground_truth<<x_true,y_true,vx_true,vy_true;
            }
            else{
                std::cout<<"Unknown sensor type"<<std::endl;
                return -1;
            }
            measurementList.push_back(measurement);
            ground_truthList.push_back(ground_truth);

        }

    }
    fileIn.close();
    UKF track;
    Tools mytool;

    for(int i=0;i<measurementList.size();i++){
        //std::cout<<"Test"<<std::endl;
        track.ProcessMeasurement(measurementList[i]);
        //std::cout<<"hi";
        //std::cout<<track.x_<<std::endl;
        double p_x=track.x_(0);
        double p_y=track.x_(1);
        double v=track.x_(2);
        double phi=track.x_(3);
        double v_x=v*cos(phi);
        double v_y=v*sin(phi);
        VectorXd estimation(4);
        estimation<<p_x,p_y,v_x,v_y;
        estimations.push_back(estimation);


    }
    //std::cout<<"Estimations are "<<estimations[0]<<std::endl;
    //std::cout<<"Ground truth is "<<ground_truthList.size()<<std::endl;

    mytool.CalculateRMSE(estimations,ground_truthList);

    return 0;
}
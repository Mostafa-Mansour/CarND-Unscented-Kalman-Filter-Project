#include <iostream>
#include "tools.h"

using Eigen::VectorXd;
using Eigen::MatrixXd;
using std::vector;

Tools::Tools() {}

Tools::~Tools() {}

VectorXd Tools::CalculateRMSE(const vector<VectorXd> &estimations,
                              const vector<VectorXd> &ground_truth) {
  /**
  TODO:
    * Calculate the RMSE here.
  */
    std::cout<<"testing";
    VectorXd rmse=VectorXd(4);
    rmse<<0,0,0,0;
    //std::cout<<"Entered"<<std::endl;
    if(estimations.size()==0 || estimations.size()!=ground_truth.size()){
        std::cout<<"Check your vectors"<<std::endl;
        return rmse;
    }
    VectorXd residual(4);
    residual.fill(0.0);
    for(unsigned int i=0;i<estimations.size();i++){
        residual=estimations[i]-ground_truth[i];
        residual=residual.array()*residual.array();
        rmse+=residual;
    }

    rmse=rmse/estimations.size();
    rmse=rmse.array().sqrt();
    std::cout<<"RMSE is "<<rmse<<std::endl;

    return rmse;
}

VectorXd Tools::Nis(VectorXd z, VectorXd z_pred, MatrixXd S) {
    return (z-z_pred).transpose()*S.inverse()*(z-z_pred);
}
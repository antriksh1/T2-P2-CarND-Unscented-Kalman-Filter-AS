#include <iostream>

#include "tools.h"

using Eigen::VectorXd;
using std::vector;

using namespace std;

Tools::Tools() {}

Tools::~Tools() {}

VectorXd Tools::CalculateRMSE(const vector<VectorXd> &estimations,
                              const vector<VectorXd> &ground_truth) {
  /**
   * TODO: Calculate the RMSE here.
   */
 	VectorXd rmse(4);
  rmse << 0,0,0,0;

  // check the validity of the following inputs:
  //  * the estimation vector size should not be zero
  //  * the estimation vector size should equal ground truth vector size
  // ... your code here
  if(estimations.size() == 0) {
      cout << "The estimation vector size should not be zero" << endl;
      return rmse;
  }

  if(estimations.size() != ground_truth.size()) {
      cout << "The estimation vector size NOT EQUAL TO ground truth vector size" << endl;
      return rmse;
  }

  //accumulate squared residuals
  for(int i=0; i < estimations.size(); ++i){
      for(int j = 0; j < rmse.size(); j++) {
          float diff = estimations[i](j) - ground_truth[i](j);
          rmse(j) = rmse(j) + diff*diff;
      }
  }

  //calculate the mean
  // ... your code here
  for(int j = 0; j < rmse.size(); j++) {
      rmse(j) = rmse(j)/estimations.size();
  }

  //calculate the squared root
  // ... your code here
  for(int j = 0; j < rmse.size(); j++) {
    rmse(j) = sqrt(rmse(j));
  }

  //return the result
  return rmse;
}
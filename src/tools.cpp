#include <iostream>
#include "tools.h"

#define SMALL_NUM 0.0001f
#define SMALL_NUM_SQ 0.0000001f

using Eigen::VectorXd;
using Eigen::MatrixXd;
using std::vector;

Tools::Tools() {}

Tools::~Tools() {}

VectorXd Tools::CalculateRMSE(const vector<VectorXd> &estimations,
                              const vector<VectorXd> &ground_truth) {
	VectorXd rmse(4);
	rmse << 0,0,0,0;

	if (estimations.size() == 0 || estimations.size() != ground_truth.size()) {
		cout << "CalculateRMSE(): Div by 0 error" << endl;
		return rmse;
	}

	// accumulate squared residuals
	for (int i = 0; i < estimations.size(); ++i) {
		VectorXd residual = estimations[i] - ground_truth[i];
		residual = residual.array() * residual.array();
		rmse += residual;
	}

	// calculate the mean
	rmse = rmse / estimations.size();

	// calculate the sqrt
	rmse = rmse.array().sqrt();

	return rmse;
}

MatrixXd Tools::CalculateJacobian(const VectorXd& x_state) {

	MatrixXd Hj(3,4);

	float px = x_state(0);
	float py = x_state(1);
	if (fabs(px) < SMALL_NUM && fabs(py) < SMALL_NUM) {
		px = SMALL_NUM;
		py = SMALL_NUM;
	}
	float vx = x_state(2);
	float vy = x_state(3);

	float c1 = px*px + py*py;
	if (c1 < SMALL_NUM_SQ)
		c1 = SMALL_NUM_SQ;

	float c1sqrt = sqrt(c1);
	float c132 = pow(c1, 3/2.0);

	if (fabs(c1) < 0.0001) {
		cout << "CalculateJacobian(): Div by 0 error" << endl;
	} else {
		Hj << px / c1sqrt, py / c1sqrt, 0.0f, 0.0f,
			-py / c1, px / c1, 0.0f, 0.0f,
			py*(vx*py-vy*px) / c132, px*(vy*px-vx*py) / c132, px / c1sqrt, py / c1sqrt;
	}
	return Hj;
}

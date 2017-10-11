#include "PID.h"

using namespace std;

PID::PID() 
{
	p_error = 0.0;
	i_error = 0.0;
	d_error = 0.0;
}

PID::~PID() {}

void PID::Init(double Kp_in, double Ki_in, double Kd_in) 
{
	Kp = Kp_in;
	Ki = Ki_in;
	Kd = Kd_in;
	std::cout << "Kp = "<<Kp<<" Ki = "<<Ki<<" Kd = "<<Kd<<std::endl;
}

void PID::UpdateError(double cte) 
{
	if(p_error == 0.0){//For the first CTE, init p_error
		p_error = cte;
	}
	d_error = cte - p_error;
    p_error = cte;
    i_error += cte;
}

double PID::getSteerVal() 
{
	double steerVal = -1 * (Kp * p_error + Ki * i_error + Kd * d_error);
	if(steerVal > 1){
		steerVal = 1;
	}
	if(steerVal < -1){
		steerVal = -1;
	}
	return steerVal;
}


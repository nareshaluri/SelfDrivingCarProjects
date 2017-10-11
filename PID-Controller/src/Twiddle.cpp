#include "Twiddle.h"
#include <cmath>
using namespace std;

Twiddle::Twiddle()
{
	total_error = 0.0;
	cur_param_update_idx = 2;
	is_after_first_update = false;
	best_error = std::numeric_limits<double>::max();
	std::cout << "total_error = "<<total_error<<" best_error = "<<best_error<<std::endl;
	increment_dp = decrement_dp = false;
}

Twiddle::~Twiddle()
{

}

void Twiddle::init(double kp, double ki, double kd){
	K = {kp, ki, kd};
	dp = {0.1*kp, 0.1*ki, 0.1*kd};
}

void Twiddle::calculateTotalError(double cte)
{
	total_error += pow(cte,2);
}

void Twiddle::updateParameters()
{	
	std::cout << "Twiddle::updateParameters"<<std::endl;
	std::cout << "total_error = "<<total_error<<" best_error = "<<best_error<<std::endl;
	if (total_error < best_error) 
	{
        cout << "improvement!" << endl;
        best_error = total_error;
        if (is_after_first_update) 
        {
            // Skip initial update, as this is always an improvement because initial best_error
            // is always greater than total_error
            dp[cur_param_update_idx] *= 1.1;            
        }
        // next parameter
        cur_param_update_idx = (cur_param_update_idx + 1) % 3;
        increment_dp = decrement_dp = false;
    }
    if (!increment_dp && !decrement_dp) 
    {
        // try adding dp[i] to params[i]
        K[cur_param_update_idx] += dp[cur_param_update_idx];
        increment_dp = true;
    }
    else if (increment_dp && !decrement_dp) 
    {
        // try subtracting dp[i] from params[i]
        K[cur_param_update_idx] += -2 * dp[cur_param_update_idx];
        decrement_dp = true;         
    }
    else 
    {
        // set it back, reduce dp[i], move on to next parameter
        K[cur_param_update_idx] += dp[cur_param_update_idx];  
        dp[cur_param_update_idx] *= 0.9;

        // next parameter
        cur_param_update_idx = (cur_param_update_idx + 1) % 3;
        increment_dp = decrement_dp = false;
    }
    total_error = 0.0;//reset to use in next cycle
    is_after_first_update = true;
    cout << "new parameters" << endl;
    cout << "P: " << K[0] << ", I: " << K[1] << ", D: " << K[2] << endl; 
}

std::vector<double> Twiddle::getUpdatedParams()
{
	return K;
}
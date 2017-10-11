#ifndef TWIDDLE_H
#define TWIDDLE_H
#include <iostream>
#include <vector>
class Twiddle{
public:
  /*
  * Constructor
  */
  Twiddle();

  /*
  * Destructor.
  */
  virtual ~Twiddle();

  void init(double kp, double ki, double kd);
  void calculateTotalError(double cte);
  void updateParameters();
  std::vector<double> getUpdatedParams();

private:
	double total_error;
	double best_error;
	bool is_after_first_update;
	bool increment_dp;
	bool decrement_dp;
	int cur_param_update_idx;
	std::vector<double> K;
	std::vector<double> dp;
};
#endif /* TWIDDLE_H */
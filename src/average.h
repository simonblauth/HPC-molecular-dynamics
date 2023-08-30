#ifndef __AVERAGE_H_
#define __AVERAGE_H_

#include "stddef.h"

// Holds a cumulative or exponential average over a value
class Average {
  protected:
    double value_;
    Average() : value_(0) {}
    Average(double init) : value_(init) {}

  public:
    double get() { return value_; }
};

// Holds a cumulative average over a value
class CumulativeAverage: public Average {
  private:
    size_t period_;

  public:
    CumulativeAverage(size_t period) : Average(), period_(period) {}
    CumulativeAverage(size_t period, double init) : Average(init), period_(period) {}
    // update for cumulative average
    void update(double new_val, size_t timestep) {
        size_t step = timestep % period_ + 1;
        if (step == 1) {
            value_ = 0;
        }
        value_ += (new_val - value_) / step;
    }
};

// Holds an exponential average over a value
class ExponentialAverage: public Average {
  private:
    double alpha_;

  public:
    ExponentialAverage(double alpha) : Average(), alpha_(alpha) {}
    ExponentialAverage(double alpha, double init) : Average(init), alpha_(alpha) {}
    // update for exponential average
    void update(double new_val) {
        value_ *= (1 - alpha_);
        value_ += alpha_ * new_val;
    }
};

#endif // __AVERAGE_H_

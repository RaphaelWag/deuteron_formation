//
// Created by rwagn on 20.03.2019.
//

#ifndef MYSIMULATION_HISTOGRAM_H
#define MYSIMULATION_HISTOGRAM_H

#include <iostream>
#include <fstream>
#include <vector>
#include <chrono>
#include <random>

using namespace std;

class histogram {
private:

    double *error_;
    double *boundaries_;
    double *rescale_factors_;

    int N_bins_;
    int size_;
    bool normalize_bin_;
    bool rescale_;
    double rescale_factor_;
    bool rescale_array_;
    bool rescale_probability_;
    bool is_scaled_;

    void set_boundaries();

    void generate_arrays(int N_bins);

    void evaluate_error();

    void init(int N_bins);

    void populate(const histogram &other);


public:

    double *x_axis_;
    double *data_;

    histogram();//default constructor
    histogram(const histogram &other); // copy constructor
    ~histogram(); //destructor
    void set_bins(int N);

    int N_bins();

    void set_xaxis(const double *x_axis, int N);

    void set_xaxis(vector<double> x_axis);

    void set_boundaries(vector<double> boundaries);

    void fill(double value, double weight = 1.);

    void set_data(vector<double> &source);

    void set_error(vector<double> &source);

    void rescale();

    void not_scale();

    void print(const string &filename);

    void rescale_data(double value);

    void rescale_data(vector<double> values);

    void normalize(); //normalizes data by bin width

    void probability();

    void reset();

    void reset_scaling();

    double chi_squared_raw(histogram &data_hist);

    histogram &operator=(const histogram &other); // copy assignment operator
    histogram &operator+(const histogram &other); // combine histograms for more mc loops

};

#endif //MYSIMULATION_HISTOGRAM_H

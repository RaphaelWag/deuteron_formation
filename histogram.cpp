//
// Created by rwagn on 20.03.2019.
//

#include "histogram.h"
#include <cmath>
#include <iostream>
#include <cstring>
#include <vector>

histogram::histogram() {
    //generate arrays
    init(0);
}

void histogram::set_bins(int N) {
    N_bins_ = N;
    generate_arrays(N_bins_);
}

void histogram::set_xaxis(const double *x_axis_init, int N) {
    N_bins_ = N;
    generate_arrays(N_bins_);

    for (int i = 0; i < N_bins_; ++i) {
        x_axis_[i] = x_axis_init[i];
    }
    set_boundaries();

}

void histogram::set_xaxis(vector<double> x_axis) {
    N_bins_ = x_axis.size();
    generate_arrays(N_bins_);

    for (int i = 0; i < N_bins_; ++i) {
        x_axis_[i] = x_axis[i];
    }
    set_boundaries();
}

void histogram::set_boundaries() {

    boundaries_[1] = 0.5 * (x_axis_[0] + x_axis_[1]);
    boundaries_[0] = x_axis_[0] - (boundaries_[1] - x_axis_[0]);

    for (int i = 2; i < N_bins_ + 1; ++i) {
        boundaries_[i] = x_axis_[i - 1] + (x_axis_[i - 1] - boundaries_[i - 1]);
    }

}

void histogram::fill(double value, double weight) {
    for (int i = 0; i < N_bins_; ++i) {
        if (value < boundaries_[i + 1]) {
            if (value > boundaries_[i]) {
                data_[i] += weight;
                error_[i] += weight * weight;
            }
        }
    }
}

void histogram::print(const string &filename) {
    //rescaling
    rescale();

    //printing
    ofstream myfile(filename);
    if (!myfile.is_open()) {
        cout << "Error in: histogram -> print" << endl;
        cout << "Unable to open file" << endl;
    } else {
        for (int i = 0; i < N_bins_; ++i) {
            myfile << x_axis_[i] << " " << data_[i] << " " << error_[i] << endl;
        }
        myfile.close();
    }
}

void histogram::generate_arrays(int N_bins) {

    delete[] x_axis_;
    delete[] boundaries_;
    delete[] rescale_factors_;
    delete[] error_;
    delete[] data_;

    x_axis_ = new double[N_bins];
    boundaries_ = new double[N_bins + 1];
    rescale_factors_ = new double[N_bins];
    error_ = new double[N_bins];
    data_ = new double[N_bins];
    for (int i = 0; i < N_bins; ++i) {
        data_[i] = 0;
        error_[i] = 0;
    }
}

void histogram::rescale_data(double value) {
    rescale_ = true;
    rescale_factor_ = value;
    is_scaled_ = false;
}

void histogram::normalize() {
    normalize_bin_ = true;
    is_scaled_ = false;
}

void histogram::set_boundaries(vector<double> boundaries) {
    N_bins_ = boundaries.size() - 1;
    generate_arrays(N_bins_);

    for (int i = 0; i < boundaries.size(); ++i) {
        boundaries_[i] = boundaries[i];
    }
    for (int j = 0; j < N_bins_; ++j) {
        x_axis_[j] = 0.5 * (boundaries_[j + 1] + boundaries_[j]);
    }
}

void histogram::rescale() {

    if (!is_scaled_) {

        //evaluate error
        evaluate_error();

        //noramlization per bin
        if (normalize_bin_) {
            for (int i = 0; i < N_bins_; ++i) {
                data_[i] = data_[i] / (boundaries_[i + 1] - boundaries_[i]);
                error_[i] = error_[i] / (boundaries_[i + 1] - boundaries_[i]);
            }
        }
        //rescaling each bin
        if (rescale_array_) {
            for (int i = 0; i < N_bins_; ++i) {
                data_[i] = data_[i] * rescale_factors_[i];
                error_[i] = error_[i] * rescale_factors_[i];
            }
        }
        //rescaling all bins
        if (rescale_) {
            for (int i = 0; i < N_bins_; ++i) {
                data_[i] = data_[i] * rescale_factor_;
                error_[i] = error_[i] * rescale_factor_;
            }
        }
        //normalize histogram to 1
        if (rescale_probability_) {
            double S = 0;
            for (int i = 0; i < N_bins_; ++i) {
                S += data_[i] * (boundaries_[i + 1] - boundaries_[i]);
            }
            for (int j = 0; j < N_bins_; ++j) {
                data_[j] = data_[j] / S;
                error_[j] = error_[j] / S;
            }
        }

        is_scaled_ = true;
    }
}

void histogram::rescale_data(vector<double> values) {
    rescale_array_ = true;
    is_scaled_ = false;
    for (int i = 0; i < N_bins_; ++i) {
        rescale_factors_[i] = values[i];
    }
}

void histogram::probability() {
    rescale_probability_ = true;
    is_scaled_ = false;
}

void histogram::evaluate_error() {
    for (int i = 0; i < N_bins_; ++i) {
        error_[i] = sqrt(error_[i]);
    }
}

histogram::histogram(const histogram &other) {
    init(other.N_bins_);
    populate(other);
}

histogram::~histogram() {

    delete[] x_axis_;
    delete[] data_;
    delete[] error_;
    delete[] boundaries_;
    delete[] rescale_factors_;

}

histogram &histogram::operator=(const histogram &other) {

    if (this != &other) {
        init(other.N_bins_);
        populate(other);
    }
    return *this;
}

void histogram::init(int N_bins) {

    x_axis_ = new double[N_bins];
    boundaries_ = new double[N_bins + 1];
    rescale_factors_ = new double[N_bins];
    error_ = new double[N_bins];
    data_ = new double[N_bins];
    for (int i = 0; i < N_bins; ++i) {
        data_[i] = 0;
        error_[i] = 0;
    }

    normalize_bin_ = false;
    rescale_ = false;
    rescale_factor_ = 1;
    rescale_array_ = false;
    rescale_probability_ = false;
    is_scaled_ = false;
    N_bins_ = N_bins;
    size_ = N_bins_ * sizeof(double);

}

void histogram::populate(const histogram &other) {

    memcpy(x_axis_, other.x_axis_, size_);
    memcpy(boundaries_, other.boundaries_, size_ + sizeof(double));
    memcpy(rescale_factors_, other.rescale_factors_, size_);
    memcpy(error_, other.error_, size_);
    memcpy(data_, other.data_, size_);

    normalize_bin_ = other.normalize_bin_;
    rescale_ = other.rescale_;
    rescale_factor_ = other.rescale_factor_;
    rescale_array_ = other.rescale_array_;
    rescale_probability_ = other.rescale_probability_;
    N_bins_ = other.N_bins_;
    size_ = other.N_bins_ * sizeof(double);

}

void histogram::set_data(vector<double> &source) {
    copy(source.begin(), source.end(), data_);
}

void histogram::set_error(vector<double> &source) {
    copy(source.begin(), source.end(), error_);
}

int histogram::N_bins() {
    return N_bins_;
}

void histogram::reset() {

    for (int i = 0; i < N_bins_; ++i) {
        data_[i] = 0;
        error_[i] = 0;
    }
    is_scaled_ = false;

}

void histogram::not_scale() {
    is_scaled_ = true;
}

double histogram::chi_squared_raw(histogram &data_hist) {

    double chi_squared = 0;
    double error = 0;
    double value;

    if (data_hist.N_bins() == N_bins_) {
        for (int i = 0; i < data_hist.N_bins_; ++i) {
            error = data_hist.error_[i] * data_hist.error_[i];
            value = data_hist.data_[i] - data_[i];
            value = value * value;
            chi_squared += value / error;
        }

    } else {
        cout << "Histograms do not have the same amount of bins" << endl;
    }
    return chi_squared;
}

void histogram::reset_scaling() {

}

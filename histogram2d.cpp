//
// Created by rwagn on 06.08.2019.
//

#include "histogram2d.h"
#include <cstring>
#include <iostream>
#include <fstream>

using namespace std;

histogram2d::histogram2d() {
    init(0, 0);
}

histogram2d::~histogram2d() {
    delete[] boundaries_x_;
    delete[] boundaries_y_;
    delete[] bin_centre_x_;
    delete[] bin_centre_y_;
    for (int i = 0; i < N_bins_x_; ++i) {
        delete[] data_[i];
    }
    delete[] data_;
}

histogram2d::histogram2d(const histogram2d &other) {
    init(other.N_bins_x_, other.N_bins_y_);
    populate(other);
}

histogram2d &histogram2d::operator=(const histogram2d &other) {
    if (this != &other) {
        init(other.N_bins_x_, other.N_bins_y_);
        populate(other);
    }
    return *this;
}

void histogram2d::set_bins(int bins_x, int bins_y) {
    N_bins_x_ = bins_x;
    N_bins_y_ = bins_y;
    init(N_bins_x_, N_bins_y_);
}

void histogram2d::set_boundaries(double x_min, double x_max, double y_min, double y_max) {
    boundaries_x_[0] = x_min;
    boundaries_y_[0] = y_min;
    boundaries_x_[N_bins_x_] = x_max;
    boundaries_y_[N_bins_y_] = y_max;
}

void histogram2d::fill(double value_x, double value_y, double weight_x, double weight_y) {
//TODO check this
    for (int i = 0; i < N_bins_x_; ++i) {
        if (value_x < boundaries_x_[i]) { continue; }
        if (value_x >= boundaries_x_[i + 1]) { continue; }
        for (int j = 0; j < N_bins_y_; ++j) {
            if (value_y < boundaries_y_[j]) { continue; }
            if (value_y >= boundaries_y_[j + 1]) { continue; }
            data_[i][j] += weight_x * weight_y;
            break;
        }
    }
}

void histogram2d::init(int N_bins_x, int N_bins_y) {

    N_bins_x_ = N_bins_x;
    N_bins_y_ = N_bins_y;
    x_log_ = y_log_ = probability_ = rescaled_ = false;

    data_ = new double *[N_bins_x];
    for (int i = 0; i < N_bins_x; ++i) {
        data_[i] = new double[N_bins_y];
    }

    for (int j = 0; j < N_bins_x; ++j) {
        for (int i = 0; i < N_bins_y; ++i) {
            data_[j][i] = 0;
        }
    }

    boundaries_x_ = new double[N_bins_x + 1];
    boundaries_y_ = new double[N_bins_y + 1];
    bin_centre_x_ = new double[N_bins_x];
    bin_centre_y_ = new double[N_bins_y];
}

void histogram2d::populate(const histogram2d &other) {

    N_bins_x_ = other.N_bins_x_;
    N_bins_y_ = other.N_bins_y_;
    x_log_ = other.x_log_;
    y_log_ = other.y_log_;
    probability_ = other.probability_;
    rescaled_ = other.rescaled_;

    double size_y = N_bins_y_ * sizeof(double);
    double b_size_x = (N_bins_x_ + 1) * sizeof(double);
    double b_size_y = (N_bins_y_ + 1) * sizeof(double);
    double bin_size_x = N_bins_x_ * sizeof(double);
    double bin_size_y = N_bins_y_ * sizeof(double);

    for (int i = 0; i < N_bins_x_; ++i) {
        memcpy(data_[i], other.data_[i], size_y);
    }
    memcpy(boundaries_x_, other.boundaries_x_, b_size_x);
    memcpy(boundaries_y_, other.boundaries_y_, b_size_y);
    memcpy(bin_centre_x_, other.bin_centre_x_, bin_size_x);
    memcpy(bin_centre_y_, other.bin_centre_y_, bin_size_y);
}

void histogram2d::set_axis_lin(double *boundaries, double *bin_centre, int N_bins) {
    double stepsize = (boundaries[N_bins] - boundaries[0]) / N_bins;
    for (int i = 1; i < N_bins; ++i) {
        boundaries[i] = boundaries[i - 1] + stepsize;
    }
    for (int j = 0; j < N_bins; ++j) {
        bin_centre[j] = 0.5 * (boundaries[j + 1] + boundaries[j]);
    }
}

void histogram2d::set_x_lin() {
    set_axis_lin(boundaries_x_, bin_centre_x_, N_bins_x_);
    x_log_ = false;
}

void histogram2d::set_y_lin() {
    set_axis_lin(boundaries_y_, bin_centre_y_, N_bins_y_);
    y_log_ = false;
}

void histogram2d::set_axis_log(double *boundaries, double *bin_centre, int N_bins) {
    for (int i = 0; i < N_bins; ++i) {
        boundaries[i + 1] = boundaries[i] * 10;
        bin_centre[i] = 0;
    }
}

void histogram2d::set_x_log() {
    set_axis_log(boundaries_x_, bin_centre_x_, N_bins_x_);
    x_log_ = true;
}

void histogram2d::set_y_log() {
    set_axis_log(boundaries_y_, bin_centre_y_, N_bins_y_);
    y_log_ = true;
}

void histogram2d::reset_data() {
    rescaled_=false;
    for (int i = 0; i < N_bins_x_; ++i) {
        for (int j = 0; j < N_bins_y_; ++j) {
            data_[i][j] = 0;
        }
    }
}

void histogram2d::print(const string &filename) {
    if (!rescaled_) { rescale(); }
    ofstream myfile(filename);
    if (!myfile.is_open()) {
        cout << "Error in: histogram2d -> print" << endl;
        cout << "Unable to open file" << endl;
    } else {
        if (x_log_) {
            for (int k = 0; k < N_bins_x_; ++k) {
                myfile << boundaries_x_[k] << " ";
            }
        } else {
            for (int i = 0; i < N_bins_x_; ++i) {
                myfile << bin_centre_x_[i] << " ";
            }
        }
        myfile << "\n";

        if (y_log_) {
            for (int k = 0; k < N_bins_y_; ++k) {
                myfile << boundaries_y_[k] << " ";
            }
        } else {
            for (int i = 0; i < N_bins_y_; ++i) {
                myfile << bin_centre_y_[i] << " ";
            }
        }
        myfile << endl;

        for (int i = 0; i < N_bins_x_; ++i) {
            for (int j = 0; j < N_bins_y_; ++j) {
                myfile << data_[i][j] << " ";
            }
            myfile << "\n";
        }
        myfile.close();
    }
}

void histogram2d::probability() {
    probability_ = true;
}

void histogram2d::rescale() {

    if (!rescaled_) {

        if (probability_) {
            double norm;
            norm =0;
            //calculate total area
            double len_x = boundaries_x_[N_bins_x_] - boundaries_x_[0];
            double len_y = boundaries_y_[N_bins_y_] - boundaries_y_[0];
            double area = len_x * len_y;

            //calculate total events
            for (int i = 0; i < N_bins_x_; ++i) {
                len_x = boundaries_x_[i + 1] - boundaries_x_[i];
                for (int j = 0; j < N_bins_y_; ++j) {
                    len_y = boundaries_y_[j + 1] - boundaries_y_[j];
                    norm += data_[i][j]*len_x*len_y;
                }
            }

            norm = 1. / norm;

            //rescale to probability distribution
            for (int k = 0; k < N_bins_x_; ++k) {
                len_x = boundaries_x_[k + 1] - boundaries_x_[k];
                for (int i = 0; i < N_bins_y_; ++i) {
                    len_y = boundaries_y_[i + 1] - boundaries_y_[i];
                    data_[k][i] = data_[k][i] * norm * len_x * len_y;
                }
            }
            double check = 0;
            for (int l = 0; l < N_bins_x_; ++l) {
                for (int i = 0; i < N_bins_y_; ++i) {
                    check += data_[l][i];
                }
            }
            if((check-1.)>0.00000001){
                cout << "probability normalization failed" << endl;
            }
        }
    }
    rescaled_ = true;
}

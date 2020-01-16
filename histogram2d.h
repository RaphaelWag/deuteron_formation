//
// Created by rwagn on 06.08.2019.
//

#ifndef MYSIMULATION_HISTOGRAM2D_H
#define MYSIMULATION_HISTOGRAM2D_H

#include <iostream>
using namespace std;

class histogram2d {

private:
    double **data_;
    double *boundaries_x_;
    double *boundaries_y_;
    double *bin_centre_x_;
    double *bin_centre_y_;

    int N_bins_x_;
    int N_bins_y_;
    bool x_log_;
    bool y_log_;
    bool probability_;
    bool rescaled_;


    void init(int N_bins_x, int N_bins_y);

    void populate(const histogram2d &other);

    static void set_axis_lin(double *boundaries, double *bin_centre, int N_bins);

    static void set_axis_log(double *boundaries, double *bin_centre, int N_bins);


public:
    histogram2d(); //default constructor
    ~histogram2d(); //destructor
    histogram2d(const histogram2d &other); //copy constructor
    histogram2d &operator=(const histogram2d &other); // copy assignment operator

    void set_bins(int bins_x, int bins_y);

    void set_boundaries(double x_min, double x_max, double y_min, double y_max);

    void fill(double value_x, double value_y, double weight_x = 1., double weight_y = 1.);

    void set_x_lin();

    void set_y_lin();

    void set_x_log();

    void set_y_log();

    void reset_data();

    void print(const string& filename);

    void probability();

    void rescale();
};


#endif //MYSIMULATION_HISTOGRAM2D_H

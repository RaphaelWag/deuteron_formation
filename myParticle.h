//
// Created by rwagn on 07.03.2019.
//

#ifndef MYSIMULATION_MYPARTICLE_H
#define MYSIMULATION_MYPARTICLE_H

#include <iostream>
#include <vector>

using namespace std;

class myParticle {

private:
    int id_; //particle id
    vector<double> p_; //4 momentum
    vector<double> x_; //4 vertex
    double wA_; //event weight for ALICE weighting
    double wLT_; //event weight for Levy Tsallis weighting
    double Nch_eta_; //number of charges particles per rapidity
    bool used_;
    double y_n_;
    double y_p_;

    void init();

public:

    myParticle();

    explicit myParticle(string *tokens);

    int id();

    double p(int i);

    void set_p(double E, double px, double py, double pz);

    double y();

    double x(int i);

    void set_x(double xt, double xx, double xy, double xz);

    double pT();

    double wA();
    double wLT();

    void set_wA(double w_init);
    void set_wLT(double w_init);

    bool is_used();

    void used();

    void set_parent_y(double y_n, double y_p);

    double y_p();

    double y_n();

    void set_Nch_eta(double Nch_eta);

    double Nch_eta();

    double momentum_difference(myParticle &other);

    double spacetime_difference(myParticle &other);

    void bstback(const myParticle& boost_particle);

    //overload + for deuteron = neutron + proton
    myParticle operator+(myParticle &proton);

};

#endif //MYSIMULATION_MYPARTICLE_H

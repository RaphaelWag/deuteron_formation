//
// Created by rwagn on 07.03.2019.
//

#include "myParticle.h"
#include <cmath>

int myParticle::id() {
    return id_;
}

double myParticle::p(int i) {
    return p_[i];
}

double myParticle::x(int i) {
    return x_[i];
}

myParticle::myParticle() {
    init();
}

myParticle::myParticle(string *tokens) {

    init();
    id_ = stoi(tokens[1]);
    p_[0] = stod(tokens[2]);
    p_[1] = stod(tokens[3]);
    p_[2] = stod(tokens[4]);
    p_[3] = stod(tokens[5]);
    x_[0] = stod(tokens[6]);
    x_[1] = stod(tokens[7]);
    x_[2] = stod(tokens[8]);
    x_[3] = stod(tokens[9]);
    Nch_eta_ = stod(tokens[10]);

    wA_ = 1.;
    wLT_ = 1.;
    used_ = false;

}

double myParticle::pT() {
    return sqrt(p_[1] * p_[1] + p_[2] * p_[2]);
}

void myParticle::set_p(double E, double px, double py, double pz) {
    p_[0] = E;
    p_[1] = px;
    p_[2] = py;
    p_[3] = pz;
}

double myParticle::y() {
    return 0.5 * log((p_[0] + p_[3]) / (p_[0] - p_[3]));
}

void myParticle::set_x(double xt, double xx, double xy, double xz) {
    x_[0] = xt;
    x_[1] = xx;
    x_[2] = xy;
    x_[3] = xz;
}

myParticle myParticle::operator+(myParticle &proton) {
    myParticle Deuteron;
    double E = this->p_[0] + proton.p_[0];
    double px = this->p_[1] + proton.p_[1];
    double py = this->p_[2] + proton.p_[2];
    double pz = this->p_[3] + proton.p_[3];
    double wA = this->wA() * proton.wA();
    double wLT = this->wLT() * proton.wLT();

    Deuteron.set_p(E, px, py, pz);
    Deuteron.set_wA(wA);
    Deuteron.set_wLT(wLT);
    Deuteron.set_parent_y(this->y(), proton.y());
    Deuteron.set_Nch_eta(proton.Nch_eta_);

    return Deuteron;
}

double myParticle::wA() {
    return wA_;
}

void myParticle::set_wA(double w_init) {
    wA_ = w_init;
}

double myParticle::wLT() {
    return wLT_;
}

void myParticle::set_wLT(double w_init) {
    wLT_ = w_init;
}


void myParticle::init() {
    p_ = {0, 0, 0, 0};
    x_ = {0, 0, 0, 0};

    Nch_eta_ = 0;
    id_ = 0;
    wA_ = 1.;
    wLT_ = 1.;
    used_ = false;
    y_n_ = 0;
    y_p_ = 0;
}

bool myParticle::is_used() {
    return used_;
}

void myParticle::used() {
    used_ = true;
}

double myParticle::momentum_difference(myParticle &other) {
    double E = other.p_[0] - p_[0];
    double px = other.p_[1] - p_[1];
    double py = other.p_[2] - p_[2];
    double pz = other.p_[3] - p_[3];
    return -E * E + px * px + py * py + pz * pz;
}

double myParticle::spacetime_difference(myParticle &other) {
    double t = other.x_[0] - x_[0];
    double x = other.x_[1] - x_[1];
    double y = other.x_[2] - x_[2];
    double z = other.x_[3] - x_[3];
    return -t * t + x * x + y * y + z * z;
}

void myParticle::set_parent_y(double y_n, double y_p) {
    y_n_ = y_n;
    y_p_ = y_p;
}

double myParticle::y_p() {
    return y_p_;
}

double myParticle::y_n() {
    return y_n_;
}

void myParticle::set_Nch_eta(double Nch_eta) {
    Nch_eta_ = Nch_eta;
}

double myParticle::Nch_eta() {
    return Nch_eta_;
}

void myParticle::bstback(const myParticle& boost_particle) {
    if (abs(boost_particle.p_[0]) < 1e-20) {cout << "boost failed in p[0]< 1e-20" << endl; return;}
    double betaX = -boost_particle.p_[1] / boost_particle.p_[0];
    double betaY = -boost_particle.p_[2] / boost_particle.p_[0];
    double betaZ = -boost_particle.p_[3] / boost_particle.p_[0];
    double beta2 = betaX*betaX + betaY*betaY + betaZ*betaZ;
    if (beta2 >= 1.){ cout << "boost failed in beta2>=1." << endl; return;}

    //boost momentum components
    double gamma = 1. / sqrt(1. - beta2);
    double prod1 = betaX * p_[1] + betaY * p_[2] + betaZ * p_[3];
    double prod2 = gamma * (gamma * prod1 / (1. + gamma) + p_[0]);
    p_[1]          += prod2 * betaX;
    p_[2]          += prod2 * betaY;
    p_[3]          += prod2 * betaZ;
    p_[0]           = gamma * (p_[0] + prod1);

    //boost position components
    prod1 = betaX * x_[1] + betaY * x_[2] + betaZ * x_[3];
    prod2 = gamma * (gamma * prod1 / (1. + gamma) + x_[0]);
    x_[1]          += prod2 * betaX;
    x_[2]          += prod2 * betaY;
    x_[3]          += prod2 * betaZ;
    x_[0]           = gamma * (x_[0] + prod1);

}


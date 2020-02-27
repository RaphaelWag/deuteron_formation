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

void myParticle::bstback(const myParticle &boost_particle) {
    if (abs(boost_particle.p_[0]) < 1e-20) {
        cout << "boost failed in p[0]< 1e-20" << endl;
        return;
    }
    double betaX = -boost_particle.p_[1] / boost_particle.p_[0];
    double betaY = -boost_particle.p_[2] / boost_particle.p_[0];
    double betaZ = -boost_particle.p_[3] / boost_particle.p_[0];
    double beta2 = betaX * betaX + betaY * betaY + betaZ * betaZ;
    if (beta2 >= 1.) {
        cout << "boost failed in beta2>=1." << endl;
        return;
    }

    //boost momentum components
    double gamma = 1. / sqrt(1. - beta2);
    double prod1 = betaX * p_[1] + betaY * p_[2] + betaZ * p_[3];
    double prod2 = gamma * (gamma * prod1 / (1. + gamma) + p_[0]);
    p_[1] += prod2 * betaX;
    p_[2] += prod2 * betaY;
    p_[3] += prod2 * betaZ;
    p_[0] = gamma * (p_[0] + prod1);

    //boost position components
    prod1 = betaX * x_[1] + betaY * x_[2] + betaZ * x_[3];
    prod2 = gamma * (gamma * prod1 / (1. + gamma) + x_[0]);
    x_[1] += prod2 * betaX;
    x_[2] += prod2 * betaY;
    x_[3] += prod2 * betaZ;
    x_[0] = gamma * (x_[0] + prod1);

}

void myParticle::bst_E(double betaX, double betaY, double betaZ) {

    double beta2 = betaX * betaX + betaY * betaY + betaZ * betaZ;
    if (beta2 >= 1.) {
        cout << "boost failed in bst_E\n";
        return;
    }
    double gamma = 1. / sqrt(1. - beta2);
    double prod1 = betaX * p_[1] + betaY * p_[2] + betaZ * p_[3];
    //double prod2 = gamma * (gamma * prod1 / (1. + gamma) + p_[0]);

    p_[0] = gamma * (p_[0] + prod1);

    //p_[1] += prod2 * betaX;
    //p_[2] += prod2 * betaY;
    //p_[3] += prod2 * betaZ;
}


void myParticle::bst(double betaX, double betaY, double betaZ) {
    double beta2 = betaX * betaX + betaY * betaY + betaZ * betaZ;
    if (beta2 >= 1.) {
        cout << "boost failed in bst\n";
        return;
    }
    double gamma = 1. / sqrt(1. - beta2);
    double prod1 = betaX * p_[1] + betaY * p_[2] + betaZ * p_[3];
    double prod2 = gamma * (gamma * prod1 / (1. + gamma) + p_[0]);
    p_[0] = gamma * (p_[0] + prod1);
    p_[1] += prod2 * betaX;
    p_[2] += prod2 * betaY;
    p_[3] += prod2 * betaZ;
}

void myParticle::emit_single(double m, double r1, double r2) {

    double md2 = 3.51792391; //deuteron mass squared
    double m2 = m * m;
    double betaX, betaY, betaZ;

    //boost energy to COM frame

    betaX = -p_[1] / p_[0];
    betaY = -p_[2] / p_[0];
    betaZ = -p_[3] / p_[0];
    bst_E(betaX, betaY, betaZ);


    //in COM frame calculate new energy and momentum
    double E2 = p_[0] * p_[0];
    double p2 = -0.5 * (md2 + m2) - 0.5 * md2 * m2 / E2 + 0.25 * (E2 + (md2 * md2 + m2 * m2) / E2);

    //set momentum along x axis
    set_p(sqrt(md2 + p2), sqrt(p2), 0, 0);

    //rotate randomly
    double phi = 2. * M_PI * r1;
    double theta = acos(1. - 2. * r2);
    rotate(phi, theta);

    //boost to lab frame
    bst(-betaX, -betaY, -betaZ);

}

void myParticle::emit_double(double m1, double m2, double r1, double r2, double r3, double r4) {

    //1,2 = pi , 3 = d
    double md = 1.87561294257;
    double md2 = 3.51792391;

    double betaX = -p_[1] / p_[0];
    double betaY = -p_[2] / p_[0];
    double betaZ = -p_[3] / p_[0];
    bst_E(betaX, betaY, betaZ);

    double a = m1 + m2;
    double mpp2_min = a * a;
    a = p_[0] - md;
    double mpp2_max = a * a;

    double mpp2 = mpp2_min + r1 * (mpp2_max - mpp2_min);
    double mpp = sqrt(mpp2);

    double m12 = m1 * m1;
    double m22 = m2 * m2;
    double e = 0.5 * (mpp2 - m12 + m22) / mpp; //E_2
    double E = 0.5 * (p_[0] * p_[0] - mpp2 - md2) / mpp;// E_3
    a = e + E;
    double sqrt2 = sqrt(e * e - m22);
    double sqrt3 = sqrt(E * E - md2);
    double b = sqrt2 + sqrt3;
    double b1 = b * b;
    b = sqrt2 - sqrt3;
    double b2 = b * b;
    double mpd2 = a * a - b1 + r2 * (b1 - b2);
    double s = p_[0] * p_[0];
    b = p_[0] * p_[0] + md2 - mpd2;

    double pd2 = 0.25 / s * b * b - md2;
    set_p(sqrt(md2 + pd2), sqrt(pd2), 0, 0);

    double phi = 2. * M_PI * r3;
    double theta = acos(1. - 2. * r4);
    rotate(phi, theta);

    //boost to lab frame
    bst(-betaX, -betaY, -betaZ);
}

void myParticle::rotate(double phi, double theta) {
    double cthe = cos(theta);
    double sthe = sin(theta);
    double cphi = cos(phi);
    double sphi = sin(phi);
    double tmpx = cthe * cphi * p_[1] - sphi * p_[2] + sthe * cphi * p_[3];
    double tmpy = cthe * sphi * p_[1] + cphi * p_[2] + sthe * sphi * p_[3];
    double tmpz = -sthe * p_[1] + cthe * p_[3];
    p_[1] = tmpx;
    p_[2] = tmpy;
    p_[3] = tmpz;
}

bool myParticle::operator==(myParticle &other) {
    if(this->id_!=other.id_){ return false;}
    if(this->p_!=other.p_){ return false;}
    if(this->x_!=other.x_){ return false;}
    if(this->wA_!=other.wA_){ return false;}
    if(this->wLT_!=other.wLT_){ return false;}
    if(this->Nch_eta_!=other.Nch_eta_){ return false;}
    if(this->used_!=other.used_){ return false;}
    if(this->y_n_!=other.y_n_){ return false;}
    return this->y_p_ == other.y_p_;
}

void myParticle::move(double t) {
    double betaX, betaY, betaZ;
    betaX = p_[1]/p_[0];
    betaY = p_[2]/p_[0];
    betaZ = p_[3]/p_[0];

    x_[0] += t;
    x_[1] += t*betaX;
    x_[2] += t*betaY;
    x_[3] += t*betaZ;
}






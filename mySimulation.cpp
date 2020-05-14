//
// Created by rwagn on 07.03.2019.
//

#include "mySimulation.h"
#include <istream>
#include <sstream>

using namespace std;

//seed rng
random_device rd;
mt19937 gen(rd());
uniform_real_distribution<> random01(0, 1);

mySimulation::mySimulation() {
    N_events_ = 0;
    //Particles = new vector<myParticle> *[N_events_];
    cutoff_momentum_2_ = 0.195 * 0.195;
    cutoff_position_2_ = pow(10, -20);


}

void mySimulation::load_txt(const string &file, int N_Events, bool reduced_data) {
    cout << file << endl;

    auto myParticle_constructor_ = new string[11];
    string line;
    string token;

    //destruct old arrays
    Event.clear();
    //construct new particle array
    N_events_ = N_Events;


    //read file and construct particles
    ifstream myfile(file);
    bool firstline = true;
    int event_ = 0;
    event temp_event;

    if (!myfile.is_open()) {
        cout << "Error in: mySimulation -> load_txt" << "\n";
        cout << "Unable to open file" << "\n";
    } else {
        while (getline(myfile, line)) {

            if (firstline and reduced_data) {
                firstline = false;
            } else {
                int i = 0;//loop variable
                stringstream ss(line);
                while (getline(ss, token, ' ')) {
                    myParticle_constructor_[i] = token;
                    i++;
                }
                int new_event = stoi(myParticle_constructor_[0]);
                if (new_event == event_) {
                    event_ = new_event;
                    //create new particle
                    myParticle new_particle(myParticle_constructor_);

                    //add particle to Particles array
                    if (abs(new_particle.id()) == p_id_) {
                        temp_event.protons.push_back(new_particle);
                    }
                    if (abs(new_particle.id()) == n_id_) {
                        temp_event.neutrons.push_back(new_particle);
                    }
                }
                if (new_event != event_) {
                    event_ = new_event;
                    Event.push_back(temp_event);
                    temp_event.protons.clear();
                    temp_event.neutrons.clear();
                    myParticle new_particle(myParticle_constructor_);
                    //add particle to Particles array
                    if (abs(new_particle.id()) == p_id_) {
                        temp_event.protons.push_back(new_particle);
                    }
                    if (abs(new_particle.id()) == n_id_) {
                        temp_event.neutrons.push_back(new_particle);
                    }
                }
            }

        }
        myfile.close();
        Event.push_back(temp_event);
    }
    delete[]myParticle_constructor_;
}

void mySimulation::coalescence() {
    //delete old particles
    deuteron.clear();
    //form new ones
    for (auto &element:Event) {
        form_deuterons(element.protons, element.neutrons);
    }
}

void mySimulation::coalescence_px() {
    //delete old particles
    deuteron.clear();
    //form new ones
    for (auto &element:Event) {
        form_deuterons_px(element.protons, element.neutrons);
    }
}


void mySimulation::cs_model_formation() {
    //delete old particles
    deuteron.clear();
    //form new ones
    for (auto &element:Event) {
        form_deuterons_cs(element.protons, element.neutrons);
    }
}

bool mySimulation::check_momentum(myParticle &neutron, myParticle &proton) {
    double delta = proton.momentum_difference(neutron);
    return delta <= cutoff_momentum_2_;
}

bool mySimulation::check_position(myParticle neutron, myParticle proton) {

    //boost to CoM frame
    myParticle boostparticle = neutron + proton;
    neutron.bstback(boostparticle);
    proton.bstback(boostparticle);

    double d2;
    double t = neutron.x(0) - proton.x(0);
    if (t > 0) {
        proton.move(t);
    } else {
        neutron.move(-t);
    }

    // in CoM frame get closest distance
    //check if they are coming closer or not
    if (coming_closer(proton, neutron)) {
        //get closest distance
        // cross product of direction vectors (= momentum vectors)
        double nx = proton.p(2) * neutron.p(3) - proton.p(3) * neutron.p(2);
        double ny = proton.p(3) * neutron.p(1) - proton.p(1) * neutron.p(3);
        double nz = proton.p(1) * neutron.p(2) - proton.p(2) * neutron.p(1);

        // r1 - r2 = r
        double rx = proton.x(1) - neutron.x(1);
        double ry = proton.x(2) - neutron.x(2);
        double rz = proton.x(3) - neutron.x(3);

        // d =n*r/n
        double a = abs(nx * rx + ny * ry + nz * rz);
        d2 = a * a / (nx * nx + ny * ny + nz * nz);

    } else {
        // r1 - r2 = r
        double rx = proton.x(1) - neutron.x(1);
        double ry = proton.x(2) - neutron.x(2);
        double rz = proton.x(3) - neutron.x(3);

        d2 = rx * rx + ry * ry + rz * rz;
    }

    return d2 <= cutoff_position_2_;

}

void
mySimulation::form_deuterons(vector<myParticle> Protons, vector<myParticle> Neutrons) {
    myParticle d;
    for (auto &proton:Protons) {
        for (auto &neutron:Neutrons) {
            if (proton.is_used()) { continue; }
            if (neutron.is_used()) { continue; }
            if (!check_momentum(neutron, proton)) { continue; }

            d = neutron + proton;
            d.emit_single(0, random01(gen), random01(gen));
            deuteron.push_back(d);
            proton.used();
            neutron.used();
        }
    }
}

void mySimulation::form_deuterons_px(vector<myParticle> Protons, vector<myParticle> Neutrons) {
    myParticle d;
    for (auto &proton:Protons) {
        for (auto &neutron:Neutrons) {
            if (proton.is_used()) { continue; }
            if (neutron.is_used()) { continue; }
            if (!check_momentum(neutron, proton)) { continue; }
            if (!check_position(neutron, proton)) { continue; }

            d = neutron + proton;
            d.emit_single(0, random01(gen), random01(gen));
            deuteron.push_back(d);
            proton.used();
            neutron.used();
        }
    }
}

void mySimulation::set_ALICE_weights_lt(const string &cms, bool particle_type) {
    string file;
    if (particle_type) {
        file = "/mnt/d/Uni/Lectures/thesis/ALICE_weighting/levy_tsallis_weights/weights_p_" + cms + ".txt";
    } else {
        file = "/mnt/d/Uni/Lectures/thesis/ALICE_weighting/levy_tsallis_weights/weights_pbar_" + cms + ".txt";
    }

    load_weights(file, weights_LT);

}

void mySimulation::load_weights(const string &file, weights &weight_struct) {
    string line;
    string token;

    ifstream myfile(file);

    bool initialized = false;
    int k = 0;
    int N_weights = 1;

    //read in weights
    if (!myfile.is_open()) {
        cout << "Error in: mySimulation -> set_ALICE_weights_discrete" << "\n";
        cout << "Unable to open file" << "\n";
    } else {
        while ((getline(myfile, line)) and (k < N_weights)) {
            if (!initialized) {

                N_weights = stoi(line);

                weight_struct.weight_boundaries_.clear();
                weight_struct.weight_centre_.clear();
                weight_struct.weights_.clear();
                initialized = true;
            } else {
                stringstream ss(line);
                int j = 0;
                while (getline(ss, token, ' ')) {
                    if (!j) {
                        weight_struct.weight_centre_.push_back(stod(token));
                        j++;
                    } else {
                        weight_struct.weights_.push_back(stod(token));
                    }
                }
            }
        }
        myfile.close();
    }

    //construct weight boundaries
    weight_struct.weight_boundaries_.push_back(0);
    weight_struct.weight_boundaries_.push_back(
            0.5 * (weight_struct.weight_centre_[0] + weight_struct.weight_centre_[1]));
    weight_struct.weight_boundaries_[0] =
            weight_struct.weight_centre_[0] - (weight_struct.weight_boundaries_[1] - weight_struct.weight_centre_[0]);

    for (int i = 2; i < N_weights + 1; ++i) {
        weight_struct.weight_boundaries_.push_back(weight_struct.weight_centre_[i - 1] +
                                                   (weight_struct.weight_centre_[i - 1] -
                                                    weight_struct.weight_boundaries_[i - 1]));
    }
    weight_struct.rescaled = true;
}

void mySimulation::rescale_spectrum() {

    double pT = 0;

    if (!(weights_LT.rescaled)) {
        cout << "Rescaling weights not set" << "\n";
    } else {
        for (auto &element:Event) { //loop for all events
            for (auto &pbar: element.protons) { //loop for all pbar in event
                pT = pbar.pT();
                //then rescale with LT weights
                for (int j = 0; j < weights_LT.weights_.size(); ++j) { //loop for get the correct weight
                    if (pT < weights_LT.weight_boundaries_[j]) { continue; }
                    if (pT > weights_LT.weight_boundaries_[j + 1]) { continue; }
                    pbar.set_wLT(weights_LT.weights_[j]);
                }
            }
            for (auto &nbar: element.neutrons) {//loop for all nbar in event
                pT = nbar.pT();
                //then rescale with LT weights
                for (int j = 0; j < weights_LT.weights_.size(); ++j) { //loop for get the correct weight
                    if (pT < weights_LT.weight_boundaries_[j]) { continue; }
                    if (pT > weights_LT.weight_boundaries_[j + 1]) { continue; }
                    nbar.set_wLT(weights_LT.weights_[j]);
                }
            }
        }
    }
}

void mySimulation::set_cutoff_momentum(double cutoff) {
    cutoff_momentum_2_ = cutoff * cutoff;
}

void mySimulation::set_cutoff_position(double cutoff) {
    cutoff_position_2_ = cutoff * cutoff;
}

void mySimulation::load_txt(const string &file, int N_Events, int particle_id) {
    auto myParticle_constructor_ = new string[11];
    string line;
    string token;

    //destruct old arrays
    Event.clear();
    //construct new particle array
    N_events_ = N_Events;

    int event_ = -1;
    event temp_event;
    //read file and construct particles
    ifstream myfile(file);
    cout << file << "\n";
    if (!myfile.is_open()) {
        cout << "Error in: mySimulation -> load_txt" << "\n";
        cout << "Unable to open file" << "\n";
    } else {
        while (getline(myfile, line)) {

            int i = 0;//loop variable
            stringstream ss(line);
            while (getline(ss, token, ' ')) {
                myParticle_constructor_[i] = token;
                i++;
            }
            int new_event = stoi(myParticle_constructor_[0]);
            if (new_event == event_) {
                event_ = new_event;
                //create new particle
                myParticle new_particle(myParticle_constructor_);

                //add particle to Particles array
                if (new_particle.id() == particle_id) {
                    temp_event.protons.push_back(new_particle);
                }
            } else {
                event_ = new_event;
                Event.push_back(temp_event);
                temp_event.protons.clear();
                temp_event.neutrons.clear();
                myParticle new_particle(myParticle_constructor_);
                if (new_particle.id() == particle_id) {
                    temp_event.protons.push_back(new_particle);
                }
            }

        }
        myfile.close();
        Event.push_back(temp_event);
    }

    delete[]myParticle_constructor_;
}


void
mySimulation::form_deuterons_cs(vector<myParticle> Protons, vector<myParticle> Neutrons) {

    double m1 = 0.1349770; //neutral pion mass
    double m2 = 0.13957061; // charged pion mass
    double m2_inv = 7.1648322; // inverse charged pion mass
    double md2 = 3.51792391; //deuteron mass squared
    double md = 1.87561294257; //deuteron mass
    double betaX, betaY, betaZ;
    double E, k, dx, dy, dz;
    myParticle N1, N2; //temporary particles
    temp_d process{};
    vector<temp_d> still_possible_deuterons;
    temporary_deuteron.clear();
    bool p = false;
    bool n = false;
    if (!Protons.empty()) { p = true; }
    if (!Neutrons.empty()) { n = true; }

    //check if energy is sufficient
    //loop over pn -> dX
    if (p and n) {
        for (int i = 0; i < Protons.size(); ++i) {
            for (int j = 0; j < Neutrons.size(); ++j) {
                //general computations for all processes
                N1 = Protons[i];
                N2 = Neutrons[j];
                E = N1.p(0) + N2.p(0);
                betaX = -(N1.p(1) + N2.p(1)) / E;
                betaY = -(N1.p(2) + N2.p(2)) / E;
                betaZ = -(N1.p(3) + N2.p(3)) / E;
                N1.bst(betaX, betaY, betaZ);
                N2.bst(betaX, betaY, betaZ);
                dx = N1.p(1) - N2.p(1);
                dy = N1.p(2) - N2.p(2);
                dz = N1.p(3) - N2.p(3);
                k = sqrt(dx * dx + dy * dy + dz * dz);
                process.n1 = i;
                process.n2 = j;

                E = N1.p(0) + N2.p(0); // update energy in com frame
                double E2 = E * E;
                //calculate pion momentum

                double m22 = m2 * m2;
                double q = sqrt(-0.5 * (md2 + m22) - 0.5 * md2 * m22 / E2 + 0.25 * (E2 + (md2 * md2 + m22 * m22) / E2));


                //compute individual probabilities and save informations
                process.type = 1;
                process.probability = prob_1(k);
                temporary_deuteron.push_back(process);

                if (E > md + m1) {
                    process.type = 2;
                    process.probability = prob_2(q, m2_inv);
                    temporary_deuteron.push_back(process);

                    if (E > md + m1 + m1) {
                        process.type = 3;
                        process.probability = prob_3(k);
                        temporary_deuteron.push_back(process);

                        if (E > md + m2 + m2) {
                            process.type = 4;
                            process.probability = prob_4(k);
                            temporary_deuteron.push_back(process);
                        }
                    }
                }
            }
        }
    }
    if (p) {
        //loop over pp -> dX
        for (int ii = 0; ii < Protons.size() - 1; ++ii) {
            for (int jj = ii + 1; jj < Protons.size(); ++jj) {
                //general computations for all processes
                N1 = Protons[ii];
                N2 = Protons[jj];
                E = N1.p(0) + N2.p(0);
                betaX = -(N1.p(1) + N2.p(1)) / E;
                betaY = -(N1.p(2) + N2.p(2)) / E;
                betaZ = -(N1.p(3) + N2.p(3)) / E;

                N1.bst(betaX, betaY, betaZ);
                N2.bst(betaX, betaY, betaZ);
                dx = N1.p(1) - N2.p(1);
                dy = N1.p(2) - N2.p(2);
                dz = N1.p(3) - N2.p(3);
                k = sqrt(dx * dx + dy * dy + dz * dz);
                process.n1 = ii;
                process.n2 = jj;

                E = N1.p(0) + N2.p(0); // update energy in com frame
                double E2 = E * E;

                //calculate pion momentum
                double m22 = m2 * m2;
                double q = sqrt(-0.5 * (md2 + m22) - 0.5 * md2 * m22 / E2 + 0.25 * (E2 + (md2 * md2 + m22 * m22) / E2));

                //compute individual probabilities and save informations
                if (E > md + m2) {
                    process.type = 5;
                    process.probability = prob_5(q, m2_inv);
                    temporary_deuteron.push_back(process);

                    if (E > md + m1 + m2) {
                        process.type = 6;
                        process.probability = prob_6(k);
                        temporary_deuteron.push_back(process);
                    }
                }
            }
        }
    }
    if (n) {
        //loop over nn -> dX
        for (int I = 0; I < Neutrons.size() - 1; ++I) {
            for (int J = I + 1; J < Neutrons.size(); ++J) {
                //general computations for all processes
                N1 = Neutrons[I];
                N2 = Neutrons[J];
                E = N1.p(0) + N2.p(0);
                betaX = -(N1.p(1) + N2.p(1)) / E;
                betaY = -(N1.p(2) + N2.p(2)) / E;
                betaZ = -(N1.p(3) + N2.p(3)) / E;
                N1.bst(betaX, betaY, betaZ);
                N2.bst(betaX, betaY, betaZ);
                dx = N1.p(1) - N2.p(1);
                dy = N1.p(2) - N2.p(2);
                dz = N1.p(3) - N2.p(3);
                k = sqrt(dx * dx + dy * dy + dz * dz);
                process.n1 = I;
                process.n2 = J;

                E = N1.p(0) + N2.p(0); // update energy in com frame
                double E2 = E * E;
                //calculate pion momentum
                double m22 = m2 * m2;
                double q = sqrt(-0.5 * (md2 + m22) - 0.5 * md2 * m22 / E2 + 0.25 * (E2 + (md2 * md2 + m22 * m22) / E2));

                //compute individual probabilities and save informations
                if (E > md + m2) {
                    process.type = 7;
                    process.probability = prob_7(q, m2_inv);
                    temporary_deuteron.push_back(process);
                    if (E > md + m1 + m2) {
                        process.type = 8;
                        process.probability = prob_8(k);
                        temporary_deuteron.push_back(process);
                    }
                }
            }
        }
    }

    //draw all r_i for processes and check if r_i < p_i
    vector<temp_d> possible_deuterons;
    vector<double> probabilities;

    for (auto &element:temporary_deuteron) {
        if (element.probability * sigma_0_inv_ > random01(gen)) { possible_deuterons.push_back(element); }
    }

    while (!possible_deuterons.empty()) {
        double norm = 0;
        double prob;
        probabilities = {0};
        for (auto &element:possible_deuterons) {
            prob = element.probability * sigma_0_inv_;
            if (prob > 1.) { prob = 1.; }
            norm += prob;
            probabilities.push_back(norm);
        }

        double r = random01(gen) * norm;
        // select one process from all possible processes
        for (int m = 0; m < probabilities.size() - 1; ++m) {
            if (probabilities[m + 1] <= r) { continue; }
            if (probabilities[m] >= r) { continue; }
            myParticle d;
            int i = possible_deuterons[m].n1;
            int j = possible_deuterons[m].n2;

            switch (possible_deuterons[m].type) {
                case 1: // pn -> d gamma
                    d = Protons[i] + Neutrons[j];
                    Protons[i].used();
                    Neutrons[j].used();
                    d.emit_single(0., random01(gen), random01(gen));
                    goto done;
                case 2: // pn -> d pi0
                    d = Protons[i] + Neutrons[j];
                    Protons[i].used();
                    Neutrons[j].used();
                    d.emit_single(m1, random01(gen), random01(gen));
                    goto done;
                case 3: // pn -> d pi0 pi0
                    d = Protons[i] + Neutrons[j];
                    Protons[i].used();
                    Neutrons[j].used();
                    d.emit_double(m1, m1, random01(gen), random01(gen), random01(gen), random01(gen));
                    goto done;
                case 4: // pn -> d pi+ pi-
                    d = Protons[i] + Neutrons[j];
                    Protons[i].used();
                    Neutrons[j].used();
                    d.emit_double(m2, m2, random01(gen), random01(gen), random01(gen), random01(gen));
                    goto done;
                case 5: // pp -> d pi-
                    d = Protons[i] + Protons[j];
                    Protons[i].used();
                    Protons[j].used();
                    d.emit_single(m2, random01(gen), random01(gen));
                    goto done;
                case 6: // pp -> d pi- pi0
                    d = Protons[i] + Protons[j];
                    Protons[i].used();
                    Protons[j].used();
                    d.emit_double(m1, m2, random01(gen), random01(gen), random01(gen), random01(gen));
                    goto done;
                case 7: // nn -> d pi+
                    d = Neutrons[i] + Neutrons[j];
                    Neutrons[i].used();
                    Neutrons[j].used();
                    d.emit_single(m2, random01(gen), random01(gen));
                    goto done;
                case 8: // nn -> d pi+ pi0
                    d = Neutrons[i] + Neutrons[j];
                    Neutrons[i].used();
                    Neutrons[j].used();
                    d.emit_double(m1, m2, random01(gen), random01(gen), random01(gen), random01(gen));
            }
            done:
            deuteron.push_back(d);
            // loop over possible deuterons again
            // delete all possibilities which contain alrady used partcles
            // repeat until possible deuteron container is empty

            for (auto &element:possible_deuterons) {
                switch (element.type) {
                    case 1:
                    case 2:
                    case 3:
                    case 4: // all pn processes
                        if (Protons[element.n1].is_used()) { continue; }
                        if (Neutrons[element.n2].is_used()) { continue; }
                        still_possible_deuterons.push_back(element);
                        continue;
                    case 5:
                    case 6:
                        if (Protons[element.n1].is_used()) { continue; }
                        if (Protons[element.n2].is_used()) { continue; }
                        still_possible_deuterons.push_back(element);
                        continue;
                    case 7:
                    case 8:
                        if (Neutrons[element.n1].is_used()) { continue; }
                        if (Neutrons[element.n2].is_used()) { continue; }
                        still_possible_deuterons.push_back(element);
                        continue;
                }
            }
        }
        possible_deuterons = still_possible_deuterons;
        still_possible_deuterons.clear();
    }
}

double mySimulation::prob_1(double k) { //pn -> d gamma
    double sigma = 0;

    double b1 = -5.1885;
    double b2 = 2.9196;

    vector<double> A = {2.30346, -9.366346e1, 2.565390e3,
                        -2.5594101e4, 1.43513109e5, -5.03572889e5,
                        1.14924802e6, -1.72368391e6, 1.6793476e6,
                        -1.01988855e6, 3.4984035e5, -5.1662760e4};

    if (k < 1.28) {
        double value = 1. / k;
        for (auto &a:A) {
            sigma += a * value;
            value = value * k;
        }
    } else { sigma = exp(-(b1 * k + b2 * k * k)); }

    return sigma;
}

double mySimulation::prob_2(double q, double inverse_pi_mass) {// pn -> d pi0
    return 0.5 * prob_5(q, inverse_pi_mass);
}

double mySimulation::prob_3(double k) { // pn -> d pi0 pi0
    return prob(k, 2.855e6, 1.311e1, 2.961e3, 5.572, 1.462e6);
}

double mySimulation::prob_4(double k) { // pn -> d pi+ pi-
    return prob(k, 6.465e6, 1.051e1, 1.979e3, 5.363, 6.045e5) + prob(k, 2.549e15, 1.657e1, 2.330e7, 1.119e1, 2.868e16);
}

double mySimulation::prob_5(double q, double inverse_pi_mass) { // pp -> d pi+
    double eta = q * inverse_pi_mass;
    return prob(eta, 170, 1.34, 1.77, 0.38, 0.096);
}

double mySimulation::prob_6(double k) { // pp -> d pi+ pi0
    return prob(k, 5.099e15, 1.656e1, 2.333e7, 1.133e1, 2.868e16);
}

double mySimulation::prob_7(double q, double inverse_pi_mass) { // nn -> d pi-
    return prob_5(q, inverse_pi_mass);
}

double mySimulation::prob_8(double k) { // nn -> d pi- pi0
    return prob_6(k);
}

double mySimulation::prob(double k, double a, double b, double c, double d, double e) {
    double temp = c - exp(d * k);
    return a * pow(k, b) / (temp * temp + e);
}

void mySimulation::set_sigma0(double sigma0) {
    sigma_0_inv_ = sigma0;
}

int mySimulation::N_events() {
    return N_events_;
}

bool mySimulation::coming_closer(myParticle &proton, myParticle &neutron) {
    double x, y, z;
    x = neutron.x(1) - proton.x(1);
    y = neutron.x(2) - proton.x(2);
    z = neutron.x(3) - proton.x(3);

    double sp = x * proton.p(1) + y * proton.p(2) + z * proton.p(3);

    return sp > 0;
}

void mySimulation::set_ALICE_weights_multi(const string &cms, bool particle_type, string *class_name, int N_classes) {

    multi_weights.clear();

    string file;
    for (int i = 0; i < N_classes; ++i) {
        weights new_weight;
        if (particle_type) {
            file = "/mnt/d/Uni/Lectures/thesis/ALICE_weighting/levy_tsallis_weights/weights_p_" + cms + "_" +
                   class_name[i] + ".txt";
        } else {
            file = "/mnt/d/Uni/Lectures/thesis/ALICE_weighting/levy_tsallis_weights/weights_pbar_" + cms + "_" +
                   class_name[i] + ".txt";
        }

        load_weights(file, new_weight);
        multi_weights.push_back(new_weight);
    }

}

void mySimulation::rescale_spectrum_multi(vector<double> class_limits) {
    double pT = 0;

    for (auto &element:Event) { //loop for all events
        for (auto &pbar: element.protons) { //loop for all pbar in event
            pT = pbar.pT();

            //then rescale with weights
            for (int i = 0; i < class_limits.size() - 1; ++i) {
                //check for correct dNch_deta
                if (pbar.Nch_eta() > class_limits[i]) { continue; }
                if (pbar.Nch_eta() < class_limits[i + 1]) { continue; }

                weights weight = multi_weights[i];
                if (!weight.rescaled) { cout << "weights not set" << endl; }
                for (int j = 0; j < weight.weights_.size(); ++j) { //loop for get the correct weight
                    //check for correct pT
                    if (pT < weight.weight_boundaries_[j]) { continue; }
                    if (pT > weight.weight_boundaries_[j + 1]) { continue; }
                    pbar.set_wLT(weight.weights_[j]);
                }
            }
        }
        for (auto &nbar: element.neutrons) {//loop for all nbar in event
            pT = nbar.pT();
            //then rescale with weights
            for (int i = 0; i < class_limits.size() - 1; ++i) {
                //check for correct dNch_deta
                if (nbar.Nch_eta() > class_limits[i]) { continue; }
                if (nbar.Nch_eta() < class_limits[i + 1]) { continue; }

                weights weight = multi_weights[i];
                for (int j = 0; j < weight.weights_.size(); ++j) { //loop for get the correct weight
                    //check for correct pT
                    if (pT < weight.weight_boundaries_[j]) { continue; }
                    if (pT > weight.weight_boundaries_[j + 1]) { continue; }
                    nbar.set_wLT(weight.weights_[j]);
                }
            }
        }
    }
}





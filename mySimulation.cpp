//
// Created by rwagn on 07.03.2019.
//

#include "mySimulation.h"
#include <istream>
#include <sstream>

using namespace std;

mySimulation::mySimulation() {
    N_events_ = 0;
    Particles = new vector<myParticle> *[N_events_];

    cutoff_momentum_2_ = 0.195 * 0.195;
    cutoff_position_2_ = pow(10, -20);

}

void mySimulation::load_txt(const string &file, int N_Events) {

    auto myParticle_constructor_ = new string[11];
    string line;
    string token;

    //destruct old arrays

    for (int l = 0; l < N_events_; ++l) {
        for (int i = 0; i < id_array_.size(); ++i) {
            Particles[l][i].clear();
        }
    }

    for (int k = 0; k < N_events_; ++k) {
        delete[]Particles[k];
    }
    delete[]Particles;

    //construct new particle array
    N_events_ = N_Events;


    Particles = new vector<myParticle> *[N_events_];
    for (int i = 0; i < N_events_; ++i) {
        Particles[i] = new vector<myParticle>[id_array_.size()];
    }

    //read file and construct particles
    ifstream myfile(file);

    if (!myfile.is_open()) {
        cout << "Error in: mySimulation -> load_txt" << "\n";
        cout << "Unable to open file" << "\n";
    } else {
        while (getline(myfile, line)) {

            int i = 0;//loop variable
            int event;
            stringstream ss(line);
            while (getline(ss, token, ' ')) {
                myParticle_constructor_[i] = token;
                i++;
            }
            event = stoi(myParticle_constructor_[0]);

            //create new particle
            myParticle new_particle(myParticle_constructor_);

            //add particle to Particles array
            for (unsigned j = 0; j < id_array_.size(); ++j) {
                if (abs(new_particle.id()) == id_array_[j]) {
                    Particles[event][j].push_back(new_particle);
                }
            }
        }
        myfile.close();
    }

    delete[]myParticle_constructor_;
}

void mySimulation::coalescence() {

    //delete old particles
    deuteron.clear();
    //form new ones
    for (int k = 0; k < N_events_; ++k) {
        form_deuterons(Particles[k][0], Particles[k][1], deuteron);
    }
}

bool mySimulation::cascade_momentum(myParticle &neutron, myParticle &proton) {
    double delta = proton.momentum_difference(neutron);
    return delta <= cutoff_momentum_2_;
}

bool mySimulation::cascade_position(myParticle &neutron, myParticle &proton) {

    double delta = proton.spacetime_difference(neutron);

    //return delta <= cutoff_position_2_;

    return delta < 0;

    //return true;
}

void
mySimulation::form_deuterons(vector<myParticle> Protons, vector<myParticle> Neutrons, vector<myParticle> &Deuterons) {

    for (auto &proton:Protons) {
        for (auto &neutron:Neutrons) {
            if (proton.is_used()) { continue; }
            if (neutron.is_used()) { continue; }
            if (!cascade_momentum(neutron, proton)) { continue; }
            // if (!cascade_position(neutron, proton)) { continue; }
            Deuterons.push_back(neutron + proton);
            proton.used();
            neutron.used();
        }
    }
}

void mySimulation::set_ALICE_weights_discrete(const string &cms, bool particle_type) {

    string file;
    if (particle_type) {
        file = "/mnt/d/Uni/Lectures/thesis/ALICE_weighting/discrete_weights/weights_p_" + cms + ".txt";
    } else {
        file = "/mnt/d/Uni/Lectures/thesis/ALICE_weighting/discrete_weights/weights_pbar_" + cms + ".txt";
    }

    load_weights(file, weights_A);

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
    weight_struct.weight_boundaries_.push_back(0.5 * (weight_struct.weight_centre_[0] + weight_struct.weight_centre_[1]));
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

    if (!(weights_LT.rescaled and weights_A.rescaled)) {
        cout << "Rescaling weights not set" << "\n";
    } else {
        for (int i = 0; i < N_events_; ++i) { //loop for all events
            for (auto &pbar: Particles[i][0]) { //loop for all pbar in event
                pT = pbar.pT();
                //first recale with ALICE weights
                for (int j = 0; j < weights_A.weights_.size(); ++j) { //loop for get the correct weight
                    if (pT < weights_A.weight_boundaries_[j]) { continue; }
                    if (pT > weights_A.weight_boundaries_[j + 1]) { continue; }
                    pbar.set_wA(weights_A.weights_[j]);
                }
                //then rescale with LT weights
                for (int j = 0; j < weights_LT.weights_.size(); ++j) { //loop for get the correct weight
                    if (pT < weights_LT.weight_boundaries_[j]) { continue; }
                    if (pT > weights_LT.weight_boundaries_[j + 1]) { continue; }
                    pbar.set_wLT(weights_LT.weights_[j]);
                }
            }
            for (auto &nbar: Particles[i][1]) {//loop for all nbar in event
                pT = nbar.pT();
                //first recale with ALICE weights
                for (int j = 0; j < weights_A.weights_.size(); ++j) { //loop for get the correct weight
                    if (pT < weights_A.weight_boundaries_[j]) { continue; }
                    if (pT > weights_A.weight_boundaries_[j + 1]) { continue; }
                    nbar.set_wA(weights_A.weights_[j]);
                }
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

    for (int l = 0; l < N_events_; ++l) {
        for (int i = 0; i < id_array_.size(); ++i) {
            Particles[l][i].clear();
        }
    }

    for (int k = 0; k < N_events_; ++k) {
        delete[]Particles[k];
    }
    delete[]Particles;

    //construct new particle array
    N_events_ = N_Events;


    Particles = new vector<myParticle> *[N_events_];
    for (int i = 0; i < N_events_; ++i) {
        Particles[i] = new vector<myParticle>[1];
    }

    //read file and construct particles
    ifstream myfile(file);

    if (!myfile.is_open()) {
        cout << "Error in: mySimulation -> load_txt" << "\n";
        cout << "Unable to open file" << "\n";
    } else {
        while (getline(myfile, line)) {

            int i = 0;//loop variable
            int event;
            stringstream ss(line);
            while (getline(ss, token, ' ')) {
                myParticle_constructor_[i] = token;
                i++;
            }
            event = stoi(myParticle_constructor_[0]);

            //create new particle
            myParticle new_particle(myParticle_constructor_);

            //add particle to Particles array

            if (new_particle.id() == particle_id) {
                Particles[event][0].push_back(new_particle);
            }
        }
        myfile.close();
    }

    delete[]myParticle_constructor_;
}

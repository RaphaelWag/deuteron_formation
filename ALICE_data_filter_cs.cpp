//
// Created by rwagn on 30.01.2020.
//

#include "ALICE_input.h"
#include <iostream>
#include <fstream>
#include "myParticle.h"
#include <vector>

using namespace std;

void load_txt(const string &file, vector <myParticle> **Particles);

void write_txt(const string &file, vector<int> &selected_events, vector <myParticle> **Particles, bool p_type,
               int events_total);

int event_selector(int N_Events, vector <myParticle> **Particles, vector<int> &selected_events, bool p_type);


bool cascade_momentum(myParticle &neutron, myParticle &proton);


int main() {

    Input input;
    input.full_data();
    double cms = 0.9;
    input.set_cms(cms);
    int max_events;

    string data_dir = "/mnt/d/Uni/Lectures/thesis/ALICE_datasets_pythia/";
    string target_np_dir = data_dir + "np_pairs_cs/";
    string target_npbar_dir = data_dir + "npbar_pairs_cs/";

    vector <myParticle> **Particles = nullptr; //Particle[Event][Type][i-th particle of this event and this typ]
    //Type: 0 proton, 1 neutron, 2 anti proton, 3 anti neutron

    vector<int> selected_events;

    for (int i = 0; i < input.N_simulations; ++i) {

        //construct new particles
        Particles = new vector <myParticle> *[input.N_events[i]];
        for (int x = 0; x < input.N_events[i]; ++x) {
            Particles[x] = new vector<myParticle>[4];
        }

        load_txt(input.dataset_folder + input.files[i], Particles);
        max_events = event_selector(input.N_events[i], Particles, selected_events, true);
        write_txt(target_np_dir + input.files[i], selected_events, Particles, true, max_events);
        max_events = event_selector(input.N_events[i], Particles, selected_events, false);
        write_txt(target_npbar_dir + input.files[i], selected_events, Particles, false, max_events);

        //delete old Particles
        for (int l = 0; l < input.N_events[i]; ++l) {
            for (int k = 0; k < 4; ++k) {
                Particles[l][k].clear();
            }
        }

        for (int m = 0; m < input.N_events[i]; ++m) {
            delete[]Particles[m];
        }
        delete[]Particles;
    }

    return 0;
}

void load_txt(const string &file, vector <myParticle> **Particles) {

    auto myParticle_constructor_ = new string[11];
    string line;
    string token;
    int n_id_ = 2112; // neutron
    int N_id_ = -n_id_; // anti neutron
    int p_id_ = 2212; // proton
    int P_id_ = -p_id_; // anti proton
    int id_array_[4] = {p_id_, n_id_, P_id_, N_id_};

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
            for (int j = 0; j < 4; ++j) {
                if (new_particle.id() == id_array_[j]) {
                    Particles[event][j].push_back(new_particle);
                }
            }
        }
        myfile.close();
    }

    delete[]myParticle_constructor_;
}

void write_txt(const string &file, vector<int> &selected_events, vector <myParticle> **Particles, bool p_type,
               int events_total) {

    vector<int> p_types = {}; // if true use particles if false use anti particles

    if (p_type) {
        p_types = {0, 1};
    } else {
        p_types = {2, 3};
    }

    ofstream myfile(file);
    if (!myfile.is_open()) {
        cout << "unable to open file in write_txt" << "\n";
    } else {
        myfile << events_total << endl;
        int j = 0;
        for (auto &i:selected_events) {
            for (auto &type:p_types)
                for (auto &Particle:Particles[i][type]) {
                    myfile << j << " ";
                    myfile << Particle.id() << " ";
                    myfile << Particle.p(0) << " ";
                    myfile << Particle.p(1) << " ";
                    myfile << Particle.p(2) << " ";
                    myfile << Particle.p(3) << " ";
                    myfile << Particle.x(0) << " ";
                    myfile << Particle.x(1) << " ";
                    myfile << Particle.x(2) << " ";
                    myfile << Particle.x(3) << " ";
                    myfile << Particle.Nch_eta();
                    myfile << endl;
                }
            j++;
        }
    }
}

int event_selector(int N_Events, vector <myParticle> **Particles, vector<int> &selected_events, bool p_type) {

    int event = 0;
    vector<int> p_types = {}; // if true use particles if false use anti particles
    selected_events.clear(); //clear old events

    if (p_type) {
        p_types = {0, 1};
    } else {
        p_types = {2, 3};
    }

    bool usefull_event = false;

    for (int i = 0; i < N_Events; ++i) {
        //check if both are there for each event
        for (auto &proton:Particles[i][p_types[0]]) {
            for (auto &neutron:Particles[i][p_types[1]]) {
               // if (!cascade_momentum(neutron, proton)) { continue; }
                myParticle deuteron = neutron + proton;
                if (abs(deuteron.y()) <= 0.55) {
                    usefull_event = true;
                    continue;
                }
            }
        }
        for (auto &proton:Particles[i][p_types[0]]) {
            for (auto &neutron:Particles[i][p_types[0]]) {
                if(proton == neutron){ continue;}
                //if (!cascade_momentum(neutron, proton)) { continue; }
                myParticle deuteron = neutron + proton;
                if (abs(deuteron.y()) <= 0.55) {
                    usefull_event = true;
                    continue;
                }
            }
        }

        for (auto &proton:Particles[i][p_types[1]]) {
            for (auto &neutron:Particles[i][p_types[1]]) {
                if(proton == neutron){ continue;}
                //if (!cascade_momentum(neutron, proton)) { continue; }
                myParticle deuteron = neutron + proton;
                if (abs(deuteron.y()) <= 0.55) {
                    usefull_event = true;
                    continue;
                }
            }
        }
        if (usefull_event) {
            selected_events.push_back(i);
            event++;
        }
        usefull_event = false;
    }
    return event;
}


bool cascade_momentum(myParticle &neutron, myParticle &proton) {
    double delta = proton.momentum_difference(neutron);
    return delta <= 0.5;
}
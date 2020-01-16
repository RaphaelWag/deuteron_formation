//
// Created by rwagn on 27.08.2019.
//

#ifndef MYSIMULATION_ALICE_INPUT_H
#define MYSIMULATION_ALICE_INPUT_H

#include <iostream>
#include <utility>
#include <vector>
#include <sstream>

using namespace std;

class Input {
public:

    vector<string> files;
    vector<int> N_events;
    int N_simulations;
    double ff;
    bool settings_correct;
    string cms_string;
    string ALICE_data;
    string particle;
    double N_event_total;
    string dataset_folder;
    bool reduced_data;
    bool particle_type; //true for particle false for anti

    Input() = default;

    void max_simulations(int N_max) {
        if (N_simulations > N_max) { N_simulations = N_max; }
        N_event_total = 0;
        for (int i =0; i < N_max; ++i) {
            N_event_total += N_events[i];
        }
    }

    void set_particle(string particle_, bool particle_type_) {
            particle = move(particle_);
            particle_type = particle_type_;
    }

    void reduce_data() {reduced_data = true;}
    void full_data() {reduced_data = false;}

    void set_cms(double cms_init) {
        settings_correct = false;
        N_event_total = 0;

        if(reduced_data) {
            if (particle_type) {
                dataset_folder = "/mnt/d/Uni/Lectures/thesis/ALICE_datasets_pythia/np_pairs/";
            } else{
                dataset_folder = "/mnt/d/Uni/Lectures/thesis/ALICE_datasets_pythia/npbar_pairs/";
            }
        } else {
            dataset_folder = "/mnt/d/Uni/Lectures/thesis/ALICE_datasets_pythia/";
        }



        if (cms_init == 7.) {
            // settings for 7 TeV
            N_simulations = 16;
            cms_string = "7";
            files = {"min_bias_cms_7000_events_1800000_seed_1.txt",
                     "min_bias_cms_7000_events_3600000_seed_2.txt",
                     "min_bias_cms_7000_events_3000000_seed_3.txt",
                     "min_bias_cms_7000_events_4800000_seed_4.txt",
                     "min_bias_cms_7000_events_2400000_seed_5.txt",
                     "min_bias_cms_7000_events_3000000_seed_6.txt",
                     "min_bias_cms_7000_events_3000000_seed_7.txt",
                     "min_bias_cms_7000_events_1800000_seed_8.txt",
                     "min_bias_cms_7000_events_1800000_seed_9.txt",
                     "min_bias_cms_7000_events_3000000_seed_10.txt",
                     "min_bias_cms_7000_events_6000000_seed_11.txt",
                     "min_bias_cms_7000_events_8400000_seed_12.txt",
                     "min_bias_cms_7000_events_8400000_seed_13.txt",
                     "min_bias_cms_7000_events_9600000_seed_14.txt",
                     "min_bias_cms_7000_events_9600000_seed_15.txt",
                     "min_bias_cms_7000_events_8400000_seed_16.txt"};
            N_events = {1800000,3600000,3000000,4800000,2400000,
                        3000000,3000000,1800000,1800000,3000000,
                        6000000,8400000,8400000,9600000,9600000,
                        8400000};
            ff = 0.861; //fudge factor
            settings_correct = true;
        }

        if (cms_init == 2.76) {

            //settings for 2.76 TeV
            N_simulations = 4;
            cms_string = "2.76";
            files = {"min_bias_cms_2760_events_6000000_seed_1.txt",
                     "min_bias_cms_2760_events_2400000_seed_2.txt",
                     "min_bias_cms_2760_events_4800000_seed_3.txt",
                     "min_bias_cms_2760_events_4800000_seed_4.txt"};
            N_events = {6000000, 2400000,4800000,4800000};
            ff = 0.854;
            settings_correct = true;
        }


        if (cms_init == 0.9) {
            //settings for 0.9 TeV
            N_simulations = 2;
            cms_string = "0.9";
            files = {"min_bias_cms_900_events_10800000_seed_1.txt",
                     "min_bias_cms_900_events_14400000_seed_2.txt"};
            N_events = {10800000,14400000};
            ff = 0.847;
            settings_correct = true;
        }

        if (!settings_correct) {
            cout << "Wrong cms energy in input.h" << endl;
        }

        for (auto &events:N_events) {
            N_event_total += events;
        }

        ALICE_data = "/mnt/d/Uni/Lectures/thesis/ALICE_data/ALICE_"+particle+"_data_" + cms_string + ".txt";

    }
};

#endif //MYSIMULATION_ALICE_INPUT_H

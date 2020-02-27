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
    bool cs_data_;
    string cms_string;
    string ALICE_data;
    string particle;
    double N_event_total;
    string dataset_folder;
    bool reduced_data;
    bool particle_type; //true for particle false for anti

    Input() = default;

    void max_simulations(int N_max) {
        N_simulations = N_max;
        N_event_total = 0;
        for (int i = 0; i < N_max; ++i) {
            N_event_total += N_events[i];
        }
    }

    void set_particle(string particle_, bool particle_type_) {
        particle = move(particle_);
        particle_type = particle_type_;
    }

    void reduce_data() {
        reduced_data = true;
        cs_data_ = false;
    }

    void cs_data() {
        cs_data_ = true;
        reduced_data = false;
    }

    void full_data() { cs_data_ = reduced_data = false; }

    void set_cms(double cms_init) {
        settings_correct = false;
        N_event_total = 0;

        if (reduced_data) {
            if (particle_type) {
                dataset_folder = "/mnt/d/Uni/Lectures/thesis/ALICE_datasets_pythia/np_pairs/";
            } else {
                dataset_folder = "/mnt/d/Uni/Lectures/thesis/ALICE_datasets_pythia/npbar_pairs/";
            }
        } else {
            dataset_folder = "/mnt/d/Uni/Lectures/thesis/ALICE_datasets_pythia/";
        }

        if (cs_data_) {
            if (particle_type) {
                dataset_folder = "/mnt/d/Uni/Lectures/thesis/ALICE_datasets_pythia/np_pairs_cs/";
            } else {
                dataset_folder = "/mnt/d/Uni/Lectures/thesis/ALICE_datasets_pythia/npbar_pairs_cs/";
            }
        }

        if (cms_init == 7.) {
            // settings for 7 TeV
            N_simulations = 54;
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
                     "min_bias_cms_7000_events_8400000_seed_16.txt",
                     "min_bias_cms_7000_events_7200000_seed_17.txt",
                     "min_bias_cms_7000_events_7200000_seed_18.txt",
                     "min_bias_cms_7000_events_7200000_seed_19.txt",
                     "min_bias_cms_7000_events_7200000_seed_20.txt",
                     "min_bias_cms_7000_events_7200000_seed_21.txt",
                     "min_bias_cms_7000_events_7200000_seed_22.txt",
                     "min_bias_cms_7000_events_8400000_seed_23.txt",
                     "min_bias_cms_7000_events_8400000_seed_24.txt",
                     "min_bias_cms_7000_events_8400000_seed_25.txt",
                     "min_bias_cms_7000_events_8400000_seed_26.txt",
                     "min_bias_cms_7000_events_8400000_seed_27.txt",
                     "min_bias_cms_7000_events_8400000_seed_28.txt",
                     "min_bias_cms_7000_events_8400000_seed_29.txt",

                     "min_bias_cms_7000_events_8400000_seed_100.txt",
                     "min_bias_cms_7000_events_8400000_seed_101.txt",
                     "min_bias_cms_7000_events_8400000_seed_102.txt",
                     "min_bias_cms_7000_events_8400000_seed_103.txt",

                     "min_bias_cms_7000_events_8400000_seed_200.txt",
                     "min_bias_cms_7000_events_8400000_seed_202.txt",
                     "min_bias_cms_7000_events_8400000_seed_203.txt",
                     "min_bias_cms_7000_events_8400000_seed_204.txt",
                     "min_bias_cms_7000_events_8400000_seed_205.txt",

                     "min_bias_cms_7000_events_8400000_seed_301.txt",
                     "min_bias_cms_7000_events_8400000_seed_302.txt",
                     "min_bias_cms_7000_events_8400000_seed_303.txt",
                     "min_bias_cms_7000_events_8400000_seed_304.txt",
                     "min_bias_cms_7000_events_8400000_seed_305.txt",

                     "min_bias_cms_7000_events_8400000_seed_400.txt",
                     "min_bias_cms_7000_events_8400000_seed_401.txt",
                     "min_bias_cms_7000_events_8400000_seed_402.txt",
                     "min_bias_cms_7000_events_8400000_seed_403.txt",
                     "min_bias_cms_7000_events_8400000_seed_404.txt",
                     "min_bias_cms_7000_events_8400000_seed_405.txt",

                     "min_bias_cms_7000_events_8400000_seed_500.txt",
                     "min_bias_cms_7000_events_8400000_seed_501.txt",
                     "min_bias_cms_7000_events_8400000_seed_502.txt",
                     "min_bias_cms_7000_events_8400000_seed_503.txt",
                     "min_bias_cms_7000_events_8400000_seed_504.txt"

            };
            N_events = {1800000, 3600000, 3000000, 4800000, 2400000,
                        3000000, 3000000, 1800000, 1800000, 3000000,
                        6000000, 8400000, 8400000, 9600000, 9600000,
                        8400000, 7200000, 7200000, 7200000, 7200000,
                        7200000, 7200000, 8400000, 8400000, 8400000,
                        8400000, 8400000, 8400000, 8400000, 8400000,
                        8400000, 8400000, 8400000, 8400000, 8400000,
                        8400000, 8400000, 8400000, 8400000, 8400000,
                        8400000, 8400000, 8400000, 8400000, 8400000,
                        8400000, 8400000, 8400000, 8400000, 8400000,
                        8400000, 8400000, 8400000, 8400000};

            ff = 0.861; //fudge factor
            settings_correct = true;
        }

        if (cms_init == 2.76) {

            //settings for 2.76 TeV
            N_simulations = 26;
            cms_string = "2.76";
            files = {"min_bias_cms_2760_events_6000000_seed_1.txt",
                     "min_bias_cms_2760_events_2400000_seed_2.txt",
                     "min_bias_cms_2760_events_4800000_seed_3.txt",
                     "min_bias_cms_2760_events_4800000_seed_4.txt",
                     "min_bias_cms_2760_events_8400000_seed_5.txt",
                     "min_bias_cms_2760_events_8400000_seed_6.txt",
                     "min_bias_cms_2760_events_8400000_seed_7.txt",
                     "min_bias_cms_2760_events_8400000_seed_8.txt",
                     "min_bias_cms_2760_events_8400000_seed_9.txt",

                     "min_bias_cms_2760_events_8400000_seed_100.txt",
                     "min_bias_cms_2760_events_8400000_seed_101.txt",
                     "min_bias_cms_2760_events_8400000_seed_102.txt",
                     "min_bias_cms_2760_events_8400000_seed_103.txt",
                     "min_bias_cms_2760_events_8400000_seed_104.txt",

                     "min_bias_cms_2760_events_8400000_seed_200.txt",
                     "min_bias_cms_2760_events_8400000_seed_201.txt",
                     "min_bias_cms_2760_events_8400000_seed_202.txt",
                     "min_bias_cms_2760_events_8400000_seed_203.txt",
                     "min_bias_cms_2760_events_8400000_seed_204.txt",
                     "min_bias_cms_2760_events_8400000_seed_205.txt",

                     "min_bias_cms_2760_events_8400000_seed_300.txt",
                     "min_bias_cms_2760_events_8400000_seed_301.txt",
                     "min_bias_cms_2760_events_8400000_seed_302.txt",
                     "min_bias_cms_2760_events_8400000_seed_303.txt",
                     "min_bias_cms_2760_events_8400000_seed_304.txt",
                     "min_bias_cms_2760_events_8400000_seed_305.txt",
            };
            N_events = {6000000, 2400000, 4800000, 4800000, 8400000,
                        8400000, 8400000, 8400000, 8400000, 8400000,
                        8400000, 8400000, 8400000, 8400000, 8400000,
                        8400000, 8400000, 8400000, 8400000, 8400000,
                        8400000, 8400000, 8400000, 8400000, 8400000,
                        8400000};
            ff = 0.854;
            settings_correct = true;
        }


        if (cms_init == 0.9) {
            //settings for 0.9 TeV
            N_simulations = 19;
            cms_string = "0.9";
            files = {"min_bias_cms_900_events_10800000_seed_1.txt",
                     "min_bias_cms_900_events_14400000_seed_2.txt",

                     "min_bias_cms_900_events_8400000_seed_300.txt",

                     "min_bias_cms_900_events_8400000_seed_400.txt",
                     "min_bias_cms_900_events_8400000_seed_401.txt",
                     "min_bias_cms_900_events_8400000_seed_402.txt",
                     "min_bias_cms_900_events_8400000_seed_403.txt",
                     "min_bias_cms_900_events_8400000_seed_404.txt",
                     "min_bias_cms_900_events_8400000_seed_405.txt",
                     "min_bias_cms_900_events_8400000_seed_406.txt",
                     "min_bias_cms_900_events_8400000_seed_407.txt",

                     "min_bias_cms_900_events_8400000_seed_500.txt",
                     "min_bias_cms_900_events_8400000_seed_501.txt",
                     "min_bias_cms_900_events_8400000_seed_502.txt",
                     "min_bias_cms_900_events_8400000_seed_503.txt",
                     "min_bias_cms_900_events_8400000_seed_504.txt",
                     "min_bias_cms_900_events_8400000_seed_505.txt",
                     "min_bias_cms_900_events_8400000_seed_506.txt",
                     "min_bias_cms_900_events_8400000_seed_507.txt",

            };
            N_events = {10800000, 14400000, 8400000, 8400000, 8400000,
                        8400000, 8400000, 8400000, 8400000, 8400000,
                        8400000, 8400000, 8400000, 8400000, 8400000,
                        8400000, 8400000, 8400000, 8400000};
            ff = 0.847;
            settings_correct = true;
        }

        if (!settings_correct) {
            cout << "Wrong cms energy in input.h" << endl;
        }

        for (auto &events:N_events) {
            N_event_total += events;
        }

        ALICE_data = "/mnt/d/Uni/Lectures/thesis/ALICE_data/ALICE_" + particle + "_data_" + cms_string + ".txt";

    }
};

#endif //MYSIMULATION_ALICE_INPUT_H

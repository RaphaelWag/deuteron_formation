//
// Created by rwagn on 21.10.2019.
//

#include "mySimulation.h"
#include <vector>
#include "histogram2d.h"
#include "ALICE_input.h"

using namespace std;

int main() {

    Input input;
    //simulation parameters
    double cms = 7;
    double cutoff_momentum = 206./1000.;
    double cutoff_momentum_w = 215./1000.;

    histogram2d nd;

    nd.set_bins(34,34);
    nd.set_boundaries(-0.85,0.85,-0.85,0.85);
    nd.set_x_lin();
    nd.set_y_lin();
    nd.probability();

    histogram2d pd(nd);
    histogram2d nd_w(nd);
    histogram2d pd_w(nd);


    //generate input
    input.set_particle("dbar", false);
    input.reduce_data();
    input.set_cms(cms);

    auto *S = new mySimulation[input.N_simulations];

    //set weights
    for (int k = 0; k < input.N_simulations; ++k) {
        S[k].set_ALICE_weights_discrete(input.cms_string, input.particle_type);
    }

    //load files and rescale spectrum
    for (int l = 0; l < input.N_simulations; ++l) {
        S[l].load_txt(input.dataset_folder+input.files[l], input.N_events[l]);
        S[l].rescale_spectrum();
    }

    //analyse without weights
    for (int i = 0; i < input.N_simulations; ++i) {
        S[i].set_cutoff_momentum(cutoff_momentum);
        S[i].coalescence();
        for (auto &d:S[i].deuteron) {
            nd.fill(d.y_n(), d.y());
            pd.fill(d.y_p(), d.y());
        }
    }

    //analyse with weights
    for (int i = 0; i < input.N_simulations; ++i) {
        S[i].set_cutoff_momentum(cutoff_momentum_w);
        S[i].coalescence();
        for (auto &d:S[i].deuteron) {
            nd_w.fill(d.y_n(), d.y(), d.w());
            pd_w.fill(d.y_p(), d.y(), d.w());
        }
    }

    //print results
    nd.print("/mnt/d/Uni/Lectures/thesis/ALICE_results/y_kinematics/y_nd.txt");
    nd_w.print("/mnt/d/Uni/Lectures/thesis/ALICE_results/y_kinematics/y_nd_w.txt");
    pd.print("/mnt/d/Uni/Lectures/thesis/ALICE_results/y_kinematics/y_pd.txt");
    pd_w.print("/mnt/d/Uni/Lectures/thesis/ALICE_results/y_kinematics/y_pd_w.txt");



    return 0;
}

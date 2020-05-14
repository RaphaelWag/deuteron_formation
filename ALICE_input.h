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
    double nsdtoinel;
    double norm_err_above;
    double norm_err_below;
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

        if (cms_init == 13.) {
            // settings for 13 TeV
            N_simulations = 324;
            cms_string = "13";
            ff = 0.865; //fudge factor
            nsdtoinel = 0.836;
            norm_err_above = 0.05;
            norm_err_below = 0.02;
            files = {
                    "min_bias_cms_13000_events_4800000_seed_1.txt",
                    "min_bias_cms_13000_events_2400000_seed_2.txt",
                    "min_bias_cms_13000_events_2400000_seed_3.txt",
                    "min_bias_cms_13000_events_2400000_seed_4.txt",
                    "min_bias_cms_13000_events_2400000_seed_5.txt",
                    "min_bias_cms_13000_events_2400000_seed_6.txt",
                    "min_bias_cms_13000_events_2400000_seed_7.txt",
                    "min_bias_cms_13000_events_2400000_seed_8.txt",
                    "min_bias_cms_13000_events_2400000_seed_9.txt",
                    "min_bias_cms_13000_events_2400000_seed_10.txt",
                    "min_bias_cms_13000_events_2400000_seed_11.txt",
                    "min_bias_cms_13000_events_2400000_seed_12.txt",
                    "min_bias_cms_13000_events_2400000_seed_13.txt",
                    "min_bias_cms_13000_events_2400000_seed_14.txt",
                    "min_bias_cms_13000_events_2400000_seed_15.txt",
                    "min_bias_cms_13000_events_2400000_seed_16.txt",
                    "min_bias_cms_13000_events_2400000_seed_17.txt",
                    "min_bias_cms_13000_events_2400000_seed_18.txt",
                    "min_bias_cms_13000_events_2400000_seed_19.txt",
                    "min_bias_cms_13000_events_2400000_seed_20.txt",
                    "min_bias_cms_13000_events_2400000_seed_21.txt",
                    "min_bias_cms_13000_events_2400000_seed_22.txt",
                    "min_bias_cms_13000_events_2400000_seed_23.txt",
                    "min_bias_cms_13000_events_2400000_seed_24.txt",
                    "min_bias_cms_13000_events_2400000_seed_25.txt",
                    "min_bias_cms_13000_events_2400000_seed_26.txt",
                    "min_bias_cms_13000_events_2400000_seed_27.txt",
                    "min_bias_cms_13000_events_2400000_seed_28.txt",
                    "min_bias_cms_13000_events_2400000_seed_29.txt",
                    "min_bias_cms_13000_events_2400000_seed_30.txt",
                    "min_bias_cms_13000_events_2400000_seed_31.txt",
                    "min_bias_cms_13000_events_2400000_seed_32.txt",
                    "min_bias_cms_13000_events_2400000_seed_33.txt",

                    "min_bias_cms_13000_events_2400000_seed_100.txt",
                    "min_bias_cms_13000_events_2400000_seed_101.txt",
                    "min_bias_cms_13000_events_2400000_seed_102.txt",
                    "min_bias_cms_13000_events_2400000_seed_103.txt",
                    "min_bias_cms_13000_events_2400000_seed_104.txt",
                    "min_bias_cms_13000_events_2400000_seed_105.txt",
                    "min_bias_cms_13000_events_2400000_seed_106.txt",
                    "min_bias_cms_13000_events_2400000_seed_107.txt",
                    "min_bias_cms_13000_events_2400000_seed_108.txt",
                    "min_bias_cms_13000_events_2400000_seed_109.txt",
                    "min_bias_cms_13000_events_2400000_seed_110.txt",
                    "min_bias_cms_13000_events_2400000_seed_111.txt",
                    "min_bias_cms_13000_events_2400000_seed_112.txt",
                    "min_bias_cms_13000_events_2400000_seed_113.txt",
                    "min_bias_cms_13000_events_2400000_seed_114.txt",
                    "min_bias_cms_13000_events_2400000_seed_115.txt",
                    "min_bias_cms_13000_events_2400000_seed_116.txt",
                    "min_bias_cms_13000_events_2400000_seed_117.txt",
                    "min_bias_cms_13000_events_2400000_seed_118.txt",
                    "min_bias_cms_13000_events_2400000_seed_119.txt",
                    "min_bias_cms_13000_events_2400000_seed_120.txt",
                    "min_bias_cms_13000_events_2400000_seed_121.txt",
                    "min_bias_cms_13000_events_2400000_seed_122.txt",
                    "min_bias_cms_13000_events_2400000_seed_123.txt",
                    "min_bias_cms_13000_events_2400000_seed_124.txt",
                    "min_bias_cms_13000_events_2400000_seed_125.txt",
                    "min_bias_cms_13000_events_2400000_seed_126.txt",
                    "min_bias_cms_13000_events_2400000_seed_127.txt",
                    "min_bias_cms_13000_events_2400000_seed_128.txt",
                    "min_bias_cms_13000_events_2400000_seed_129.txt",
                    "min_bias_cms_13000_events_2400000_seed_130.txt",
                    "min_bias_cms_13000_events_2400000_seed_131.txt",
                    "min_bias_cms_13000_events_2400000_seed_132.txt",

                    "min_bias_cms_13000_events_4800000_seed_200.txt",
                    "min_bias_cms_13000_events_2400000_seed_201.txt",
                    "min_bias_cms_13000_events_2400000_seed_202.txt",
                    "min_bias_cms_13000_events_2400000_seed_203.txt",
                    "min_bias_cms_13000_events_2400000_seed_204.txt",
                    "min_bias_cms_13000_events_2400000_seed_205.txt",
                    "min_bias_cms_13000_events_2400000_seed_206.txt",
                    "min_bias_cms_13000_events_2400000_seed_207.txt",
                    "min_bias_cms_13000_events_2400000_seed_208.txt",
                    "min_bias_cms_13000_events_2400000_seed_209.txt",
                    "min_bias_cms_13000_events_2400000_seed_210.txt",
                    "min_bias_cms_13000_events_2400000_seed_211.txt",
                    "min_bias_cms_13000_events_2400000_seed_212.txt",
                    "min_bias_cms_13000_events_2400000_seed_213.txt",
                    "min_bias_cms_13000_events_2400000_seed_214.txt",
                    "min_bias_cms_13000_events_2400000_seed_215.txt",
                    "min_bias_cms_13000_events_2400000_seed_216.txt",
                    "min_bias_cms_13000_events_2400000_seed_217.txt",
                    "min_bias_cms_13000_events_2400000_seed_218.txt",
                    "min_bias_cms_13000_events_2400000_seed_219.txt",
                    "min_bias_cms_13000_events_2400000_seed_220.txt",
                    "min_bias_cms_13000_events_2400000_seed_221.txt",
                    "min_bias_cms_13000_events_2400000_seed_222.txt",
                    "min_bias_cms_13000_events_2400000_seed_223.txt",
                    "min_bias_cms_13000_events_2400000_seed_224.txt",
                    "min_bias_cms_13000_events_2400000_seed_225.txt",
                    "min_bias_cms_13000_events_2400000_seed_226.txt",
                    "min_bias_cms_13000_events_2400000_seed_227.txt",
                    "min_bias_cms_13000_events_2400000_seed_228.txt",
                    "min_bias_cms_13000_events_2400000_seed_229.txt",
                    "min_bias_cms_13000_events_2400000_seed_230.txt",
                    "min_bias_cms_13000_events_2400000_seed_231.txt",
                    "min_bias_cms_13000_events_2400000_seed_232.txt",

                    "min_bias_cms_13000_events_4800000_seed_300.txt",
                    "min_bias_cms_13000_events_2400000_seed_301.txt",
                    "min_bias_cms_13000_events_2400000_seed_302.txt",
                    "min_bias_cms_13000_events_2400000_seed_303.txt",
                    "min_bias_cms_13000_events_2400000_seed_304.txt",
                    "min_bias_cms_13000_events_2400000_seed_305.txt",
                    "min_bias_cms_13000_events_2400000_seed_306.txt",
                    "min_bias_cms_13000_events_2400000_seed_307.txt",
                    "min_bias_cms_13000_events_2400000_seed_308.txt",
                    "min_bias_cms_13000_events_2400000_seed_309.txt",
                    "min_bias_cms_13000_events_2400000_seed_310.txt",
                    "min_bias_cms_13000_events_2400000_seed_311.txt",
                    "min_bias_cms_13000_events_2400000_seed_312.txt",
                    "min_bias_cms_13000_events_2400000_seed_313.txt",
                    //"min_bias_cms_13000_events_2400000_seed_314.txt", //bugged
                    "min_bias_cms_13000_events_2400000_seed_315.txt",
                    "min_bias_cms_13000_events_2400000_seed_316.txt",
                    "min_bias_cms_13000_events_2400000_seed_317.txt",
                    "min_bias_cms_13000_events_2400000_seed_318.txt",
                    "min_bias_cms_13000_events_2400000_seed_319.txt",
                    "min_bias_cms_13000_events_2400000_seed_320.txt",
                    "min_bias_cms_13000_events_2400000_seed_321.txt",
                    "min_bias_cms_13000_events_2400000_seed_322.txt",
                    "min_bias_cms_13000_events_2400000_seed_323.txt",
                    "min_bias_cms_13000_events_2400000_seed_324.txt",
                    "min_bias_cms_13000_events_2400000_seed_325.txt",
                    "min_bias_cms_13000_events_2400000_seed_326.txt",
                    "min_bias_cms_13000_events_2400000_seed_327.txt",
                    "min_bias_cms_13000_events_2400000_seed_328.txt",
                    "min_bias_cms_13000_events_2400000_seed_329.txt",
                    "min_bias_cms_13000_events_2400000_seed_330.txt",
                    "min_bias_cms_13000_events_2400000_seed_331.txt",
                    "min_bias_cms_13000_events_2400000_seed_332.txt",

                    "min_bias_cms_13000_events_4800000_seed_400.txt",
                    "min_bias_cms_13000_events_2400000_seed_401.txt",
                    "min_bias_cms_13000_events_2400000_seed_402.txt",
                    "min_bias_cms_13000_events_2400000_seed_403.txt",
                    "min_bias_cms_13000_events_2400000_seed_404.txt",
                    "min_bias_cms_13000_events_2400000_seed_405.txt",
                    "min_bias_cms_13000_events_2400000_seed_406.txt",
                    "min_bias_cms_13000_events_2400000_seed_407.txt",
                    "min_bias_cms_13000_events_2400000_seed_408.txt",
                    "min_bias_cms_13000_events_2400000_seed_409.txt",
                    "min_bias_cms_13000_events_2400000_seed_410.txt",
                    "min_bias_cms_13000_events_2400000_seed_411.txt",
                    "min_bias_cms_13000_events_2400000_seed_412.txt",
                    "min_bias_cms_13000_events_2400000_seed_413.txt",
                    "min_bias_cms_13000_events_2400000_seed_414.txt",
                    "min_bias_cms_13000_events_2400000_seed_415.txt",
                    "min_bias_cms_13000_events_2400000_seed_416.txt",
                    "min_bias_cms_13000_events_2400000_seed_417.txt",
                    "min_bias_cms_13000_events_2400000_seed_418.txt",
                    "min_bias_cms_13000_events_2400000_seed_419.txt",
                    "min_bias_cms_13000_events_2400000_seed_420.txt",
                    "min_bias_cms_13000_events_2400000_seed_421.txt",
                    "min_bias_cms_13000_events_2400000_seed_422.txt",
                    "min_bias_cms_13000_events_2400000_seed_423.txt",
                    "min_bias_cms_13000_events_2400000_seed_424.txt",
                    "min_bias_cms_13000_events_2400000_seed_425.txt",
                    "min_bias_cms_13000_events_2400000_seed_426.txt",
                    "min_bias_cms_13000_events_2400000_seed_427.txt",
                    "min_bias_cms_13000_events_2400000_seed_428.txt",
                    "min_bias_cms_13000_events_2400000_seed_429.txt",
                    "min_bias_cms_13000_events_2400000_seed_430.txt",
                    "min_bias_cms_13000_events_2400000_seed_431.txt",
                    "min_bias_cms_13000_events_2400000_seed_432.txt",

                    "min_bias_cms_13000_events_2400000_seed_500.txt",
                    "min_bias_cms_13000_events_2400000_seed_501.txt",
                    "min_bias_cms_13000_events_2400000_seed_502.txt",
                    "min_bias_cms_13000_events_2400000_seed_503.txt",
                    "min_bias_cms_13000_events_2400000_seed_504.txt",
                    "min_bias_cms_13000_events_2400000_seed_505.txt",
                    "min_bias_cms_13000_events_2400000_seed_506.txt",
                    "min_bias_cms_13000_events_2400000_seed_507.txt",
                    "min_bias_cms_13000_events_2400000_seed_508.txt",
                    "min_bias_cms_13000_events_2400000_seed_509.txt",
                    "min_bias_cms_13000_events_2400000_seed_510.txt",
                    "min_bias_cms_13000_events_2400000_seed_511.txt",
                    "min_bias_cms_13000_events_2400000_seed_512.txt",
                    "min_bias_cms_13000_events_2400000_seed_513.txt",
                    "min_bias_cms_13000_events_2400000_seed_514.txt",
                    "min_bias_cms_13000_events_2400000_seed_515.txt",
                    "min_bias_cms_13000_events_2400000_seed_516.txt",
                    "min_bias_cms_13000_events_2400000_seed_517.txt",
                    "min_bias_cms_13000_events_2400000_seed_518.txt",
                    "min_bias_cms_13000_events_2400000_seed_519.txt",
                    "min_bias_cms_13000_events_2400000_seed_520.txt",
                    "min_bias_cms_13000_events_2400000_seed_521.txt",
                    "min_bias_cms_13000_events_2400000_seed_522.txt",
                    "min_bias_cms_13000_events_2400000_seed_523.txt",
                    "min_bias_cms_13000_events_2400000_seed_524.txt",
                    "min_bias_cms_13000_events_2400000_seed_525.txt",
                    "min_bias_cms_13000_events_2400000_seed_526.txt",
                    "min_bias_cms_13000_events_2400000_seed_527.txt",
                    "min_bias_cms_13000_events_2400000_seed_528.txt",
                    "min_bias_cms_13000_events_2400000_seed_529.txt",
                    "min_bias_cms_13000_events_2400000_seed_530.txt",
                    "min_bias_cms_13000_events_2400000_seed_531.txt",
                    "min_bias_cms_13000_events_2400000_seed_532.txt",

                    "min_bias_cms_13000_events_2400000_seed_600.txt",
                    "min_bias_cms_13000_events_2400000_seed_601.txt",
                    "min_bias_cms_13000_events_2400000_seed_602.txt",
                    "min_bias_cms_13000_events_2400000_seed_603.txt",
                    "min_bias_cms_13000_events_2400000_seed_604.txt",
                    "min_bias_cms_13000_events_2400000_seed_605.txt",
                    "min_bias_cms_13000_events_2400000_seed_606.txt",
                    "min_bias_cms_13000_events_2400000_seed_607.txt",
                    "min_bias_cms_13000_events_2400000_seed_608.txt",
                    "min_bias_cms_13000_events_2400000_seed_609.txt",
                    "min_bias_cms_13000_events_2400000_seed_610.txt",
                    "min_bias_cms_13000_events_2400000_seed_611.txt",
                    "min_bias_cms_13000_events_2400000_seed_612.txt",
                    "min_bias_cms_13000_events_2400000_seed_613.txt",
                    "min_bias_cms_13000_events_2400000_seed_614.txt",
                    "min_bias_cms_13000_events_2400000_seed_615.txt",
                    "min_bias_cms_13000_events_2400000_seed_616.txt",
                    "min_bias_cms_13000_events_2400000_seed_617.txt",
                    "min_bias_cms_13000_events_2400000_seed_618.txt",
                    "min_bias_cms_13000_events_2400000_seed_619.txt",
                    "min_bias_cms_13000_events_2400000_seed_620.txt",
                    "min_bias_cms_13000_events_2400000_seed_621.txt",
                    "min_bias_cms_13000_events_2400000_seed_622.txt",
                    "min_bias_cms_13000_events_2400000_seed_623.txt",
                    "min_bias_cms_13000_events_2400000_seed_624.txt",
                    "min_bias_cms_13000_events_2400000_seed_625.txt",
                    "min_bias_cms_13000_events_2400000_seed_626.txt",
                    "min_bias_cms_13000_events_2400000_seed_627.txt",
                    "min_bias_cms_13000_events_2400000_seed_628.txt",
                    "min_bias_cms_13000_events_2400000_seed_629.txt",
                    "min_bias_cms_13000_events_2400000_seed_630.txt",
                    "min_bias_cms_13000_events_2400000_seed_631.txt",
                    "min_bias_cms_13000_events_2400000_seed_632.txt",

                    "min_bias_cms_13000_events_2400000_seed_700.txt",
                    "min_bias_cms_13000_events_2400000_seed_701.txt",
                    "min_bias_cms_13000_events_2400000_seed_702.txt",
                    "min_bias_cms_13000_events_2400000_seed_703.txt",
                    "min_bias_cms_13000_events_2400000_seed_704.txt",
                    "min_bias_cms_13000_events_2400000_seed_705.txt",
                    "min_bias_cms_13000_events_2400000_seed_706.txt",
                    "min_bias_cms_13000_events_2400000_seed_707.txt",
                    "min_bias_cms_13000_events_2400000_seed_708.txt",
                    "min_bias_cms_13000_events_2400000_seed_709.txt",
                    "min_bias_cms_13000_events_2400000_seed_710.txt",
                    "min_bias_cms_13000_events_2400000_seed_711.txt",
                    "min_bias_cms_13000_events_2400000_seed_712.txt",
                    "min_bias_cms_13000_events_2400000_seed_713.txt",
                    "min_bias_cms_13000_events_2400000_seed_714.txt",
                    "min_bias_cms_13000_events_2400000_seed_715.txt",
                    "min_bias_cms_13000_events_2400000_seed_716.txt",
                    "min_bias_cms_13000_events_2400000_seed_717.txt",
                    "min_bias_cms_13000_events_2400000_seed_718.txt",
                    "min_bias_cms_13000_events_2400000_seed_719.txt",
                    "min_bias_cms_13000_events_2400000_seed_720.txt",
                    "min_bias_cms_13000_events_2400000_seed_721.txt",
                    "min_bias_cms_13000_events_2400000_seed_722.txt",
                    "min_bias_cms_13000_events_2400000_seed_723.txt",
                    "min_bias_cms_13000_events_2400000_seed_724.txt",
                    "min_bias_cms_13000_events_2400000_seed_725.txt",
                    "min_bias_cms_13000_events_2400000_seed_726.txt",
                    "min_bias_cms_13000_events_2400000_seed_727.txt",
                    "min_bias_cms_13000_events_2400000_seed_728.txt",
                    "min_bias_cms_13000_events_2400000_seed_729.txt",
                    "min_bias_cms_13000_events_2400000_seed_730.txt",
                    "min_bias_cms_13000_events_2400000_seed_731.txt",
                    "min_bias_cms_13000_events_2400000_seed_732.txt",

                    //"min_bias_cms_13000_events_2400000_seed_800.txt",//bugged
                    "min_bias_cms_13000_events_2400000_seed_801.txt",
                    "min_bias_cms_13000_events_2400000_seed_802.txt",
                    "min_bias_cms_13000_events_2400000_seed_803.txt",
                    "min_bias_cms_13000_events_2400000_seed_804.txt",
                    "min_bias_cms_13000_events_2400000_seed_805.txt",
                    "min_bias_cms_13000_events_2400000_seed_806.txt",
                    "min_bias_cms_13000_events_2400000_seed_807.txt",
                    "min_bias_cms_13000_events_2400000_seed_808.txt",
                    "min_bias_cms_13000_events_2400000_seed_809.txt",
                    "min_bias_cms_13000_events_2400000_seed_810.txt",
                    "min_bias_cms_13000_events_2400000_seed_811.txt",
                    "min_bias_cms_13000_events_2400000_seed_812.txt",
                    "min_bias_cms_13000_events_2400000_seed_813.txt",
                    "min_bias_cms_13000_events_2400000_seed_814.txt",
                    "min_bias_cms_13000_events_2400000_seed_815.txt",
                    "min_bias_cms_13000_events_2400000_seed_816.txt",
                    "min_bias_cms_13000_events_2400000_seed_817.txt",
                    "min_bias_cms_13000_events_2400000_seed_818.txt",
                    "min_bias_cms_13000_events_2400000_seed_819.txt",
                    "min_bias_cms_13000_events_2400000_seed_820.txt",
                    "min_bias_cms_13000_events_2400000_seed_821.txt",
                    "min_bias_cms_13000_events_2400000_seed_822.txt",
                    "min_bias_cms_13000_events_2400000_seed_823.txt",
                    "min_bias_cms_13000_events_2400000_seed_824.txt",
                    "min_bias_cms_13000_events_2400000_seed_825.txt",
                    "min_bias_cms_13000_events_2400000_seed_826.txt",
                    "min_bias_cms_13000_events_2400000_seed_827.txt",
                    "min_bias_cms_13000_events_2400000_seed_828.txt",
                    "min_bias_cms_13000_events_2400000_seed_829.txt",
                    "min_bias_cms_13000_events_2400000_seed_830.txt",

                    "min_bias_cms_13000_events_2400000_seed_900.txt",
                    "min_bias_cms_13000_events_2400000_seed_901.txt",
                    "min_bias_cms_13000_events_2400000_seed_902.txt",
                    "min_bias_cms_13000_events_2400000_seed_903.txt",
                    "min_bias_cms_13000_events_2400000_seed_904.txt",
                    "min_bias_cms_13000_events_2400000_seed_905.txt",
                    "min_bias_cms_13000_events_2400000_seed_906.txt",
                    "min_bias_cms_13000_events_2400000_seed_907.txt",
                    "min_bias_cms_13000_events_2400000_seed_908.txt",
                    "min_bias_cms_13000_events_2400000_seed_909.txt",
                    "min_bias_cms_13000_events_2400000_seed_910.txt",
                    "min_bias_cms_13000_events_2400000_seed_911.txt",
                    "min_bias_cms_13000_events_2400000_seed_912.txt",
                    "min_bias_cms_13000_events_2400000_seed_913.txt",
                    "min_bias_cms_13000_events_2400000_seed_914.txt",
                    "min_bias_cms_13000_events_2400000_seed_915.txt",
                    "min_bias_cms_13000_events_2400000_seed_916.txt",
                    "min_bias_cms_13000_events_2400000_seed_917.txt",
                    "min_bias_cms_13000_events_2400000_seed_918.txt",
                    "min_bias_cms_13000_events_2400000_seed_919.txt",
                    "min_bias_cms_13000_events_2400000_seed_920.txt",
                    "min_bias_cms_13000_events_2400000_seed_921.txt",
                    "min_bias_cms_13000_events_2400000_seed_922.txt",
                    "min_bias_cms_13000_events_2400000_seed_923.txt",
                    "min_bias_cms_13000_events_2400000_seed_924.txt",
                    "min_bias_cms_13000_events_2400000_seed_925.txt",
                    "min_bias_cms_13000_events_2400000_seed_926.txt",
                    "min_bias_cms_13000_events_2400000_seed_927.txt",
                    "min_bias_cms_13000_events_2400000_seed_928.txt",
                    "min_bias_cms_13000_events_2400000_seed_929.txt",
                    "min_bias_cms_13000_events_2400000_seed_930.txt"


            };
            N_events = {
                    //0
                    4800000, 2400000, 2400000, 2400000, 2400000, 2400000, 2400000, 2400000, 2400000, 2400000,
                    2400000, 2400000, 2400000, 2400000, 2400000, 2400000, 2400000, 2400000, 2400000, 2400000,
                    2400000, 2400000, 2400000, 2400000, 2400000, 2400000, 2400000, 2400000, 2400000, 2400000,
                    2400000, 2400000, 2400000,

                    //1
                    2400000, 2400000, 2400000, 2400000, 2400000, 2400000, 2400000, 2400000, 2400000, 2400000,
                    2400000, 2400000, 2400000, 2400000, 2400000, 2400000, 2400000, 2400000, 2400000, 2400000,
                    2400000, 2400000, 2400000, 2400000, 2400000, 2400000, 2400000, 2400000, 2400000, 2400000,
                    2400000, 2400000, 2400000,

                    //2
                    4800000, 2400000, 2400000, 2400000, 2400000, 2400000, 2400000, 2400000, 2400000, 2400000,
                    2400000, 2400000, 2400000, 2400000, 2400000, 2400000, 2400000, 2400000, 2400000, 2400000,
                    2400000, 2400000, 2400000, 2400000, 2400000, 2400000, 2400000, 2400000, 2400000, 2400000,
                    2400000, 2400000, 2400000,

                    //3
                    4800000, 2400000, 2400000, 2400000, 2400000, 2400000, 2400000, 2400000, 2400000, 2400000,
                    2400000, 2400000, 2400000, 2400000, /*2400000,*/ 2400000, 2400000, 2400000, 2400000, 2400000,
                    2400000, 2400000, 2400000, 2400000, 2400000, 2400000, 2400000, 2400000, 2400000, 2400000,
                    2400000, 2400000, 2400000,

                    //4
                    4800000, 2400000, 2400000, 2400000, 2400000, 2400000, 2400000, 2400000, 2400000, 2400000,
                    2400000, 2400000, 2400000, 2400000, 2400000, 2400000, 2400000, 2400000, 2400000, 2400000,
                    2400000, 2400000, 2400000, 2400000, 2400000, 2400000, 2400000, 2400000, 2400000, 2400000,
                    2400000, 2400000, 2400000,

                    //5
                    2400000, 2400000, 2400000, 2400000, 2400000, 2400000, 2400000, 2400000, 2400000, 2400000,
                    2400000, 2400000, 2400000, 2400000, 2400000, 2400000, 2400000, 2400000, 2400000, 2400000,
                    2400000, 2400000, 2400000, 2400000, 2400000, 2400000, 2400000, 2400000, 2400000, 2400000,
                    2400000, 2400000, 2400000,

                    //6
                    2400000, 2400000, 2400000, 2400000, 2400000, 2400000, 2400000, 2400000, 2400000, 2400000,
                    2400000, 2400000, 2400000, 2400000, 2400000, 2400000, 2400000, 2400000, 2400000, 2400000,
                    2400000, 2400000, 2400000, 2400000, 2400000, 2400000, 2400000, 2400000, 2400000, 2400000,
                    2400000, 2400000, 2400000,

                    //7
                    2400000, 2400000, 2400000, 2400000, 2400000, 2400000, 2400000, 2400000, 2400000, 2400000,
                    2400000, 2400000, 2400000, 2400000, 2400000, 2400000, 2400000, 2400000, 2400000, 2400000,
                    2400000, 2400000, 2400000, 2400000, 2400000, 2400000, 2400000, 2400000, 2400000, 2400000,
                    2400000, 2400000, 2400000,

                    //8
                    /*2400000,*/ 2400000, 2400000, 2400000, 2400000, 2400000, 2400000, 2400000, 2400000, 2400000,
                    2400000, 2400000, 2400000, 2400000, 2400000, 2400000, 2400000, 2400000, 2400000, 2400000,
                    2400000, 2400000, 2400000, 2400000, 2400000, 2400000, 2400000, 2400000, 2400000, 2400000,
                    2400000,
                    //9
                    2400000, 2400000, 2400000, 2400000, 2400000, 2400000, 2400000, 2400000, 2400000, 2400000,
                    2400000, 2400000, 2400000, 2400000, 2400000, 2400000, 2400000, 2400000, 2400000, 2400000,
                    2400000, 2400000, 2400000, 2400000, 2400000, 2400000, 2400000, 2400000, 2400000, 2400000,
                    2400000
            };


            settings_correct = true;
        }
        if (cms_init == 7.) {
            // settings for 7 TeV
            N_simulations = 54;
            cms_string = "7";
            ff = 0.861; //fudge factor
            nsdtoinel = 0.742;
            norm_err_above = 0.05;
            norm_err_below = 0.02;
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


            settings_correct = true;
        }

        if (cms_init == 2.76) {

            //settings for 2.76 TeV
            N_simulations = 31;
            cms_string = "2.76";
            ff = 0.854;
            nsdtoinel = 0.760;
            norm_err_above = 0.052;
            norm_err_below = 0.028;
            files = {"min_bias_cms_2760_events_6000000_seed_1.txt",
                     "min_bias_cms_2760_events_2400000_seed_2.txt",
                     "min_bias_cms_2760_events_4800000_seed_3.txt",
                     "min_bias_cms_2760_events_4800000_seed_4.txt",
                     "min_bias_cms_2760_events_8400000_seed_5.txt",
                     "min_bias_cms_2760_events_8400000_seed_6.txt",
                     "min_bias_cms_2760_events_8400000_seed_7.txt",
                     "min_bias_cms_2760_events_8400000_seed_8.txt",
                     "min_bias_cms_2760_events_8400000_seed_9.txt",
                     "min_bias_cms_2760_events_8400000_seed_10.txt",
                     "min_bias_cms_2760_events_8400000_seed_11.txt",

                     "min_bias_cms_2760_events_8400000_seed_100.txt",
                     "min_bias_cms_2760_events_8400000_seed_101.txt",
                     "min_bias_cms_2760_events_8400000_seed_102.txt",
                     "min_bias_cms_2760_events_8400000_seed_103.txt",
                     "min_bias_cms_2760_events_8400000_seed_104.txt",
                     "min_bias_cms_2760_events_8400000_seed_105.txt",
                     "min_bias_cms_2760_events_8400000_seed_106.txt",
                     "min_bias_cms_2760_events_8400000_seed_107.txt",

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
                        8400000, 8400000, 8400000, 8400000, 8400000,
                        8400000};

            settings_correct = true;
        }


        if (cms_init == 0.9) {
            //settings for 0.9 TeV
            N_simulations = 37;
            cms_string = "0.9";
            ff = 0.847;
            nsdtoinel = 0.763;
            norm_err_above = 0.022;
            norm_err_below = 0.008;
            files = {"min_bias_cms_900_events_10800000_seed_1.txt",
                     "min_bias_cms_900_events_14400000_seed_2.txt",

                     "min_bias_cms_900_events_8400000_seed_200.txt",
                     "min_bias_cms_900_events_8400000_seed_201.txt",
                     "min_bias_cms_900_events_8400000_seed_202.txt",
                     "min_bias_cms_900_events_8400000_seed_203.txt",
                     "min_bias_cms_900_events_8400000_seed_204.txt",

                     "min_bias_cms_900_events_8400000_seed_300.txt",
                     "min_bias_cms_900_events_8400000_seed_301.txt",
                     "min_bias_cms_900_events_8400000_seed_302.txt",
                     "min_bias_cms_900_events_8400000_seed_303.txt",
                     "min_bias_cms_900_events_8400000_seed_304.txt",
                     "min_bias_cms_900_events_8400000_seed_305.txt",

                     "min_bias_cms_900_events_8400000_seed_400.txt",
                     "min_bias_cms_900_events_8400000_seed_401.txt",
                     "min_bias_cms_900_events_8400000_seed_402.txt",
                     "min_bias_cms_900_events_8400000_seed_403.txt",
                     "min_bias_cms_900_events_8400000_seed_404.txt",
                     "min_bias_cms_900_events_8400000_seed_405.txt",
                     "min_bias_cms_900_events_8400000_seed_406.txt",
                     "min_bias_cms_900_events_8400000_seed_407.txt",
                     "min_bias_cms_900_events_8400000_seed_408.txt",
                     "min_bias_cms_900_events_8400000_seed_409.txt",
                     "min_bias_cms_900_events_8400000_seed_410.txt",
                     "min_bias_cms_900_events_8400000_seed_411.txt",

                     "min_bias_cms_900_events_8400000_seed_500.txt",
                     "min_bias_cms_900_events_8400000_seed_501.txt",
                     "min_bias_cms_900_events_8400000_seed_502.txt",
                     "min_bias_cms_900_events_8400000_seed_503.txt",
                     "min_bias_cms_900_events_8400000_seed_504.txt",
                     "min_bias_cms_900_events_8400000_seed_505.txt",
                     "min_bias_cms_900_events_8400000_seed_506.txt",
                     "min_bias_cms_900_events_8400000_seed_507.txt",
                     "min_bias_cms_900_events_8400000_seed_508.txt",
                     "min_bias_cms_900_events_8400000_seed_509.txt",
                     "min_bias_cms_900_events_8400000_seed_510.txt",
                     "min_bias_cms_900_events_8400000_seed_511.txt",

            };
            N_events = {10800000, 14400000, 8400000, 8400000, 8400000,
                        8400000, 8400000, 8400000, 8400000, 8400000,
                        8400000, 8400000, 8400000, 8400000, 8400000,
                        8400000, 8400000, 8400000, 8400000, 8400000,
                        8400000, 8400000, 8400000, 8400000, 8400000,
                        8400000, 8400000, 8400000, 8400000, 8400000,
                        8400000, 8400000, 8400000, 8400000, 8400000,
                        8400000, 8400000};

            settings_correct = true;
        }

        if (!settings_correct) {
            cout << "Wrong cms energy in input.h" << endl;
        }

        for (auto &events:N_events) {
            N_event_total += events;
        }

        ALICE_data = "/mnt/d/Uni/Lectures/thesis/ALICE_data_new/ALICE_" + particle + "_data_" + cms_string + ".txt";

    }
};

#endif //MYSIMULATION_ALICE_INPUT_H

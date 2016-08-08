//
//  main.cpp
//  speed_test
//
//  Created by Ellis Sparky Hoag on 6/13/16.
//  Copyright © 2016 Ellis Sparky Hoag. All rights reserved.
//

/*
 *  Records as of 7/30/16
 *
 *  1000 sites
 *  One thread: 0.0110322
 *  Two threads: 0.00808198
 *  Four threads: 0.00971464
 *
 *  2000 sites
 *  One thread: 0.0231604
 *  Two threads: 0.0161947
 *  Four threads: 0.0200929
 *
 *  4000 sites:
 *  One thread: 0.0630737
 *  Two threads: 0.0433619
 *  Four threads: 0.0457458
 *
 *  8000 sites:
 *  One thread: 0.163324
 *  Two threads: 0.105508
 *  Four threads: 0.115384
 *
 *  16000 sites:
 *  One thread: 0.45876
 *  Two threads: 0.284353
 *  Four threads: 0.311399
 */

#include <iostream>
#include "voronoi_sphere.h"

using namespace std;
using namespace Voronoi;

const int num_sites = 4000;
const int num_trials = 100;

int main(int argc, const char * argv[]) {
    
    srand((unsigned int)time(NULL));
    
    chrono::duration<float> total_time_one_thread;
    chrono::duration<float> total_time_two_threads;
    chrono::duration<float> total_time_four_threads;
    
    chrono::duration<float> max_trial_time_one_thread = chrono::duration<float>(0);
    chrono::duration<float> max_trial_time_two_threads = chrono::duration<float>(0);
    chrono::duration<float> max_trial_time_four_threads = chrono::duration<float>(0);
    
    chrono::time_point<chrono::system_clock> start_time;
    chrono::duration<float> trial_time;
    
    for (int trial = 1; trial <= num_trials; trial++)
    {
        unsigned int seed = rand();
        
        srand(seed);
        
        cout << trial << ") Seed = " << seed << endl;
        
        vector<tuple<double, double, double>> verts;
        
        const int precision = numeric_limits<int>::max();
        for (int i = 0; i < num_sites; i++)
        {
            double x = (rand() % precision) / (double)precision - 0.5f;
            double y = (rand() % precision) / (double)precision - 0.5f;
            double z = (rand() % precision) / (double)precision - 0.5f;
            
            double r = sqrt(x*x + y*y + z*z);
            
            tuple<double, double, double> point = make_tuple(x / r, y / r, z / r);
            
            verts.push_back(point);
        }

        start_time = chrono::system_clock::now();
        
        generate_voronoi(&verts, ONE_THREAD);
        
        trial_time = chrono::system_clock::now() - start_time;
        
        total_time_one_thread += trial_time;
        
        if (trial_time > max_trial_time_one_thread)
        {
            max_trial_time_one_thread = trial_time;
        }

        start_time = chrono::system_clock::now();
        
        generate_voronoi(&verts, TWO_THREADS);
        
        trial_time = chrono::system_clock::now() - start_time;
        
        total_time_two_threads += trial_time;
        
        if (trial_time > max_trial_time_two_threads)
        {
            max_trial_time_two_threads = trial_time;
        }

        start_time = chrono::system_clock::now();
        
        generate_voronoi(&verts, FOUR_THREADS);
        
        trial_time = chrono::system_clock::now() - start_time;
        
        total_time_four_threads += trial_time;
        
        if (trial_time > max_trial_time_four_threads)
        {
            max_trial_time_four_threads = trial_time;
        }
    }
    
    cout << "\nGenerated " << num_trials << " voronoi diagrams with " << num_sites << " sites.\n";
    cout << "One thread:\n";
    cout << "Total time = " << total_time_one_thread.count() << " seconds.\n";
    cout << "Average = " << total_time_one_thread.count() / num_trials << " seconds.\n";
    cout << "Max = " << max_trial_time_one_thread.count() << " seconds.\n\n";
    cout << "Two threads:\n";
    cout << "Total time = " << total_time_two_threads.count() << " seconds.\n";
    cout << "Average = " << total_time_two_threads.count() / num_trials << " seconds.\n";
    cout << "Max = " << max_trial_time_two_threads.count() << " seconds.\n\n";
    cout << "Four threads:\n";
    cout << "Total time = " << total_time_four_threads.count() << " seconds.\n";
    cout << "Average = " << total_time_four_threads.count() / num_trials << " seconds.\n";
    cout << "Max = " << max_trial_time_four_threads.count() << " seconds.\n\n";
        
    return 0;
}

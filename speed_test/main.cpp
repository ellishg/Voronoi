//
//  main.cpp
//  speed_test
//
//  Created by Ellis Sparky Hoag on 6/13/16.
//  Copyright © 2016 Ellis Sparky Hoag. All rights reserved.
//

/*
 *  Records as of 7/18/16
 *
 *  100 sites:
 *  0.0007065 average
 *  0.000964 max
 *
 *  1000 sites:
 *  0.0118539 average
 *  0.019233 max
 *
 *  4000 sites:
 *  0.07808 average
 *  0.114018 max
 *
 *  10000 sites:
 *  0.294484 average
 *  0.3578 max
 *
 *  16000 sites:
 *  0.54561 average
 *  0.637421 max
 */

#include <iostream>
#include "voronoi_sphere.h"

using namespace std;
using namespace Voronoi;

const int num_sites = 4000;
const int num_trials = 100;

int main(int argc, const char * argv[]) {
    
    srand((unsigned int)time(NULL));
    
    clock_t t;
    
    float total_time = 0;
    
    float max_trial_time = 0;
    
    for (int trial = 1; trial <= num_trials; trial++)
    {
        
        unsigned int seed = rand();
        
        srand(seed);
        
        cout << trial << ") Seed = " << seed << endl;
        
        vector<tuple<float, float, float>> verts;
        
        const int precision = numeric_limits<int>::max();
        for (int i = 0; i < num_sites; i++)
        {
            float x = (rand() % precision) / (float)precision - 0.5f;
            float y = (rand() % precision) / (float)precision - 0.5f;
            float z = (rand() % precision) / (float)precision - 0.5f;
            
            float r = sqrt(x*x + y*y + z*z);
            
            tuple<float, float, float> point = make_tuple(x / r, y / r, z / r);
            
            verts.push_back(point);
        }

        t = clock();
        
        generate_voronoi(&verts);
        //generate_voronoi_parallelized(&verts);
        
        float trial_time = (clock() - t) / (float)CLOCKS_PER_SEC;
        
        total_time += trial_time;
        
        if (trial_time > max_trial_time)
        {
            max_trial_time = trial_time;
        }
    }
    cout << "\nGenerated " << num_trials << " voronoi diagrams with " << num_sites << " sites.\n\n";
    
    cout << "Total time = " << total_time << " seconds.\n";
    cout << "Average = " << total_time / num_trials << " seconds.\n";
    cout << "Max = " << max_trial_time << " seconds.\n\n";
        
    return 0;
}

//
//  main.cpp
//  speed_test
//
//  Created by Ellis Sparky Hoag on 6/13/16.
//  Copyright © 2016 Ellis Sparky Hoag. All rights reserved.
//

/*
 *  Records as of 6/17/16
 *
 *  100 sites:
 *  0.0008778 average
 *  0.001331    max
 *
 *  1000 sites:
 *  0.0127625 average
 *  0.02085 max
 *
 *  4000 sites:
 *  0.082358 average
 *  0.10838 max
 *
 *  10000 sites:
 *  0.298549 average
 *  0.369286 max
 *
 *  16000 sites:
 *  0.617974 average
 *  0.983193 max
 */

#include <iostream>
#include "VoronoiSphere.h"

using namespace std;
using namespace Voronoi;

const int num_sites = 1000;
const int num_trials = 100;

int main(int argc, const char * argv[]) {
    
    srand((unsigned int)time(NULL));
    
    clock_t t;
    
    VoronoiSphere voronoi;
        
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
        
        voronoi.reset();
        
        voronoi.set_sites(&verts);
        
        voronoi.generate_voronoi();
        
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

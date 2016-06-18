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

const int num_sites = 16000;
const int num_trials = 10;

int main(int argc, const char * argv[]) {
    
    srand((unsigned int)time(NULL));
    
    float points[3 * num_sites];
    
    clock_t t;
    
    VoronoiSphere voronoi;
        
    float total_time = 0;
    
    float max_trial_time = 0;
    
    for (int trial = 1; trial <= num_trials; trial++)
    {
        
        unsigned int seed = rand();
        
        srand(seed);
        
        cout << trial << ") Seed = " << seed << endl;
        
        const int precision = numeric_limits<int>::max();
        for (int i = 0; i < num_sites; i++)
        {
            float x = (rand() % precision) / (float)precision - 0.5f;
            float y = (rand() % precision) / (float)precision - 0.5f;
            float z = (rand() % precision) / (float)precision - 0.5f;
            
            float r = sqrt(x*x + y*y + z*z);
            points[3 * i] = x / r;
            points[3 * i + 1] = y / r;
            points[3 * i + 2] = z / r;
        }

        t = clock();
        
        VoronoiDiagramSphere voronoi_diagram = voronoi.generate_voronoi(points, num_sites);
        
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

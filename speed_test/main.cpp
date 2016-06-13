//
//  main.cpp
//  speed_test
//
//  Created by Ellis Sparky Hoag on 6/13/16.
//  Copyright © 2016 Ellis Sparky Hoag. All rights reserved.
//

#include <iostream>
#include "VoronoiSphere.h"

using namespace std;
using namespace Voronoi;

const int num_sites = 1000;
const int num_trials = 50;

int main(int argc, const char * argv[]) {
    
    unsigned int seed = (unsigned int)time(NULL);
    
    srand(seed);
    
    cout << "Seed = " << seed << endl;
    
    float points[3 * num_sites];
    
    clock_t t;
    
    VoronoiSphere voronoi;
        
    float total_time = 0;
    
    float max_trial_time = 0;
    
    for (int i = 0; i < num_trials; i++)
    {
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
    cout << "\nGenerated " << num_trials << " voronoi diagrams from " << num_sites << " sites.\n\n";
    
    cout << "Total time = " << total_time << " seconds.\n";
    cout << "Average = " << total_time / num_trials << " seconds.\n";
    cout << "Max = " << max_trial_time << " seconds.\n\n";
    
    return 0;
}

//
//  main.cpp
//  Voronoi
//
//  Created by Ellis Sparky Hoag on 6/1/16.
//  Copyright © 2016 Ellis Sparky Hoag. All rights reserved.
//

#include <iostream>
#include "SDL2/SDL.h"
#include "SDL2/SDL_opengl.h"

#include "Voronoi2D.h"
#include "VoronoiSphere.h"

#define SPHERICAL_MODE

using namespace std;
using namespace Voronoi;

SDL_Window * window;
SDL_Event event;

const int num_sites = 100;

#ifdef SPHERICAL_MODE
float points[num_sites * 3];
#else
float points[num_sites * 2];
#endif

const float point_size = 5.f;
const float line_width = 1.f;

void render_points(float * points, int num_points, float r, float g, float b) {
    
    glPointSize(point_size);
    
    glColor3f(r, g, b);
    
    glBegin(GL_POINTS);
        for (int i = 0; i < num_points; i++) {
#ifdef SPHERICAL_MODE
            glVertex3f(points[3 * i], points[3 * i + 1], points[3 * i + 2]);
#else
            glVertex3f(points[2 * i], points[2 * i + 1], 0.f);
#endif
        }
    glEnd();
}

void render_edges_sphere(vector<HalfEdgeSphere *> edges) {
    
    glLineWidth(line_width);
    glPointSize(point_size);
    
    //const int num_segments = 100;
    
    for (int i = 0; i < edges.size(); i++) {
        
        if (!edges[i]->is_finished) {continue;}
        if (edges[i]->start == NULL || edges[i]->end == NULL) {continue;}
        
        glColor3f(1.f, 1.f, 0.f);
        
        //float theta0 = edges[i]->start.theta;
        //float theta1 = edges[i]->end.theta;
        //float phi0 = edges[i]->start.phi;
        //float phi1 = edges[i]->end.phi;
        
        //PointCartesian start = edges[i]->start;
        //PointCartesian end = edges[i]->end.get_cartesian();
        
        
        glBegin(GL_LINE_STRIP);
            glVertex3f(edges[i]->start->x, edges[i]->start->y, edges[i]->start->z);
            glVertex3f(edges[i]->end->x, edges[i]->end->y, edges[i]->end->z);
    
            //for (int i = 0; i < num_segments; i++) {
            //    glColor3f(1.f, (float)i / num_segments, 0.f);

            //    float theta = theta0 + (theta1 - theta0) * i / num_segments;
            //    float phi = phi0 + (phi1 - phi0) * i / num_segments;
            //    glVertex3f(sin(theta) * cos(phi), sin(theta) * sin(phi), cos(theta));
            //}
        glEnd();
        
        glColor3f(0, 0, 1);
        glBegin(GL_LINES);
            glVertex3f(edges[i]->cells[0]->site.x, edges[i]->cells[0]->site.y, edges[i]->cells[0]->site.z);
            glVertex3f(edges[i]->cells[1]->site.x, edges[i]->cells[1]->site.y, edges[i]->cells[1]->site.z);
        glEnd();
    }
}

void render_edges_2d(vector<HalfEdge2D *> edges) {
    
    glLineWidth(line_width);
    glPointSize(point_size);
    
    for (int i = 0; i < edges.size(); i++) {
        
        glColor3f(0.f, 0.f, 0.f);
        
        glBegin(GL_LINES);
            glVertex3f(edges[i]->start.x, edges[i]->start.y, 0.f);
            glVertex3f(edges[i]->end.x, edges[i]->end.y, 0.f);
        glEnd();
        
        glColor3f(0.f, 1.f, 0.f);
        
        //glBegin(GL_POINTS);
        //    glVertex3f(edges[i]->start.x, edges[i]->start.y, 0.f);
        //    glVertex3f(edges[i]->end.x, edges[i]->end.y, 0.f);
        //glEnd();
    }
}

void render_sphere(vector<HalfEdgeSphere *> edges, float sweep_line = 0) {
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    glLoadIdentity();
    
    glRotatef(-90, 1, 0, 0);
    
    float angle = SDL_GetTicks() / 50.f;
    
    glRotatef(angle, 0, 0, 1);
    //glRotatef(angle / 3.f, 1, 0, 0);
    
    render_edges_sphere(edges);
    render_points(points, num_sites, 1, 0, 0);
    
    if (sweep_line != 0)
    {
        glColor3f(0, 1, 1.f);
        glBegin(GL_LINE_LOOP);
            const int num_steps = 100;
            for (int i = 0; i < num_steps; i++) {
                glVertex3f(sin(sweep_line) * cos(2.f * M_PI * i / num_steps), sin(sweep_line) * sin(2.f * M_PI * i / num_steps), cos(sweep_line));
            }
        glEnd();
    }
    
    SDL_GL_SwapWindow(window);
    
}


void render_2d(vector<HalfEdge2D *> edges, float sweep_line = 0) {
    glClear(GL_COLOR_BUFFER_BIT);
    glLoadIdentity();

    render_points(points, num_sites, 1, 0, 0);
    render_edges_2d(edges);
    
    glBegin(GL_LINES);
        glVertex3f(sweep_line, 0, 0);
        glVertex3f(sweep_line, 1, 0);
    glEnd();
    
    SDL_GL_SwapWindow(window);

}

void sleep() {
    //SDL_Delay(25);
    bool delay = true;
    
    while (delay) {
        while (SDL_PollEvent(&event)) {
            switch (event.type) {
                case SDL_KEYDOWN:
                    switch (event.key.keysym.sym) {
                        case SDLK_SPACE:
                            delay = false;
                            break;
                        case SDLK_UP:
                            glRotatef(10, 1, 0, 0);
                            break;
                        case SDLK_DOWN:
                            glRotatef(10, -1, 0, 0);
                            break;
                        case SDLK_LEFT:
                            glRotatef(10, 0, 0, 1);
                            break;
                        case SDLK_RIGHT:
                            glRotatef(10, 0, 0, -1);
                            break;
                        default:
                            break;
                    }
                    break;
                default:
                    break;
            }
        }
    }
}

int main(int argc, const char * argv[]) {
    
    /*
     *  Bad seeds:
     *  1465266993 (100 sites) Good!
     *  1465267176 (100 sites) Good!
     *  1465267270 (100 sites) Good!
     *  1465267285 (1000 sites) Good!
     *  1465483687 (1000 sites)
     *  1465483717 (1000 sites)
     *  1465495088 (500 sites)
     *  1465505027 (100 sites)
     */
    
    unsigned int seed = (unsigned int)time(NULL);
    //unsigned int seed = 1465505027;
    
    cout << "Seed = " << seed << endl;
    
    srand(seed);
    
    if(SDL_Init(SDL_INIT_VIDEO))   {
        cout << "Error initializing video.\n";
        return 0;
    }
    
    window = SDL_CreateWindow("Voronoi", SDL_WINDOWPOS_UNDEFINED, SDL_WINDOWPOS_UNDEFINED, 500, 500, SDL_WINDOW_OPENGL);
    
    if(window == NULL)   {
        cout << "Error creating window.\n";
        return 0;
    }

    SDL_GLContext context = SDL_GL_CreateContext(window);
    
    if (context == NULL) {
        cout << "Error creating context.\n";
        return 0;
    }
    
    glMatrixMode(GL_PROJECTION);
    
    glLoadIdentity();
    
#ifdef SPHERICAL_MODE
    glOrtho(-1, 1, -1, 1, 0, 10);
    glClearColor(0.f, 0.f, 0.f, 1.f);
#else
    glOrtho(0, 1, 0, 1, 0, 10);
    glClearColor(1.f, 1.f, 1.f, 1.f);
#endif
    
    glMatrixMode(GL_MODELVIEW);
    
    bool is_running = true;
    
    const int precision = numeric_limits<int>::max();
    //const int precision = 10;
    for (int i = 0; i < num_sites; i++) {
#ifdef SPHERICAL_MODE
        float x = (rand() % precision) / (float)precision - 0.5f;
        float y = (rand() % precision) / (float)precision - 0.5f;
        float z = (rand() % precision) / (float)precision - 0.5f;
        float r = sqrt(x*x + y*y + z*z);
        points[3 * i] = x / r;
        points[3 * i + 1] = y / r;
        points[3 * i + 2] = z / r;
        //cout << "(" << points[3 * i] << ", " << points[3 * i + 1] << ", " << points[3 * i + 2] << ")\n";
#else
        points[2 * i] = (rand() % precision) / (float)precision;
        points[2 * i + 1] = (rand() % precision) / (float)precision;
#endif
    }
    
#ifdef SPHERICAL_MODE
    glRotatef(-90, 1, 0, 0);

    clock_t t = clock();
    
    VoronoiSphere voronoi;
    VoronoiDiagramSphere voronoi_diagram = voronoi.generate_voronoi(points, num_sites, render_sphere, sleep);
    
    cout << "Generated voronoi from " << num_sites << " sites in " << (clock() - t) / (float)CLOCKS_PER_SEC << " seconds.\n";
    
#else
    Voronoi2D voronoi;
    voronoi.generate_voronoi_2D(points, num_sites, render_2d, sleep);
    vector<HalfEdge2D *> edges = voronoi.get_voronoi_edges();
#endif
    
    while(is_running) {
        
#ifdef SPHERICAL_MODE
        render_sphere(voronoi_diagram.edges);
#else
        render_2d(edges);
#endif
        
        while (SDL_PollEvent(&event)) {
            switch (event.type) {
                case SDL_QUIT:
                    is_running = false;
                    break;
                case SDL_KEYDOWN:
                    switch (event.key.keysym.sym) {
                        case SDLK_ESCAPE:
                            is_running = false;
                            break;
                        case SDLK_UP:
                            glRotatef(10, 1, 0, 0);
                            break;
                        case SDLK_DOWN:
                            glRotatef(10, -1, 0, 0);
                            break;
                        case SDLK_LEFT:
                            glRotatef(10, 0, 0, 1);
                            break;
                        case SDLK_RIGHT:
                            glRotatef(10, 0, 0, -1);
                            break;
                        default:
                            break;
                    }
                    break;
                default:
                    break;
            }
        }

    }
        
    SDL_DestroyWindow(window);
    
    SDL_GL_DeleteContext(context);
    
    SDL_Quit();
    
    return 0;
}

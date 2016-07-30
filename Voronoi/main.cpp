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
#include "voronoi_sphere.h"

#define SPHERICAL_MODE

using namespace std;
using namespace Voronoi;

SDL_Window * window;
SDL_Event event;
SDL_GLContext context;

const int window_size = 600;

int num_sites = 4000;

//const THREAD_NUMBER num_threads = ONE_THREAD;
const THREAD_NUMBER num_threads = TWO_THREADS;
//const THREAD_NUMBER num_threads = FOUR_THREADS;

/*
 *  Bad seeds:
 *  1741754784 (8000 sites)
 *  1252816138 (16000 sites)
 *  319410050
 *  660258475
 */

unsigned int seed = (unsigned int)time(NULL);
//unsigned int seed = 1741754784;

VoronoiDiagramSphere voronoi_diagram;

#ifndef SPHERICAL_MODE
float points[num_sites * 2];
#endif

void render_voronoi_sphere(VoronoiDiagramSphere voronoi_diagram, float sweep_line = 0)
{
    bool should_rotate = true;
    bool render_sites = true;
    bool render_voronoi = true;
    bool render_delaunay = true;
    bool render_sweep_line = true && sweep_line != 0;
 
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    
    glPointSize(2.f);
    glLineWidth(1.f);
    
    if (should_rotate)
    {
        glLoadIdentity();
        
        glRotatef(-90, 1, 0, 0);
        
        float angle = SDL_GetTicks() / 50.f;
        
        glRotatef(angle, 0, 0, 1);
        glRotatef(angle, 0, 1, 0);
    }
    
    if (render_voronoi)
    {
        glColor3f(1, 1, 0);
        glBegin(GL_LINES);
            for (auto edge : voronoi_diagram.voronoi_edges)
            {
                PointCartesian start = voronoi_diagram.voronoi_vertices[edge.vidx[0]];
                PointCartesian end = voronoi_diagram.voronoi_vertices[edge.vidx[1]];
                
                glVertex3f((float)start.x, (float)start.y, (float)start.z);
                glVertex3f((float)end.x, (float)end.y, (float)end.z);
            }
        glEnd();
    }
    
    if (render_delaunay)
    {
        glColor3f(0, 0, 1);
        glBegin(GL_LINES);
            for (auto edge : voronoi_diagram.delaunay_edges)
            {
                PointCartesian start = voronoi_diagram.sites[edge.vidx[0]];
                PointCartesian end = voronoi_diagram.sites[edge.vidx[1]];
                glVertex3f((float)start.x, (float)start.y, (float)start.z);
                glVertex3f((float)end.x, (float)end.y, (float)end.z);
            }
        glEnd();
    }
    
    if (render_sites)
    {
        glColor3f(1, 0, 0);
        glBegin(GL_POINTS);
            for (auto site : voronoi_diagram.sites)
            {
                glVertex3f((float)site.x, (float)site.y, (float)site.z);
            }
        glEnd();
    }
    
    if (render_sweep_line)
    {
        glColor3f(0, 1, 1);
        glBegin(GL_LINE_LOOP);
            const int num_steps = 20;
            for (int i = 0; i < num_steps; i++)
            {
                glVertex3f(sin(sweep_line) * cos(2.f * M_PI * i / num_steps), sin(sweep_line) * sin(2.f * M_PI * i / num_steps), cos(sweep_line));
            }
        glEnd();
    }
    
    SDL_GL_SwapWindow(window);
}

void quit()
{
    SDL_DestroyWindow(window);
    
    SDL_GL_DeleteContext(context);
    
    SDL_Quit();
    
    exit(0);
}

void render_points(float * points, int num_points, float r, float g, float b)
{
    
    glPointSize(2.f);
    
    glColor3f(r, g, b);
    
    glBegin(GL_POINTS);
        for (int i = 0; i < num_points; i++)
        {
#ifdef SPHERICAL_MODE
            glVertex3f(points[3 * i], points[3 * i + 1], points[3 * i + 2]);
#else
            glVertex3f(points[2 * i], points[2 * i + 1], 0.f);
#endif
        }
    glEnd();
}

void render_edges_2d(vector<HalfEdge2D *> edges)
{
    
    glLineWidth(1.f);
    glPointSize(2.f);
    
    for (int i = 0; i < edges.size(); i++)
    {
        
        glColor3f(0.f, 0.f, 0.f);
        
        glBegin(GL_LINES);
            glVertex3f(edges[i]->start.x, edges[i]->start.y, 0.f);
            glVertex3f(edges[i]->end.x, edges[i]->end.y, 0.f);
        glEnd();
        
        glColor3f(0.f, 1.f, 0.f);
    }
}

#ifndef SPHERICAL_MODE
void render_2d(vector<HalfEdge2D *> edges, float sweep_line = 0)
{
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
#endif

bool is_sleeping()
{
    
    while (SDL_PollEvent(&event))
    {
        switch (event.type)
        {
            case SDL_QUIT:
                quit();
                break;
            case SDL_KEYDOWN:
                switch (event.key.keysym.sym)
                {
                    case SDLK_ESCAPE:
                        quit();
                        break;
                    case SDLK_SPACE:
                        return false;
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
    
    return true;
}

int main(int argc, const char * argv[])
{
    if (argc == 2)
    {
        num_sites = atoi(argv[1]);
    }
    
    cout << "Seed = " << seed << endl;
    
    srand(seed);
    
    if(SDL_Init(SDL_INIT_VIDEO))
    {
        cout << "Error initializing video.\n";
        return 0;
    }
    
    window = SDL_CreateWindow("Voronoi", SDL_WINDOWPOS_UNDEFINED, SDL_WINDOWPOS_UNDEFINED, window_size, window_size, SDL_WINDOW_OPENGL);
    
    if(window == NULL)
    {
        cout << "Error creating window.\n";
        return 0;
    }

    context = SDL_GL_CreateContext(window);
    
    if (context == NULL)
    {
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
    
    vector<tuple<float, float, float>> verts;
    
    const int precision = numeric_limits<int>::max();
    //const int precision = 10;
    for (int i = 0; i < num_sites; i++)
    {
#ifdef SPHERICAL_MODE
        float x = (rand() % precision) / (float)precision - 0.5f;
        float y = (rand() % precision) / (float)precision - 0.5f;
        float z = (rand() % precision) / (float)precision - 0.5f;

        float r = sqrt(x*x + y*y + z*z);
        
        tuple<float, float, float> point = make_tuple(x / r, y / r, z / r);
        
        verts.push_back(point);
#else
        points[2 * i] = (rand() % precision) / (float)precision;
        points[2 * i + 1] = (rand() % precision) / (float)precision;
#endif
    }
    
#ifdef SPHERICAL_MODE
    glRotatef(-90, 1, 0, 0);

    auto start_time = chrono::system_clock::now();
    
    voronoi_diagram = generate_voronoi(&verts, num_threads, render_voronoi_sphere, is_sleeping);
    
    chrono::duration<float> run_time = chrono::system_clock::now() - start_time;
    
    cout << "Generated voronoi from " << num_sites << " sites in " << run_time.count() << " seconds.\n";
    
#else
    Voronoi2D voronoi;
    voronoi.generate_voronoi_2D(points, num_sites, render_2d, sleep);
    vector<HalfEdge2D *> edges = voronoi.get_voronoi_edges();
#endif
    
    while(true)
    {
        
#ifdef SPHERICAL_MODE
        render_voronoi_sphere(voronoi_diagram);
#else
        render_2d(edges);
#endif
        
        while (SDL_PollEvent(&event))
        {
            switch (event.type)
            {
                case SDL_QUIT:
                    quit();
                    break;
                case SDL_KEYDOWN:
                    switch (event.key.keysym.sym)
                    {
                        case SDLK_ESCAPE:
                            quit();
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
                        case SDLK_r:
                            seed = rand();
                            cout << "Seed = " << seed << endl;
                            srand(seed);
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

                            voronoi_diagram = generate_voronoi(&verts, num_threads);
                            break;
                    }
                    break;
                default:
                    break;
            }
        }

    }
    
    return 0;
}

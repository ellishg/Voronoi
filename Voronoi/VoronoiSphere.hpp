//
//  VoronoiSphere.hpp
//  Voronoi
//
//  Created by Ellis Sparky Hoag on 6/3/16.
//  Copyright © 2016 Ellis Sparky Hoag. All rights reserved.
//

#ifndef VoronoiSphere_hpp
#define VoronoiSphere_hpp

#include <iostream>
#include <vector>
#include <queue>
#include "math.h"

using namespace std;

namespace Voronoi {
    
    class VoronoiSphere;
    struct PointSphere;
    struct PointCartesian;
    struct HalfEdgeSphere;
    struct VoronoiCellSphere;
    struct VoronoiDiagramSphere;
    
    struct PointCartesian {
        
        PointCartesian(float a = 0, float b = 0, float c = 0) : x(a), y(b), z(c) {}
                
        friend PointCartesian operator-(const PointCartesian left, const PointCartesian right) {return PointCartesian(left.x - right.x, left.y - right.y, left.z - right.z);}
        
        friend ostream & operator<<(ostream & out, const PointCartesian & p) {return out << "(" << p.x << ", " << p.y << ", " << p.z << ")\n";}
        
        void normalize() {
            float r = sqrt(x*x + y*y + z*z);
            x /= r;
            y /= r;
            z /= r;
        }
        
        void to_unit_sphere(float & theta, float & phi) {
            theta = acos(z);
            phi = atan2(y, x);
        }
        
        static PointCartesian cross_product(PointCartesian left, PointCartesian right) {
            float x = left.y * right.z - left.z * right.y;
            float y = left.z * right.x - left.x * right.z;
            float z = left.x * right.y - left.y * right.x;
            return PointCartesian(x, y, z);
        }
        
        float x, y, z;
    };
    
    /*
     *  Points on the sphere are stored as (theta, phi).
     *  0 <= theta < 2PI is the angle from the north pole.
     *  -PI < phi <= PI is the angle around the sphere from the x-axis.
     */
    struct PointSphere {
        
        float theta, phi;
        
        float x, y, z;
        
        PointSphere(float t = 0, float p = 0) : theta(t), phi(p), x(0), y(0), z(0), has_cartesian(false) {}
        
        PointSphere(float _x, float _y, float _z) : theta(acos(_z)), phi(atan2(_y, _x)), x(_x), y(_y), z(_z), has_cartesian(true) {}
        
        PointSphere(PointCartesian point) : PointSphere(point.x, point.y, point.z) {}
        
        friend bool operator<(const PointSphere & left, const PointSphere & right) {return left.theta == right.theta ? left.phi < right.phi : left.theta < right.theta;}
        friend bool operator>(const PointSphere & left, const PointSphere & right) {return left.theta == right.theta ? left.phi > right.phi : left.theta > right.theta;}
        
        friend ostream & operator<<(ostream & out, const PointSphere & p) {return out << "(" << p.theta << ", " << p.phi << ")\n";}
        
        PointCartesian get_cartesian() {
            set_cartesian();
            return PointCartesian(x, y, z);
        }
        
        static float angle_between(PointSphere a, PointSphere b) {
            a.set_cartesian();
            b.set_cartesian();
            return acos(a.x * b.x + a.y * b.y + a.z * b.z);
        }
        
    private:
        
        void set_cartesian() {
            if (!has_cartesian) {
                float sin_theta = sin(theta);
                x = cos(phi) * sin_theta;
                y = sin(phi) * sin_theta;
                z = cos(theta);
                has_cartesian = true;
            }
        }
        
        bool has_cartesian;
    };
    
    struct VoronoiCellSphere {
        
        VoronoiCellSphere(PointSphere point) : site(point) {}
        
        VoronoiCellSphere(PointCartesian point) : site(point) {}
        
        PointSphere site;
        
        vector<HalfEdgeSphere *> edges;
    };
    
    struct HalfEdgeSphere {
        
        HalfEdgeSphere(PointCartesian * s, VoronoiCellSphere * a, VoronoiCellSphere * b) : start(s), end(NULL), is_finished(false) {
            cells[0] = a;
            cells[1] = b;
        }
        
        void finish(PointCartesian * e) {
            if (!is_finished) {
                end = e;
                is_finished = true;
            }
        }
        
        PointCartesian *start, *end;
        
        VoronoiCellSphere * cells[2];
        
        bool is_finished;
    };
    
    struct VoronoiDiagramSphere {
        
        vector<HalfEdgeSphere *> edges;
        
        vector<PointCartesian *> verticies;
        
        vector<VoronoiCellSphere *> cells;
    };
    
    
    class VoronoiSphere {
        
    public:
        
        VoronoiDiagramSphere generate_voronoi(float * points, int num_points, void (*render)(vector<HalfEdgeSphere *>, float), void (*sleep)());
        
    private:
        
        struct PriorityQueueCompare;
        struct CircleEventSphere;
        struct ArcSphere;
        
        /*
         *  Sort by theta. Smallest theta should be first
         */
        struct PriorityQueueCompare {
            bool operator()(VoronoiCellSphere * left, VoronoiCellSphere * right) {return left->site > right->site;}
            //bool operator()(PointSphere & left, PointSphere & right) {return left > right;}
            bool operator()(CircleEventSphere * left, CircleEventSphere * right) {return left->lowest_theta > right->lowest_theta;}
        };
        
        struct CircleEventSphere {
            
            CircleEventSphere(ArcSphere * a, PointSphere c, float l) : arc(a), circumcenter(c), lowest_theta(l), is_valid(true) {}
            
            ArcSphere * arc;
            
            PointSphere circumcenter;
            
            float lowest_theta;
            
            bool is_valid;
        };
        
        struct ArcSphere {
            
            static const int max_skiplist_height = 50;


            ArcSphere(VoronoiCellSphere * c) : cell(c), event(NULL), s0(NULL), s1(NULL) {
                
                height = (rand() % max_skiplist_height) + 1;
                
                prev = new ArcSphere*[height];
                next = new ArcSphere*[height];
                
                for (int i = 0; i < height; i++) {
                    prev[i] = next[i] = this;
                }
            }
            
            ArcSphere(VoronoiCellSphere * c, ArcSphere * left, ArcSphere * right) : cell(c), event(NULL), s0(NULL), s1(NULL) {
                
                height = (rand() % max_skiplist_height) + 1;
                
                prev = new ArcSphere*[height];
                next = new ArcSphere*[height];
                
                for (int i = 0; i < height; i++) {
                    
                    while (left->height <= i) {
                        left = left->prev[left->height - 1];
                    }
                    while (right->height <= i) {
                        right = right->next[right->height - 1];
                    }
                    
                    left->next[i] = this;
                    prev[i] = left;
                    
                    right->prev[i] = this;
                    next[i] = right;
                    
                }
                
            }
            
            ~ArcSphere() {
                
                for (int i = 0; i < height; i++) {
                    ArcSphere * left = prev[i];
                    left->next[i] = next[i];
                    next[i]->prev[i] = left;
                }
                
                delete [] prev;
                delete [] next;
            }
            
            void add_edge_left(HalfEdgeSphere * edge) {
                s0 = edge;
                cell->edges.push_back(edge);
                //maybe add cell to edge
            }
            
            void add_edge_right(HalfEdgeSphere * edge) {
                s1 = edge;
                cell->edges.push_back(edge);
                //maybe add cell to edge
            }
            VoronoiCellSphere * cell;
            
            ArcSphere **prev, **next;
            
            CircleEventSphere * event;
            
            HalfEdgeSphere *s0, *s1;
            
            int height;
        };

        void handle_site_event(VoronoiCellSphere * cell);
        
        void handle_circle_event(CircleEventSphere * event);
                
        bool parabolic_intersection(PointSphere left, PointSphere right, PointSphere & intersection);
        
        PointSphere phi_to_point(PointSphere arc, float phi);
        
        void check_circle_event(ArcSphere * arc);
        
        void make_circle(PointSphere a, PointSphere b, PointSphere c, PointSphere & circumcenter, float & lowest_theta);
        
        void finalize_diagram();
        
        priority_queue<CircleEventSphere, vector<CircleEventSphere *>, PriorityQueueCompare> circle_event_queue;
        
        ArcSphere * beach_head;
                
        float sweep_line;
        
        VoronoiDiagramSphere voronoi_diagram;
    };
}

#endif /* VoronoiSphere_hpp */

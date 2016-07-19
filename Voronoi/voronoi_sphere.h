//
//  VoronoiSphere.h
//  Voronoi
//
//  Created by Ellis Sparky Hoag on 6/3/16.
//  Copyright © 2016 Ellis Sparky Hoag. All rights reserved.
//

#ifndef VoronoiSphere_h
#define VoronoiSphere_h

#include <iostream>
#include <vector>
#include <queue>
#include <cmath>
#include <thread>
#include <assert.h>

namespace Voronoi
{
    typedef double Real;
    
    struct Edge;
    struct VoronoiDiagramSphere;
    struct PointCartesian;
    struct CircleEventSphere;
    struct ArcSphere;
    struct PointSphere;
    struct VoronoiCellSphere;
    struct HalfEdgeSphere;
    struct CompareTopDown;
    struct CompareBottomUp;
    
    VoronoiDiagramSphere generate_voronoi(std::vector<std::tuple<float, float, float>> * verts, void (*render)(VoronoiDiagramSphere, float) = NULL, bool (*is_sleeping)() = NULL);
    
    VoronoiDiagramSphere generate_voronoi_parallelized(std::vector<std::tuple<float, float, float>> * verts);
    
    struct Edge
    {
        Edge(int start_idx = 0, int end_idx = 0)
        {
            vidx[0] = start_idx;
            vidx[1] = end_idx;
        }
        int vidx[2];
    };
    
    struct VoronoiDiagramSphere
    {
        std::vector<PointCartesian> sites, voronoi_verticies;
        
        std::vector<Edge> voronoi_edges, delaunay_edges;
    };
    
    struct PointCartesian
    {
        PointCartesian(Real a = 0, Real b = 0, Real c = 0) : x(a), y(b), z(c) {}
                
        inline friend PointCartesian operator-(const PointCartesian left, const PointCartesian right) {return PointCartesian(left.x - right.x, left.y - right.y, left.z - right.z);}
        
        friend std::ostream & operator<<(std::ostream & out, const PointCartesian & p) {return out << "(" << p.x << ", " << p.y << ", " << p.z << ")\n";}
        
        inline void normalize()
        {
            Real r = sqrt(x*x + y*y + z*z);
            if (r == 0)
            {
                std::cout << "r = 0\n";
                x = y = z = 0;
            }
            else
            {
                x /= r;
                y /= r;
                z /= r;
            }
        }
        
        inline static PointCartesian cross_product(const PointCartesian & left, const PointCartesian & right)
        {
            Real x = left.y * right.z - left.z * right.y;
            Real y = left.z * right.x - left.x * right.z;
            Real z = left.x * right.y - left.y * right.x;
            return PointCartesian(x, y, z);
        }
        
        Real x, y, z;
    };
    
    /*
     *  Points on the sphere are stored as (theta, phi).
     *  0 <= theta < 2PI is the angle from the north pole.
     *  -PI < phi <= PI is the angle around the sphere from the x-axis.
     */
    struct PointSphere
    {
        Real theta, phi;
        
        Real x, y, z;
        
        PointSphere(Real t = 0, Real p = 0) : theta(t), phi(p), x(0), y(0), z(0), has_cartesian(false) {}
        
        PointSphere(Real _x, Real _y, Real _z) : theta(acos(_z)), phi(atan2(_y, _x)), x(_x), y(_y), z(_z), has_cartesian(true) {}
        
        PointSphere(PointCartesian point) : PointSphere(point.x, point.y, point.z) {}
        
        inline friend bool operator<(const PointSphere & left, const PointSphere & right) {return left.theta == right.theta ? left.phi < right.phi : left.theta < right.theta;}
        inline friend bool operator>(const PointSphere & left, const PointSphere & right) {return left.theta == right.theta ? left.phi > right.phi : left.theta > right.theta;}
        friend bool operator==(const PointSphere & left, const PointSphere & right) {return left.theta == right.theta && left.phi == right.phi;}
        
        friend std::ostream & operator<<(std::ostream & out, const PointSphere & p) {return out << "(" << p.theta << ", " << p.phi << ")\n";}
        
        inline PointCartesian get_cartesian()
        {
            set_cartesian();
            return PointCartesian(x, y, z);
        }
        
        inline static Real angle_between(PointSphere & a, PointSphere & b)
        {
            a.set_cartesian();
            b.set_cartesian();
            return acos(a.x * b.x + a.y * b.y + a.z * b.z);
        }
        
    private:
        
        inline void set_cartesian()
        {
            if (!has_cartesian)
            {
                Real sin_theta = sin(theta);
                x = cos(phi) * sin_theta;
                y = sin(phi) * sin_theta;
                z = cos(theta);
                has_cartesian = true;
            }
        }
        
        bool has_cartesian;
    };
    
    struct ArcSphere
    {
        ArcSphere(int _cell_idx, int _height) : cell_idx(_cell_idx), height(_height), event(NULL), left_edge_idx(-1), right_edge_idx(-1) {}
        
        int cell_idx;
        
        int height;
        
        ArcSphere **prev, **next;
        
        CircleEventSphere * event;
        
        int left_edge_idx, right_edge_idx;
    };
    
    struct CircleEventSphere
    {
        CircleEventSphere(ArcSphere * a, PointSphere c, Real l) : arc(a), circumcenter(c), lowest_theta(l), is_valid(true) {}
        
        PointSphere circumcenter;
        
        bool is_valid;
        
        ArcSphere * arc;
        
        Real lowest_theta;
    };

    struct VoronoiCellSphere
    {
        VoronoiCellSphere(PointSphere point, int _idx) : site(point), cell_idx(_idx) {}
        
        VoronoiCellSphere(PointCartesian point, int _idx) : site(point), cell_idx(_idx) {}
        
        PointSphere site;
        
        int cell_idx;
        
        std::vector<int> edge_ids;
    };
    
    struct HalfEdgeSphere
    {
        HalfEdgeSphere(int vidx) : start_idx(vidx), is_finished(false) {}
        
        int start_idx, end_idx;
        
        bool is_finished;
    };
    
    /*
     *  Sort by theta. Smallest theta should be first.
     */
    struct PriorityQueueCompare
    {
        bool operator()(VoronoiCellSphere & left, VoronoiCellSphere & right) {return left.site > right.site;}
        bool operator()(CircleEventSphere * left, CircleEventSphere * right) {return left->lowest_theta > right->lowest_theta;}
    };
    
    
    void compute_priority_queues(VoronoiDiagramSphere * voronoi_diagram, std::vector<std::tuple<float, float, float>> * verts, bool is_top_down);
    
    void handle_site_event(VoronoiCellSphere cell, VoronoiDiagramSphere * voronoi_diagram, std::vector<VoronoiCellSphere> * cells, std::vector<HalfEdgeSphere> * half_edges, std::priority_queue<CircleEventSphere *, std::vector<CircleEventSphere *>, PriorityQueueCompare> * circle_event_queue_ptr, ArcSphere * & beach_head, Real sweep_line, Real sin_sweep_line, Real cos_sweep_line);
    
    void handle_circle_event(CircleEventSphere * event, VoronoiDiagramSphere * voronoi_diagram, std::vector<VoronoiCellSphere> * cells, std::vector<HalfEdgeSphere> * half_edges, std::priority_queue<CircleEventSphere *, std::vector<CircleEventSphere *>, PriorityQueueCompare> * circle_event_queue_ptr, ArcSphere * & beach_head);
    
    bool parabolic_intersection(PointSphere left, PointSphere right, Real & phi_intersection, ArcSphere * & beach_head, Real sweep_line, Real sin_sweep_line, Real cos_sweep_line);
    
    inline PointSphere phi_to_point(PointSphere arc, Real phi, Real sweep_line, Real sin_sweep_line, Real cos_sweep_line);
    
    void check_circle_event(ArcSphere * arc, std::priority_queue<CircleEventSphere *, std::vector<CircleEventSphere *>, PriorityQueueCompare> * circle_event_queue_ptr, std::vector<VoronoiCellSphere> * cells);
    
    inline void make_circle(PointSphere a, PointSphere b, PointSphere c, PointSphere & circumcenter, Real & lowest_theta);
    
    void add_half_edge_sphere(VoronoiDiagramSphere * voronoi_diagram, std::vector<VoronoiCellSphere> * cells, std::vector<HalfEdgeSphere> * half_edges, PointCartesian start, ArcSphere * left, ArcSphere * right);
    
    void finish_half_edge_sphere(VoronoiDiagramSphere * voronoi_diagram, std::vector<HalfEdgeSphere> * half_edges, int edge_idx, PointCartesian end);
    
    void add_initial_arc_sphere(int cell_id, ArcSphere * & beach_head);
    
    void add_arc_sphere(int cell_id, ArcSphere * left, ArcSphere * right, ArcSphere * & beach_head);
    
    int random_height();
    
    void remove_arc_sphere(ArcSphere * arc, ArcSphere * & beach_head);
    
    ArcSphere * traverse_skiplist_to_site(ArcSphere * arc, Real phi, std::vector<VoronoiCellSphere> * cells);
}

#endif /* VoronoiSphere_h */
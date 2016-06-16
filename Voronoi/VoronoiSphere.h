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
#include "math.h"

#include "assert.h"

namespace Voronoi
{
    typedef double Real;
    
    class VoronoiSphere;
    struct PointSphere;
    struct PointCartesian;
    struct HalfEdgeSphere;
    struct VoronoiCellSphere;
    struct VoronoiDiagramSphere;
    struct ArcSphere;
    struct CircleEventSphere;
    
    struct PointCartesian
    {
        PointCartesian(Real a = 0, Real b = 0, Real c = 0) : x(a), y(b), z(c) {}
                
        friend PointCartesian operator-(const PointCartesian left, const PointCartesian right) {return PointCartesian(left.x - right.x, left.y - right.y, left.z - right.z);}
        
        friend std::ostream & operator<<(std::ostream & out, const PointCartesian & p) {return out << "(" << p.x << ", " << p.y << ", " << p.z << ")\n";}
        
        void normalize()
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
        
        static PointCartesian cross_product(const PointCartesian & left, const PointCartesian & right)
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
        
        friend bool operator<(const PointSphere & left, const PointSphere & right) {return left.theta == right.theta ? left.phi < right.phi : left.theta < right.theta;}
        friend bool operator>(const PointSphere & left, const PointSphere & right) {return left.theta == right.theta ? left.phi > right.phi : left.theta > right.theta;}
        friend bool operator==(const PointSphere & left, const PointSphere & right) {return left.theta == right.theta && left.phi == right.phi;}
        
        friend std::ostream & operator<<(std::ostream & out, const PointSphere & p) {return out << "(" << p.theta << ", " << p.phi << ")\n";}
        
        PointCartesian get_cartesian()
        {
            set_cartesian();
            return PointCartesian(x, y, z);
        }
        
        static Real angle_between(PointSphere & a, PointSphere & b)
        {
            a.set_cartesian();
            b.set_cartesian();
            return acos(a.x * b.x + a.y * b.y + a.z * b.z);
        }
        
    private:
        
        void set_cartesian()
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
    
    struct VoronoiCellSphere
    {
        VoronoiCellSphere(PointSphere point, int _id) : site(point), cell_id(_id) {}
        
        VoronoiCellSphere(PointCartesian point, int _id) : site(point), cell_id(_id) {}
        
        PointSphere site;
        
        std::vector<int> edge_ids;
        
        int cell_id;
    };
    
    struct ArcSphere
    {        
        ArcSphere(int _cell_id, int _height) : cell_id(_cell_id), height(_height), event(NULL), left_edge_id(-1), right_edge_id(-1) {}
        
        int cell_id;
        
        int height;
        
        ArcSphere **prev, **next;
        
        CircleEventSphere * event;
        
        int left_edge_id, right_edge_id;
    };
    
    struct CircleEventSphere
    {
        CircleEventSphere(ArcSphere * a, PointSphere c, Real l) : arc(a), circumcenter(c), lowest_theta(l), is_valid(true) {}
        
        ArcSphere * arc;
        
        PointSphere circumcenter;
        
        Real lowest_theta;
        
        bool is_valid;
    };
    
    struct HalfEdgeSphere
    {
        HalfEdgeSphere(PointCartesian s, int left_cell_id, int right_cell_id) : start(s), is_finished(false)
        {
            cell_ids[0] = left_cell_id;
            cell_ids[1] = right_cell_id;
        }
        
        void finish(PointCartesian e)
        {
            if (!is_finished)
            {
                end = e;
                is_finished = true;
            }
        }
        
        PointCartesian start, end;
        
        int cell_ids[2];
        
        bool is_finished;
    };
    
    struct VoronoiDiagramSphere
    {
        std::vector<HalfEdgeSphere> edges;
        
        std::vector<VoronoiCellSphere> cells;
    };
    
    
    class VoronoiSphere
    {
    public:
        
        VoronoiDiagramSphere generate_voronoi(float * points, int num_points, void (*render)(std::vector<HalfEdgeSphere>, float, float *, int) = NULL, void (*sleep)() = NULL);
        
    private:
        
        const int max_skiplist_height = 15;
        const int skip_list_height_probability = 2;
        
        struct PriorityQueueCompare;
        
        /*
         *  Sort by theta. Smallest theta should be first.
         */
        struct PriorityQueueCompare
        {
            bool operator()(VoronoiCellSphere & left, VoronoiCellSphere & right) {return left.site > right.site;}
            bool operator()(CircleEventSphere * left, CircleEventSphere * right) {return left->lowest_theta > right->lowest_theta;}
        };
        
        ArcSphere * beach_head;
        
        Real sweep_line;
        
        VoronoiDiagramSphere voronoi_diagram;
                
        std::priority_queue<VoronoiCellSphere, std::vector<VoronoiCellSphere>, PriorityQueueCompare> site_event_queue;

        std::priority_queue<CircleEventSphere *, std::vector<CircleEventSphere *>, PriorityQueueCompare> circle_event_queue;

        
        
        void handle_site_event(VoronoiCellSphere cell);
        
        void handle_circle_event(CircleEventSphere * event);
                
        bool parabolic_intersection(PointSphere left, PointSphere right, PointSphere & intersection);
        
        PointSphere phi_to_point(PointSphere arc, Real phi);
        
        void check_circle_event(ArcSphere * arc);
        
        void make_circle(PointSphere a, PointSphere b, PointSphere c, PointSphere & circumcenter, Real & lowest_theta);
        
        void initialize_diagram(float * points, int num_points);
        
        void finalize_diagram();
        
        void add_half_edge_sphere(PointCartesian p, ArcSphere * left, ArcSphere * right);
        
        void add_initial_arc_sphere(int cell_id);
        
        void add_arc_sphere(int cell_id, ArcSphere * left, ArcSphere * right);
        
        int random_height();
        
        void remove_arc_sphere(ArcSphere * arc);
        
        ArcSphere * traverse_skiplist_to_site(ArcSphere * arc, Real phi);

    };
}

#endif /* VoronoiSphere_h */

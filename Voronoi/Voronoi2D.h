//
//  Voronoi2D.h
//  Voronoi
//
//  Created by Ellis Sparky Hoag on 6/1/16.
//  Copyright © 2016 Ellis Sparky Hoag. All rights reserved.
//

#ifndef Voronoi2D_h
#define Voronoi2D_h

#include <iostream>
#include <vector>
#include <queue>
#include "math.h"
#include <assert.h>
#include <limits>

//#define INF numeric_limits<float>::max()
#define INF 100000000.f

//using namespace std;

namespace Voronoi {

    class Voronoi2D;
    struct Point2D;
    struct CircleEvent2D;
    struct Arc2D;
    struct HalfEdge2D;
    
    struct Point2D {
        
        Point2D(float _x = 0, float _y = 0) : x(_x), y(_y) {}
        
        const friend bool operator<(const Point2D & left, const Point2D & right) {return left.x == right.x ? left.y < right.y : left.x < right.x;}
        const friend bool operator>(const Point2D & left, const Point2D & right) {return left.x == right.x ? left.y > right.y : left.x > right.x;}
        
        friend std::ostream & operator<<(std::ostream & out, const Point2D & p) {return out << "(" << p.x << ", " << p.y << ")\n";}
        
        float x, y;
    };
    
    struct CircleEvent2D {
        
        CircleEvent2D(Arc2D * a, Point2D c, float x) : arc(a), circumcenter(c), right_most_x(x), is_valid(true) {}
        
        float right_most_x;
        
        Point2D circumcenter;
        
        Arc2D * arc;
        
        bool is_valid;
        
    };
    
    class Voronoi2D {
        
    public:
        
        void generate_voronoi_2D(float * points, int num_points, void (*render)(std::vector<HalfEdge2D *>, float), void (*sleep)());
        
        std::vector<HalfEdge2D *> get_voronoi_edges() {return voronoi_edges;}
        
    private:
        
        struct PriorityQueueCompare {
            bool operator()(const Point2D & left, const Point2D & right) {return left > right;}
            bool operator()(const CircleEvent2D * left, const CircleEvent2D * right) {return left->right_most_x > right->right_most_x;}
        };

        
        void handle_site_event(Point2D site);
        
        void handle_circle_event(CircleEvent2D * event);
        
        bool does_intersect(Point2D point, Arc2D * arc, Point2D * ret);
        
        Point2D parabolic_intersection(Point2D left, Point2D right);
        
        void check_circle_event(Arc2D * arc);
        
        bool make_circle(Point2D a, Point2D b, Point2D c, Point2D & circumcenter, float & right_most_x);
        
        //vector<Point2D> voronoi_verticies;
        
        std::vector<HalfEdge2D *> voronoi_edges;
        
        std::priority_queue<CircleEvent2D, std::vector<CircleEvent2D *>, PriorityQueueCompare> circle_event_queue;
        
        Arc2D * beach_head;
        
        float sweep_line;
        
        float X0, X1, Y0, Y1;
        
    };

    struct HalfEdge2D {
        
        HalfEdge2D(Point2D s) : start(s), end(), is_done(false) {}
        
        void finish(Point2D e) {
            if (!is_done) {
                end = e;
                is_done = true;
            }
        }
        
        Point2D start, end;
        
        bool is_done;
    };
    
    struct Arc2D {
        
        Arc2D(Point2D l, Arc2D * p = NULL, Arc2D * n = NULL) : location(l), prev(p), next(n), event(NULL), s0(NULL), s1(NULL) {}
        
        Point2D location;
        Arc2D * prev, * next;
        CircleEvent2D * event;
        HalfEdge2D * s0, * s1;
    };
}

#endif /* Voronoi2D_h */

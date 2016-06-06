//
//  Voronoi2D.cpp
//  Voronoi
//
//  Created by Ellis Sparky Hoag on 6/1/16.
//  Copyright © 2016 Ellis Sparky Hoag. All rights reserved.
//

#include "Voronoi2D.hpp"

using namespace std;
using namespace Voronoi;

void Voronoi2D::generate_voronoi_2D(float * points, int num_points, void (*render)(vector<HalfEdge2D *>, float), void (*sleep)()) {
    
    bool should_render = true;
    
    X0 = Y0 = 1;
    X1 = Y1 = 0;
    
    beach_head = NULL;
    //voronoi_verticies.clear();
    voronoi_edges.clear();
    while (!circle_event_queue.empty()) {
        circle_event_queue.pop();
    }
    sweep_line = 0;
    
    priority_queue<Point2D, vector<Point2D>, PriorityQueueCompare> site_event_queue;
    
    for (int i = 0; i < num_points; i++) {
        
        if (points[2 * i] < X0) {X0 = points[2 * i];}
        if (points[2 * i] > X1) {X1 = points[2 * i];}
        if (points[2 * i + 1] < Y0) {Y0 = points[2 * i + 1];}
        if (points[2 * i + 1] > Y1) {Y1 = points[2 * i + 1];}
        
        Point2D site = Point2D(points[2 * i], points[2 * i + 1]);
        site_event_queue.push(site);
    }
    
    // Add margins to the bounding box.
    float dx = (X1-X0+1)/5.0, dy = (Y1-Y0+1)/5.0;
    X0 -= dx;  X1 += dx;  Y0 -= dy;  Y1 += dy;
    
    //cout << X0 << ", " << X1 << ", " << Y0 << ", " << Y1 << "\n";
    
    while (!site_event_queue.empty() || !circle_event_queue.empty()) {
        
        //cout << "[" << site_event_queue.size() << ", " << circle_event_queue.size() << "]\n";
        
        if (site_event_queue.empty() || (!circle_event_queue.empty() && circle_event_queue.top()->right_most_x < site_event_queue.top().x)) {
            
            CircleEvent2D * circle = circle_event_queue.top();
            sweep_line = circle->right_most_x;
            handle_circle_event(circle);
            circle_event_queue.pop();
        }
        else {
            Point2D site = site_event_queue.top();
            sweep_line = site.x;
            handle_site_event(site);
            site_event_queue.pop();
        }
        
        if (should_render) {
            render(voronoi_edges, sweep_line);
            sleep();
        }
        
    }
    
    sweep_line = 2 * (X1 + (X1 - X0) + (Y1 - Y0));
    
    for (Arc2D * arc = beach_head; arc->next != NULL; arc = arc->next) {
        if (arc->s1 != NULL) {
            //cout << arc->location << arc->next->location << sweep_line << "\n";
            arc->s1->finish(parabolic_intersection(arc->location, arc->next->location));
            
            if (should_render) {
                render(voronoi_edges, sweep_line);
                sleep();
            }
        }
    }
}

void Voronoi2D::handle_site_event(Point2D site) {
    
    //cout << "Handle site event " << site;
    
    if (beach_head == NULL) {
        beach_head = new Arc2D(site);
        return;
    }
    
    Arc2D * i;
    
    Point2D z, zz;
    
    for (i = beach_head; i != NULL; i = i->next) {
        
        //cout << "Looking at " << i->location;
        
        if (does_intersect(site, i, &z)) {
            
            if (i->next != NULL && !does_intersect(site, i->next, &zz)) {
                //cout << "Split arc.\n";
                i->next->prev = new Arc2D(i->location, i, i->next);
                i->next = i->next->prev;
            }
            else {//if (i->next == NULL) {
                //cout << "Inbetween arcs.\n";
                i->next = new Arc2D(i->location, i);
            }
            
            i->next->s1 = i->s1;
            
            i->next->prev = new Arc2D(site, i, i->next);
            i->next = i->next->prev;
            
            i = i->next;
            
            i->prev->s1 = i->s0 = new HalfEdge2D(z);
            voronoi_edges.push_back(i->s0);
            i->next->s0 = i->s1 = new HalfEdge2D(z);
            voronoi_edges.push_back(i->s1);
            
            check_circle_event(i);
            check_circle_event(i->prev);
            check_circle_event(i->next);
            
            return;
        }
        
    }
    
    for(i = beach_head; i->next != NULL; i = i->next){}
    
    //we append arc to end
    i->next = new Arc2D(site, i);
    
    Point2D start(X0, (i->location.y + i->next->location.y) / 2.f);
    i->s1 = i->next->s0 = new HalfEdge2D(start);
    voronoi_edges.push_back(i->s1);
}

void Voronoi2D::handle_circle_event(CircleEvent2D * event) {
    
    //cout << "Handle circle event." << event->circumcenter;
    
    if (event->is_valid) {
        
        HalfEdge2D * edge = new HalfEdge2D(event->circumcenter);
        voronoi_edges.push_back(edge);
        
        Arc2D * a = event->arc;
        
        if (a->prev != NULL) {
            a->prev->next = a->next;
            a->prev->s1 = edge;
        }
        if (a->next != NULL) {
            a->next->prev = a->prev;
            a->next->s0 = edge;
        }
        
        if (a->s0 != NULL) {
            a->s0->finish(event->circumcenter);
        }
        if (a->s1 != NULL) {
            a->s1->finish(event->circumcenter);
        }
        
        if (a->prev != NULL) {
            check_circle_event(a->prev);
        }
        if (a->next != NULL) {
            check_circle_event(a->next);
        }
        
    }
    delete event;
}

void Voronoi2D::check_circle_event(Arc2D * arc) {
    
    //cout << "Check circle event.\n";
    
    if (arc->event != NULL && arc->event->right_most_x != sweep_line) {
        arc->event->is_valid = false;
    }
    arc->event = NULL;
    
    if (arc->prev == NULL || arc->next == NULL) {
        return;
    }
    
    float x;
    Point2D circumcenter;
    
    if (make_circle(arc->prev->location, arc->location, arc->next->location, circumcenter, x) && x > sweep_line) {
        
        arc->event = new CircleEvent2D(arc, circumcenter, x);
        circle_event_queue.push(arc->event);
    }
    
}

bool Voronoi2D::make_circle(Point2D a, Point2D b, Point2D c, Point2D & circumcenter, float & right_most_x) {
    
    //cout << "Make circle.\n" << a << b << c;
    
    // Check that bc is a "right turn" from ab.
    if ((b.x - a.x) * (c.y - a.y) - (c.x - a.x) * (b.y - a.y) > 0) {
        //cout << "Inverted circle.\n";
        return false;
    }
    
    // Algorithm from O'Rourke 2ed p. 189.
    float A = b.x - a.x,  B = b.y - a.y,
    C = c.x - a.x,  D = c.y - a.y,
    E = A*(a.x+b.x) + B*(a.y+b.y),
    F = C*(a.x+c.x) + D*(a.y+c.y),
    G = 2*(A*(c.y-b.y) - B*(c.x-b.x));
    
    // Points are co-linear.
    if (G == 0) {
        //cout << "Co-linear.\n";
        
        /*
        *right_most_x = INF;
        
        circumcenter->x = (a.y - b.y) * INF + a.x;
        circumcenter->y = (b.x - a.x) * INF + a.y;
        
        cout << *circumcenter;
        
        return true;*/
        return false;
    }
    
    // Point o is the center of the circle.
    circumcenter.x = (D*E-B*F)/G;
    circumcenter.y = (A*F-C*E)/G;
    
    // o.x plus radius equals max x coordinate.
    right_most_x = circumcenter.x + sqrt( pow(a.x - circumcenter.x, 2) + pow(a.y - circumcenter.y, 2) );
    
    return true;
}


// Will a new parabola at point p intersect with arc i?
bool Voronoi2D::does_intersect(Point2D point, Arc2D * arc, Point2D * ret) {
    
    if (point.x == arc->location.x) {return false;}
    
    float down_y = 0;
    float up_y = 0;
    
    if (arc->prev != NULL) {
        down_y = parabolic_intersection(arc->prev->location, arc->location).y;
    }
    if (arc->next != NULL) {
        up_y = parabolic_intersection(arc->location, arc->next->location).y;
    }
    if ((arc->prev == NULL || down_y <= point.y) && (arc->next == NULL || up_y >= point.y)) {
        
        ret->y = point.y;
        
        ret->x = (arc->location.x * arc->location.x + (arc->location.y - ret->y) * (arc->location.y - ret->y) - point.x * point.x) / (2 * arc->location.x - 2 * point.x);
        
        return true;
    }
    
    
    return false;
}

// Where do two parabolas intersect?
Point2D Voronoi2D::parabolic_intersection(Point2D left, Point2D right) {
    
    Point2D ret;
    
    Point2D p = left;
    
    if (left.y == right.y) {
        return (left < right) ? left : right;
    }
    if (left.x == right.x) {
        ret.y = (left.y + right.y) / 2.f;
    }
    else if (right.x == sweep_line) {
        ret.y = right.y;
    }
    else if (left.x == sweep_line) {
        ret.y = left.y;
        p = right;
    }
    else {
        float z0 = 2 * (left.x - sweep_line);
        float z1 = 2 * (right.x - sweep_line);
        
        float a = 1 / z0 - 1 / z1;
        float b = -2 * (left.y / z0 - right.y / z1);
        float c = (left.y * left.y + left.x * left.x - sweep_line * sweep_line) / z0 - (right.y * right.y + right.x * right.x - sweep_line * sweep_line) / z1;
        
        ret.y = (-b - sqrt(b*b - 4*a*c)) / (2 * a);
    }
    
    ret.x = (p.x * p.x + (p.y - ret.y) * (p.y - ret.y) - sweep_line * sweep_line) / (2 * p.x - 2 * sweep_line);
    
    //cout << ret;
    
    return ret;
}



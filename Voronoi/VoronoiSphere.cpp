//
//  VoronoiSphere.cpp
//  Voronoi
//
//  Created by Ellis Sparky Hoag on 6/3/16.
//  Copyright © 2016 Ellis Sparky Hoag. All rights reserved.
//

#include "VoronoiSphere.hpp"

using namespace std;
using namespace Voronoi;

VoronoiDiagramSphere VoronoiSphere::generate_voronoi(float * points, int num_points, void (*render)(vector<HalfEdgeSphere *>, float), void (*sleep)()) {
    
    bool should_render = false;
    
    beach_head = NULL;
    sweep_line = 0;
    voronoi_diagram.edges.clear();
    while (!circle_event_queue.empty()) {
        circle_event_queue.pop();
    }
    
    //priority_queue<PointSphere, vector<PointSphere>, PriorityQueueCompare> site_event_queue;
    priority_queue<VoronoiCellSphere *, vector<VoronoiCellSphere *>, PriorityQueueCompare> site_event_queue;

    for (int i = 0; i < num_points; i++) {
        
        PointCartesian point_cartesian = PointCartesian(points[3 * i], points[3 * i + 1], points[3 * i + 2]);
        
        VoronoiCellSphere * cell = new VoronoiCellSphere(point_cartesian);
        
        //site_event_queue.push(PointSphere(point_cartesian));
        site_event_queue.push(cell);
        
        voronoi_diagram.cells.push_back(cell);
        
    }
    
    while (!site_event_queue.empty() || !circle_event_queue.empty()) {
        
        //cout << "[" << site_event_queue.size() << ", " << circle_event_queue.size() << "]\n";
        
        if (site_event_queue.empty() || (!circle_event_queue.empty() && circle_event_queue.top()->lowest_theta < site_event_queue.top()->site.theta)) {
            
            CircleEventSphere * circle = circle_event_queue.top();
            sweep_line = circle->lowest_theta;
            handle_circle_event(circle);
            circle_event_queue.pop();
            
            if (should_render) {
                render(voronoi_diagram.edges, sweep_line);
                ArcSphere * cur = beach_head;
                do {
                    cout << cur->cell->site;
                    cur = cur->next[0];
                } while (cur != beach_head);
                cout << endl;
                sleep();
            }
        }
        else {
            
            VoronoiCellSphere * cell = site_event_queue.top();
            sweep_line = cell->site.theta;
            handle_site_event(cell);
            site_event_queue.pop();
            
            if (should_render) {
                render(voronoi_diagram.edges, sweep_line);
                ArcSphere * cur = beach_head;
                do {
                    cout << cur->cell->site;
                    cur = cur->next[0];
                } while (cur != beach_head);
                cout << endl;
                sleep();
            }
        }
    }
    
    
    /*sweep_line = 2 * M_PI;
    ArcSphere * arc = beach_head;
    do {
        cout << arc->location;
        PointSphere intersection;
        parabolic_intersection(arc->location, arc->next[0]->location, intersection);
        arc->s1->finish(intersection);
        
        arc = arc->next[0];
    } while (arc != beach_head);
    */
    
    //finalize_diagram();
    
    return voronoi_diagram;
}

void VoronoiSphere::finalize_diagram() {
    
    //last two arcs in the beach need to be connected?
    
    for (auto cell : voronoi_diagram.cells) {
        HalfEdgeSphere * left = NULL;
        for (auto edge : cell->edges) {
            if (!edge->is_finished) {
                if (left == NULL) {
                    left = edge;
                }
                else {
                    if ((left->cells[0] == edge->cells[0] && left->cells[1] == edge->cells[1]) || (left->cells[0] == edge->cells[1] && left->cells[1] == edge->cells[0])) {
                        cout << "fixed\n";
                        left->finish(edge->start);
                        edge->finish(left->start);
                    }
                }
            }
        }
    }
    
}

void VoronoiSphere::handle_site_event(VoronoiCellSphere * cell) {
    
    //cout << "Handle site event.\n";
    
    if (beach_head == NULL) {
        beach_head = new ArcSphere(cell);
        return;
    }
    if (beach_head->next[0] == beach_head) {
        new ArcSphere(cell, beach_head, beach_head);
        beach_head->next[0]->s0 = beach_head->s1;
        beach_head->next[0]->s1 = beach_head->s0;
        return;
    }
    
    ArcSphere * cur = beach_head;
    
    int level = cur->height - 1;

    do {
        
        PointSphere left_intersection, right_intersection;
        
        float phi_start, phi_end;
        
        if (cur->prev[0] != cur && parabolic_intersection(cur->prev[0]->cell->site, cur->cell->site, left_intersection)) {
            phi_start = left_intersection.phi;
        }
        else {
            phi_start = cur->cell->site.phi - M_PI;
        }
        if (cur != cur->next[0] && parabolic_intersection(cur->cell->site, cur->next[0]->cell->site, right_intersection)) {
            phi_end = right_intersection.phi;
        }
        else {
            phi_end = cur->cell->site.phi + M_PI;
        }
        
        if ((phi_start < phi_end && phi_start <= cell->site.phi && cell->site.phi <= phi_end) || (phi_start > phi_end && (phi_start < cell->site.phi || cell->site.phi < phi_end)) ) {
            //arc is found
            
            if (cur->event != NULL) {
                cur->event->is_valid = false;
            }
            
            //duplicate
            new ArcSphere(cur->cell, cur, cur->next[0]);
            cur->next[0]->s1 = cur->s1;

            //insert new site
            new ArcSphere(cell, cur, cur->next[0]);
            
            cur = cur->next[0];
            //cur is now the new arc in between two identical arcs
            
            //insert new vertex
            // This is not really a vertex. It is most likely in the center of some edge.
            PointCartesian * vertex = new PointCartesian();
            *vertex = phi_to_point(cur->prev[0]->cell->site, cell->site.phi).get_cartesian();
            
            voronoi_diagram.verticies.push_back(vertex);
            
            HalfEdgeSphere * left_edge = new HalfEdgeSphere(vertex, cur->cell, cur->prev[0]->cell);
            HalfEdgeSphere * right_edge = new HalfEdgeSphere(vertex, cur->cell, cur->prev[0]->cell);
            
            voronoi_diagram.edges.push_back(left_edge);
            voronoi_diagram.edges.push_back(right_edge);
            
            //add new edges
            cur->prev[0]->add_edge_right(left_edge);
            cur->add_edge_left(left_edge);
            cur->add_edge_right(right_edge);
            cur->next[0]->add_edge_left(right_edge);
            
            //check for new circle events
            check_circle_event(cur->prev[0]);
            check_circle_event(cur->next[0]);
            
            return;
        }
        
        //cur = cur->next[0];
        
        //if (cur->location.phi < site.phi) {
        //    cur = cur->next[level - 1];
        //}
        //else {
        //    cur = cur->prev[level - 1];
        //}
        cur = cur->next[level];
        level = max(0, level - 1);
        
    } while (true);
}

void VoronoiSphere::handle_circle_event(CircleEventSphere * event) {
    
    if (event->is_valid) {
    
        //cout << "Handle circle event " << event->arc->cell->site;

        ArcSphere * left = event->arc->prev[0];
        ArcSphere * right = event->arc->next[0];
        
        //add new vertex
        PointCartesian * vertex = new PointCartesian();
        *vertex = event->circumcenter.get_cartesian();
        
        voronoi_diagram.verticies.push_back(vertex);
        
        //add new edge
        HalfEdgeSphere * edge = new HalfEdgeSphere(vertex, left->cell, right->cell);
        voronoi_diagram.edges.push_back(edge);
        
        left->add_edge_right(edge);
        right->add_edge_left(edge);

        //finish old edges
        if (event->arc->s0 != NULL) {
            event->arc->s0->finish(vertex);
        }
        if (event->arc->s1 != NULL) {
            event->arc->s1->finish(vertex);
        }
        
        //remove from beach;
        if (beach_head == event->arc) {
            beach_head = event->arc->next[0];
        }
        delete event->arc;
        
        //invalidate old circle events
        if (left->event != NULL) {
            left->event->is_valid = false;
        }
        if (right->event != NULL) {
            right->event->is_valid = false;
        }
        
        //check for new circle events
        check_circle_event(left);
        check_circle_event(right);
    }
    
    delete event;
}


void VoronoiSphere::check_circle_event(ArcSphere * arc) {
    
    if (arc == NULL || arc->prev[0] == NULL || arc->next[0] == NULL || arc->prev[0] == arc->next[0] || arc == arc->next[0] || arc->prev[0] == arc) {
        return;
    }
    
    //arc->event = NULL;
    
    PointSphere center;
    float lowest_theta;
    
    make_circle(arc->prev[0]->cell->site, arc->cell->site, arc->next[0]->cell->site, center, lowest_theta);
    
    if (lowest_theta > sweep_line) {
    
        arc->event = new CircleEventSphere(arc, center, lowest_theta);
        circle_event_queue.push(arc->event);
    }
}

void VoronoiSphere::make_circle(PointSphere a, PointSphere b, PointSphere c, PointSphere & circumcenter, float & lowest_theta) {
    
    PointCartesian i = a.get_cartesian();
    PointCartesian j = b.get_cartesian();
    PointCartesian k = c.get_cartesian();
    
    PointCartesian center = PointCartesian::cross_product(i - j, k - j);
    center.normalize();
    
    center.to_unit_sphere(circumcenter.theta, circumcenter.phi);
    
    float radius = PointSphere::angle_between(circumcenter, a);
    
    lowest_theta = circumcenter.theta + radius;
}


bool VoronoiSphere::parabolic_intersection(PointSphere left, PointSphere right, PointSphere & intersection) {
    
    if (left.theta >= sweep_line && right.theta >= sweep_line) {
        return false;
    }
    else if (left.theta >= sweep_line && right.theta < sweep_line) {
        intersection = phi_to_point(right, left.phi);
        return true;
    }
    else if (right.theta >= sweep_line && left.theta < sweep_line) {
        intersection = phi_to_point(left, right.phi);
        return true;
    }
    
    float cos_sweep_line = cos(sweep_line);
    
    float cos_left_theta = cos(left.theta);
    float cos_right_theta = cos(right.theta);
    
    float cos_minus_cos_right = cos_sweep_line - cos_right_theta;
    float cos_minus_cos_left = cos_sweep_line - cos_left_theta;
    
    PointCartesian u = left.get_cartesian();
    PointCartesian v = right.get_cartesian();
    
    float a = cos_minus_cos_right * u.x - cos_minus_cos_left * v.x;
    float b = cos_minus_cos_right * u.y - cos_minus_cos_left * v.y;
    
    float e = (cos_left_theta - cos_right_theta) * sin(sweep_line);
    
    float sqrt_a_b = sqrt(a*a + b*b);
    
    if (fabs(e) > sqrt_a_b) {
        // out of domain of sin
        return false;
    }
    
    float sin_phi_int_plus_gamma = e / sqrt_a_b;
    
    float gamma = atan2(a, b);
    
    float phi_int = asin(sin_phi_int_plus_gamma) - gamma;
    
    intersection = phi_to_point(left, phi_int);
    
    if (intersection.phi >= M_PI) {
        intersection.phi -= 2 * M_PI;
    }
    else if (intersection.phi <= -M_PI) {
        intersection.phi += 2 * M_PI;
    }
    return true;
}

PointSphere VoronoiSphere::phi_to_point(PointSphere arc, float phi) {
    
    if (arc.theta >= sweep_line) {
        return PointSphere(sweep_line, phi);
    }
    float a = cos(sweep_line) - cos(arc.theta);
    float b = sin(arc.theta) * cos(phi - arc.phi) - sin(sweep_line);
    return PointSphere(atan2(-a, -b), phi);
}


//
//  VoronoiSphere.cpp
//  Voronoi
//
//  Created by Ellis Sparky Hoag on 6/3/16.
//  Copyright © 2016 Ellis Sparky Hoag. All rights reserved.
//

#include "VoronoiSphere.h"

using namespace std;
using namespace Voronoi;

VoronoiDiagramSphere VoronoiSphere::generate_voronoi(float * points, int num_points, void (*render)(VoronoiDiagramSphere, float), bool (*is_sleeping)())
{
    //bool should_render = false && render != NULL && is_sleeping != NULL;
    
    initialize_diagram(points, num_points);

    while (!site_event_queue.empty() || !circle_event_queue.empty())
    {
        //cout << "[" << site_event_queue.size() << ", " << circle_event_queue.size() << "]\n";
        
        if (!circle_event_queue.empty() && !circle_event_queue.top()->is_valid)
        {
            circle_event_queue.pop();
        }
        else if (site_event_queue.empty() || (!circle_event_queue.empty() && circle_event_queue.top()->lowest_theta < site_event_queue.top().site.theta))
        {
            
            CircleEventSphere * circle = circle_event_queue.top();
            sweep_line = circle->lowest_theta;
            cos_sweep_line = cos(sweep_line);
            sin_sweep_line = sin(sweep_line);
            handle_circle_event(circle);
            circle_event_queue.pop();
            delete circle;
            
//            while (should_render && is_sleeping())
//            {
//                render(voronoi_diagram, sweep_line);
//            }
        }
        else
        {
            VoronoiCellSphere cell = site_event_queue.top();
            sweep_line = cell.site.theta;
            cos_sweep_line = cos(sweep_line);
            sin_sweep_line = sin(sweep_line);
            handle_site_event(cell);
            site_event_queue.pop();
            
//            while (should_render && is_sleeping())
//            {
//                render(voronoi_diagram, sweep_line);
//            }
        }
    }
    
    finalize_diagram();
    
    return voronoi_diagram;
}

void VoronoiSphere::initialize_diagram(float * points, int num_points)
{    
    beach_head = NULL;
    sweep_line = 0;
    voronoi_diagram.edges.clear();
    voronoi_diagram.cells.clear();
    while (!circle_event_queue.empty())
    {
        circle_event_queue.pop();
    }
    while (!site_event_queue.empty())
    {
        site_event_queue.pop();
    }
    
    for (int i = 0; i < num_points; i++)
    {
        PointCartesian point_cartesian = PointCartesian(points[3 * i], points[3 * i + 1], points[3 * i + 2]);
        
        VoronoiCellSphere cell = VoronoiCellSphere(point_cartesian, (int)voronoi_diagram.cells.size());
        
        site_event_queue.push(cell);
        
        voronoi_diagram.cells.push_back(cell);
    }
}

void VoronoiSphere::finalize_diagram()
{
    // The final two edges need to be connected.
    
    for (int i = 0; i < voronoi_diagram.edges.size(); i++)
    {
        if (!voronoi_diagram.edges[i].is_finished)
        {
            for (int k = i + 1; k < voronoi_diagram.edges.size(); k++)
            {
                if (!voronoi_diagram.edges[k].is_finished)
                {
                    voronoi_diagram.edges[i].finish(voronoi_diagram.edges[k].start);
                    voronoi_diagram.edges[k].finish(voronoi_diagram.edges[i].start);
                    return;
                }
            }
        }
    }
}

void VoronoiSphere::handle_site_event(VoronoiCellSphere cell)
{
    //cout << "Handle site event.\n";
    
    if (beach_head == NULL)
    {
        add_initial_arc_sphere(cell.cell_id);
        return;
    }
    if (beach_head == beach_head->next[0])
    {
        add_arc_sphere(cell.cell_id, beach_head, beach_head);
        
        beach_head->next[0]->left_edge_id = beach_head->right_edge_id;
        beach_head->next[0]->right_edge_id = beach_head->left_edge_id;
        
        PointSphere cur_site = voronoi_diagram.cells[beach_head->cell_id].site;

        // This is not really a vertex. It is in the middle of some edge.
        PointCartesian vertex = phi_to_point(cur_site, cell.site.phi).get_cartesian();
        
        add_half_edge_sphere(vertex, beach_head, beach_head->next[0]);
        add_half_edge_sphere(vertex, beach_head->prev[0], beach_head);
        
        return;
    }
    
    bool move_right = true;
    
    ArcSphere * arc = traverse_skiplist_to_site(beach_head, cell.site.phi);
    
    ArcSphere * left = arc->prev[0];
    ArcSphere * right = arc->next[0];
    
    while (true)
    {
        PointSphere cur_site = voronoi_diagram.cells[arc->cell_id].site;
        PointSphere prev_site = voronoi_diagram.cells[arc->prev[0]->cell_id].site;
        PointSphere next_site = voronoi_diagram.cells[arc->next[0]->cell_id].site;

        Real phi_start, phi_end;
        
        bool valid_arc = parabolic_intersection(prev_site, cur_site, phi_start) && parabolic_intersection(cur_site, next_site, phi_end);
        
        if (!valid_arc && ((prev_site.phi < next_site.phi && prev_site.phi <= cell.site.phi && cell.site.phi <= next_site.phi) || (prev_site.phi > next_site.phi && (prev_site.phi <= cell.site.phi || cell.site.phi <= next_site.phi))))
        {
            /*
             *  Here we have two sites that are probably on the sweep line. The arcs only intersect at the north pole.
             *  Also, our query site is in between them. So it only needs to be inserted in the beachline.
             *  This is extreamly unlikely and has not happened yet. But we account for it anyway.
             */
            
            add_arc_sphere(cell.cell_id, arc, arc->next[0]);
            cout << "Two arcs interesect at the north pole.\n";
            //assert(0);
            return;
        }
        else if (valid_arc && ((phi_start < phi_end && phi_start <= cell.site.phi && cell.site.phi <= phi_end) || (phi_start > phi_end && (phi_start <= cell.site.phi || cell.site.phi <= phi_end))))
        {
            // The arc is found!
            //cout << phi_start << ", " << phi_end << endl;
            //cout << voronoi_diagram.cells[arc->cell_id].site << cell.site << voronoi_diagram.cells[arc->next[0]->cell_id].site << sweep_line << "\n\n";
            
            if (arc->event != NULL)
            {
                arc->event->is_valid = false;
            }
            
            //duplicate arc
            add_arc_sphere(arc->cell_id, arc, arc->next[0]);
            arc->next[0]->right_edge_id = arc->right_edge_id;

            //insert new site
            add_arc_sphere(cell.cell_id, arc, arc->next[0]);
            
            // This is not really a vertex. It is in the middle of some edge.
            PointCartesian vertex = phi_to_point(cur_site, cell.site.phi).get_cartesian();
            
            add_half_edge_sphere(vertex, arc, arc->next[0]);
            add_half_edge_sphere(vertex, arc->next[0], arc->next[0]->next[0]);
            
            //check for new circle events
            check_circle_event(arc);
            check_circle_event(arc->next[0]->next[0]);
            
            return;
        }
        
        if (move_right)
        {
            arc = right;
            right = right->next[0];
        }
        else
        {
            arc = left;
            left = left->prev[0];
        }
        
        move_right = !move_right;
    }
}

void VoronoiSphere::handle_circle_event(CircleEventSphere * event)
{
    //cout << "Handle circle event " << event->circumcenter;

    ArcSphere * left = event->arc->prev[0];
    ArcSphere * right = event->arc->next[0];
    
    //add new vertex
    PointCartesian vertex = event->circumcenter.get_cartesian();
    
    //voronoi_diagram.verticies.push_back(vertex);
    
    //add new edge
    add_half_edge_sphere(vertex, left, right);
    
    //finish old edges
    int left_id = event->arc->left_edge_id;
    int right_id = event->arc->right_edge_id;
    if (left_id != -1)
    {
        voronoi_diagram.edges[left_id].finish(vertex);
    }
    if (right_id != -1)
    {
        voronoi_diagram.edges[right_id].finish(vertex);
    }
    
    //invalidate old circle events
    if (left->event != NULL)
    {
        left->event->is_valid = false;
    }
    if (right->event != NULL)
    {
        right->event->is_valid = false;
    }
    
    remove_arc_sphere(event->arc);
    
    //check for new circle events
    check_circle_event(left);
    check_circle_event(right);
}

void VoronoiSphere::check_circle_event(ArcSphere * arc)
{
    if (arc == NULL || arc->prev[0] == NULL || arc->next[0] == NULL || arc->prev[0] == arc->next[0] || arc == arc->next[0] || arc->prev[0] == arc)
    {
        //cout << "Invalid circle event.\n";
        return;
    }
    
    PointSphere circumcenter;
    Real lowest_theta;
    
    PointSphere cur_site = voronoi_diagram.cells[arc->cell_id].site;
    PointSphere prev_site = voronoi_diagram.cells[arc->prev[0]->cell_id].site;
    PointSphere next_site = voronoi_diagram.cells[arc->next[0]->cell_id].site;
    
    make_circle(prev_site, cur_site, next_site, circumcenter, lowest_theta);
    
    //This if statement causes errors for some reason...
    //if (lowest_theta > sweep_line)
    {
        arc->event = new CircleEventSphere(arc, circumcenter, lowest_theta);
        circle_event_queue.push(arc->event);
    }
}

void VoronoiSphere::make_circle(PointSphere a, PointSphere b, PointSphere c, PointSphere & circumcenter, Real & lowest_theta)
{    
    PointCartesian i = a.get_cartesian();
    PointCartesian j = b.get_cartesian();
    PointCartesian k = c.get_cartesian();
    
    PointCartesian center = PointCartesian::cross_product(i - j, k - j);
    center.normalize();
    
    circumcenter = PointSphere(center);

    Real radius = PointSphere::angle_between(circumcenter, a);
    
    lowest_theta = circumcenter.theta + radius;
}

bool VoronoiSphere::parabolic_intersection(PointSphere left, PointSphere right, Real & phi_intersection)
{
    if (left.theta == sweep_line && right.theta == sweep_line)
    {
        /*
         *  Both sites are on the sweep line so our intersection is at the north pole.
         *  We need to handle this special case.
         */
        //cout << "Intersection is at north pole.\n";
        return false;
    }
    else if (left.theta == sweep_line && right.theta < sweep_line)
    {
        /*
         *  The left site is on our sweep line so it contains our intersection phi.
         */
        phi_intersection = left.phi;
        //cout << "Left site is on sweep line.\n";
        return true;
    }
    else if (right.theta == sweep_line && left.theta < sweep_line)
    {
        /*
         *  The right site is on our sweep line so it contains our intersection phi.
         *  If the beachline contains exactly two arcs then our intersection phi is
         *  on the other side of our sphere.
         */
        if (beach_head == beach_head->next[0]->next[0])
        {
            phi_intersection = right.phi + ((right.phi > 0) ? M_PI : -M_PI);
            //cout << "Right site is on sweep line and there are exactly two arcs in beach.\n";
        }
        else
        {
            phi_intersection = right.phi;
            //cout << "Right site is on sweep line.\n";
        }
        return true;
    }
    
    Real cos_left_theta = cos(left.theta);
    Real cos_right_theta = cos(right.theta);
    
    Real cos_minus_cos_right = cos_sweep_line - cos_right_theta;
    Real cos_minus_cos_left = cos_sweep_line - cos_left_theta;
    
    PointCartesian u = left.get_cartesian();
    PointCartesian v = right.get_cartesian();
    
    Real a = cos_minus_cos_right * u.x - cos_minus_cos_left * v.x;
    Real b = cos_minus_cos_right * u.y - cos_minus_cos_left * v.y;
    
    Real e = (cos_left_theta - cos_right_theta) * sin_sweep_line;
    
    Real sqrt_a_b = sqrt(a*a + b*b);
    
    if (fabs(e) > sqrt_a_b)
    {
        /*
         *  Theoretically this is impossible becasue as the sweep line
         *  approaches either the left or right theta e / sqrt_a_b approaches
         *  1. In reality, e can be greater than sqrt_a_b due to rounding.
         *  This excecutes when a site is close enought to the sweep line
         *  to cause problems. So one of the parabolas is actually a line.
         *  The lower site contains the phi of our intersection.
         */
        phi_intersection = (left.theta > right.theta) ? left.phi : right.phi;
        //cout << "Arc is close to sweep line. " << setprecision(50) << fabs(e) << " > " << sqrt_a_b << endl << left << right << sweep_line << endl;
        return true;
    }
    
    Real sin_phi_int_plus_gamma = e / sqrt_a_b;
    
    Real gamma = atan2(a, b);
    
    phi_intersection = asin(sin_phi_int_plus_gamma) - gamma;
    
    if (phi_intersection > M_PI)
    {
        phi_intersection -= 2 * M_PI;
    }
    else if (phi_intersection < -M_PI)
    {
        phi_intersection += 2 * M_PI;
    }
    return true;
}

PointSphere VoronoiSphere::phi_to_point(PointSphere arc, Real phi)
{
    Real a = cos(arc.theta) - cos_sweep_line;
    Real b = sin_sweep_line - sin(arc.theta) * cos(phi - arc.phi);
    return PointSphere(atan2(a, b), phi);
}

void VoronoiSphere::add_initial_arc_sphere(int cell_id)
{    
    int height = random_height();
    
    beach_head = new ArcSphere(cell_id, height);
    
    beach_head->prev = new ArcSphere*[height];
    beach_head->next = new ArcSphere*[height];
    
    for (int i = 0; i < height; i++)
    {
        beach_head->prev[i] = beach_head->next[i] = beach_head;
    }
}

void VoronoiSphere::add_arc_sphere(int cell_id, ArcSphere * left, ArcSphere * right)
{
    int height = random_height();
    
    ArcSphere * arc = new ArcSphere(cell_id, height);

    arc->prev = new ArcSphere*[height];
    arc->next = new ArcSphere*[height];
    
    for (int i = 0; i < height; i++)
    {
        while (left->height <= i)
        {
            left = left->prev[left->height - 1];
        }
        while (right->height <= i)
        {
            right = right->next[right->height - 1];
        }
        
        left->next[i] = arc;
        arc->prev[i] = left;
        
        right->prev[i] = arc;
        arc->next[i] = right;
    }
}

void VoronoiSphere::remove_arc_sphere(ArcSphere * arc)
{
    if (beach_head == arc)
    {
        beach_head = arc->next[0];
    }
    
    if (arc->event != NULL)
    {
        arc->event->is_valid = false;
    }
    
    for (int i = 0; i < arc->height; i++)
    {
        ArcSphere * left = arc->prev[i];
        left->next[i] = arc->next[i];
        arc->next[i]->prev[i] = left;
    }
    
    delete [] arc->prev;
    delete [] arc->next;
    
    delete arc;
}

void VoronoiSphere::add_half_edge_sphere(PointCartesian start, ArcSphere * left, ArcSphere *right)
{
    int edge_id = (int)voronoi_diagram.edges.size();
    
    HalfEdgeSphere edge(start, left->cell_id, right->cell_id);
    voronoi_diagram.edges.push_back(edge);
    
    voronoi_diagram.cells[left->cell_id].edge_ids.push_back(edge_id);
    voronoi_diagram.cells[right->cell_id].edge_ids.push_back(edge_id);
    
    left->right_edge_id = right->left_edge_id = edge_id;
}

ArcSphere * VoronoiSphere::traverse_skiplist_to_site(ArcSphere * arc, Real phi)
{
    int level = arc->height - 1;
    
    bool move_right = true;
    
    while (level >= 0)
    {
        Real delta = fabs(phi - voronoi_diagram.cells[arc->cell_id].site.phi);
        if (delta > M_PI) {delta = 2.f * M_PI - delta;}
        
        Real next_delta = fabs(phi - voronoi_diagram.cells[arc->next[level]->cell_id].site.phi);
        if (next_delta > M_PI) {delta = 2.f * M_PI - next_delta;}

        Real prev_delta = fabs(phi - voronoi_diagram.cells[arc->prev[level]->cell_id].site.phi);
        if (prev_delta > M_PI) {delta = 2.f * M_PI - prev_delta;}

        if (next_delta < delta && prev_delta > delta)
        {
            arc = arc->next[level];
            if (!move_right) {level--;}
            move_right = true;
        }
        else if (prev_delta < delta && next_delta > delta)
        {
            arc = arc->prev[level];
            if (move_right) {level--;}
            move_right = false;
        }
        else// if (next_delta >= delta && prev_delta >= delta)
        {
            level--;
        }
    }
    return arc;
}

int VoronoiSphere::random_height()
{
    int r = rand() % 32767;
    
    if (r < 16384)
    {
        return 1;
    }
    else if (r < 24576)
    {
        return 2;
    }
    else if (r < 28672)
    {
        return 3;
    }
    else if (r < 30720)
    {
        return 4;
    }
    else if (r < 31744)
    {
        return 5;
    }
    else if (r < 32256)
    {
        return 6;
    }
    else if (r < 32512)
    {
        return 7;
    }
    else if (r < 32640)
    {
        return 8;
    }
    else if (r < 32704)
    {
        return 9;
    }
    else if (r < 32736)
    {
        return 10;
    }
    else if (r < 32752)
    {
        return 11;
    }
    else if (r < 32760)
    {
        return 12;
    }
    else if (r < 32764)
    {
        return 13;
    }
    else if (r < 32766)
    {
        return 14;
    }
    return 15;
}
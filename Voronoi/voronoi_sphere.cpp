//
//  VoronoiSphere.cpp
//  Voronoi
//
//  Created by Ellis Sparky Hoag on 6/3/16.
//  Copyright © 2016 Ellis Sparky Hoag. All rights reserved.
//

#include "voronoi_sphere.h"

using namespace std;

namespace Voronoi {

    VoronoiDiagramSphere generate_voronoi_parallelized(vector<tuple<float, float, float>> * verts)
    {
        VoronoiDiagramSphere diagram_top_down, diagram_bottom_up;
     
        thread diagram_top_down_thread(compute_priority_queues, &diagram_top_down, verts, true);
        
        compute_priority_queues(&diagram_bottom_up, verts, false);

        diagram_top_down_thread.join();
        
        unsigned long site_length = diagram_top_down.sites.size();
        unsigned long verticies_length = diagram_top_down.voronoi_verticies.size();
        
        for (auto site : diagram_bottom_up.sites)
        {
            diagram_top_down.sites.push_back(site);
        }
        
        for (auto vertex : diagram_bottom_up.voronoi_verticies)
        {
            diagram_top_down.voronoi_verticies.push_back(vertex);
        }
        
        for (auto voronoi_edge : diagram_bottom_up.voronoi_edges)
        {
            voronoi_edge.vidx[0] += verticies_length;
            voronoi_edge.vidx[1] += verticies_length;
            
            diagram_top_down.voronoi_edges.push_back(voronoi_edge);
        }
        
        for (auto delaunay_edge : diagram_bottom_up.delaunay_edges)
        {
            delaunay_edge.vidx[0] += site_length;
            delaunay_edge.vidx[1] += site_length;
            
            diagram_top_down.delaunay_edges.push_back(delaunay_edge);
        }
        
        return diagram_top_down;
        
    }

    void compute_priority_queues(VoronoiDiagramSphere * voronoi_diagram, vector<tuple<float, float, float>> * verts, bool is_top_down)
    {
        ArcSphere * beach_head = NULL;
        
        Real sweep_line = 0;
                
        vector<HalfEdgeSphere> half_edges;
        
        vector<VoronoiCellSphere> cells;

        priority_queue<VoronoiCellSphere, vector<VoronoiCellSphere>, PriorityQueueCompare> site_event_queue;

        priority_queue<CircleEventSphere *, vector<CircleEventSphere *>, PriorityQueueCompare> circle_event_queue;
        
        for (auto point : *verts)
        {
            Real z_value = get<2>(point);
            
            if ((is_top_down && z_value >= 0) || (!is_top_down && z_value < 0))
            {
                PointCartesian point_cartesian = PointCartesian(get<0>(point), get<1>(point), z_value);

                if (!is_top_down) {point_cartesian.z *= -1;}
       
                VoronoiCellSphere cell = VoronoiCellSphere(point_cartesian, (int)cells.size());
                
                site_event_queue.push(cell);
                
                cells.push_back(cell);
                
                voronoi_diagram->sites.push_back(point_cartesian);
            }
        }
        
        while (sweep_line <= M_PI_2 && (!site_event_queue.empty() || !circle_event_queue.empty()))
        {
            if (!circle_event_queue.empty() && !circle_event_queue.top()->is_valid)
            {
                circle_event_queue.pop();
            }
            else if (site_event_queue.empty() || (!circle_event_queue.empty() && site_event_queue.top().site.theta > circle_event_queue.top()->lowest_theta))
            {
                CircleEventSphere * circle = circle_event_queue.top();
                sweep_line = circle->lowest_theta;
                handle_circle_event(circle, voronoi_diagram, &cells, &half_edges, &circle_event_queue, beach_head);
                circle_event_queue.pop();
                delete circle;
            }
            else
            {
                VoronoiCellSphere cell = site_event_queue.top();
                sweep_line = cell.site.theta;
                handle_site_event(cell, voronoi_diagram, &cells, &half_edges, &circle_event_queue, beach_head, sweep_line, sin(sweep_line), cos(sweep_line));
                site_event_queue.pop();
            }
        }
        
        if (!is_top_down)
        {
            for (int i = 0; i < voronoi_diagram->sites.size(); i++)
            {
                voronoi_diagram->sites[i].z *= -1;
            }
            for (int i = 0; i < voronoi_diagram->voronoi_verticies.size(); i++)
            {
                voronoi_diagram->voronoi_verticies[i].z *= -1;
            }
        }
    }

    VoronoiDiagramSphere generate_voronoi(vector<tuple<float, float, float>> * verts, void (*render)(VoronoiDiagramSphere, float), bool (*is_sleeping)())
    {
        bool should_render = (false) && render != NULL && is_sleeping != NULL;

        VoronoiDiagramSphere voronoi_diagram;
        
        ArcSphere * beach_head = NULL;
        
        Real sweep_line;
        
        vector<HalfEdgeSphere> half_edges;
        
        vector<VoronoiCellSphere> cells;

        priority_queue<VoronoiCellSphere, vector<VoronoiCellSphere>, PriorityQueueCompare> site_event_queue;
        
        priority_queue<CircleEventSphere *, vector<CircleEventSphere *>, PriorityQueueCompare> circle_event_queue;
        
        for (auto point : *verts)
        {
            PointCartesian point_cartesian = PointCartesian(get<0>(point), get<1>(point), get<2>(point));
            
            VoronoiCellSphere cell = VoronoiCellSphere(point_cartesian, (int)cells.size());
            
            site_event_queue.push(cell);
            
            cells.push_back(cell);
            
            voronoi_diagram.sites.push_back(point_cartesian);
        }
        
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
                handle_circle_event(circle, &voronoi_diagram, &cells, &half_edges, &circle_event_queue, beach_head);
                circle_event_queue.pop();
                delete circle;
                
                while (should_render && is_sleeping())
                {
                    render(voronoi_diagram, sweep_line);
                }
            }
            else
            {
                VoronoiCellSphere cell = site_event_queue.top();
                sweep_line = cell.site.theta;
                handle_site_event(cell, &voronoi_diagram, &cells, &half_edges, &circle_event_queue, beach_head, sweep_line, sin(sweep_line), cos(sweep_line));
                site_event_queue.pop();
                
                while (should_render && is_sleeping())
                {
                    render(voronoi_diagram, sweep_line);
                }
            }
        }
        
        // The final two edges need to be connected.
        for (int i = 0; i < half_edges.size(); i++)
        {
            if (!half_edges[i].is_finished)
            {
                for (int k = i + 1; k < half_edges.size(); k++)
                {
                    if (!half_edges[k].is_finished)
                    {
                        finish_half_edge_sphere(&voronoi_diagram, &half_edges, i, voronoi_diagram.voronoi_verticies[half_edges[k].start_idx]);
                        finish_half_edge_sphere(&voronoi_diagram, &half_edges, k, voronoi_diagram.voronoi_verticies[half_edges[i].start_idx]);
                        
                        return voronoi_diagram;
                    }
                }
            }
        }
        
        //This should never happen
        assert(0);
        return VoronoiDiagramSphere();
    }

    void handle_site_event(VoronoiCellSphere cell, VoronoiDiagramSphere * voronoi_diagram, vector<VoronoiCellSphere> * cells, vector<HalfEdgeSphere> * half_edges, priority_queue<CircleEventSphere *, vector<CircleEventSphere *>, PriorityQueueCompare> * circle_event_queue_ptr, ArcSphere * & beach_head, Real sweep_line, Real sin_sweep_line, Real cos_sweep_line)
    {
        //cout << "Handle site event.\n" << cell.site;
        
        if (beach_head == NULL)
        {
            add_initial_arc_sphere(cell.cell_idx, beach_head);
            return;
        }
        if (beach_head == beach_head->next[0])
        {
            add_arc_sphere(cell.cell_idx, beach_head, beach_head, beach_head);
            
            beach_head->next[0]->left_edge_idx = beach_head->right_edge_idx;
            beach_head->next[0]->right_edge_idx = beach_head->left_edge_idx;
            
            PointSphere cur_site = (*cells)[beach_head->cell_idx].site;

            // This is not really a vertex. It is in the middle of some edge.
            PointCartesian vertex = phi_to_point(cur_site, cell.site.phi, sweep_line, sin_sweep_line, cos_sweep_line).get_cartesian();
            
            add_half_edge_sphere(voronoi_diagram, cells, half_edges, vertex, beach_head, beach_head->next[0]);
            add_half_edge_sphere(voronoi_diagram, cells, half_edges, vertex, beach_head->prev[0], beach_head);
            
            return;
        }
        
        bool move_right = true;
        
        ArcSphere * arc = traverse_skiplist_to_site(beach_head, cell.site.phi, cells);
        
        ArcSphere * left = arc->prev[0];
        ArcSphere * right = arc->next[0];
        
        while (true)
        {
            PointSphere cur_site = (*cells)[arc->cell_idx].site;
            PointSphere prev_site = (*cells)[arc->prev[0]->cell_idx].site;
            PointSphere next_site = (*cells)[arc->next[0]->cell_idx].site;

            Real phi_start, phi_end;
            
            bool valid_arc = parabolic_intersection(prev_site, cur_site, phi_start, beach_head, sweep_line, sin_sweep_line, cos_sweep_line) && parabolic_intersection(cur_site, next_site, phi_end, beach_head, sweep_line, sin_sweep_line, cos_sweep_line);
            
            if (!valid_arc && ((prev_site.phi < next_site.phi && prev_site.phi <= cell.site.phi && cell.site.phi <= next_site.phi) || (prev_site.phi > next_site.phi && (prev_site.phi <= cell.site.phi || cell.site.phi <= next_site.phi))))
            {
                /*
                 *  Here we have two sites that are probably on the sweep line. The arcs only intersect at the north pole.
                 *  Also, our query site is in between them. So it only needs to be inserted in the beachline.
                 *  This is extreamly unlikely and has not happened yet. But we account for it anyway.
                 */
                
                add_arc_sphere(cell.cell_idx, arc, arc->next[0], beach_head);
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
                add_arc_sphere(arc->cell_idx, arc, arc->next[0], beach_head);
                arc->next[0]->right_edge_idx = arc->right_edge_idx;

                //insert new site
                add_arc_sphere(cell.cell_idx, arc, arc->next[0], beach_head);
                
                // This is not really a vertex. It is in the middle of some edge.
                PointCartesian vertex = phi_to_point(cur_site, cell.site.phi, sweep_line, sin_sweep_line, cos_sweep_line).get_cartesian();
                
                add_half_edge_sphere(voronoi_diagram, cells, half_edges, vertex, arc, arc->next[0]);
                add_half_edge_sphere(voronoi_diagram, cells, half_edges, vertex, arc->next[0], arc->next[0]->next[0]);
                
                //check for new circle events
                check_circle_event(arc, circle_event_queue_ptr, cells);
                check_circle_event(arc->next[0]->next[0], circle_event_queue_ptr, cells);
                
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

    void handle_circle_event(CircleEventSphere * event, VoronoiDiagramSphere * voronoi_diagram, vector<VoronoiCellSphere> * cells, vector<HalfEdgeSphere> * half_edges, priority_queue<CircleEventSphere *, vector<CircleEventSphere *>, PriorityQueueCompare> * circle_event_queue_ptr, ArcSphere * & beach_head)
    {
        //cout << "Handle circle event " << event->circumcenter;

        ArcSphere * left = event->arc->prev[0];
        ArcSphere * right = event->arc->next[0];
        
        //add new vertex
        PointCartesian vertex = event->circumcenter.get_cartesian();
        
        //add new edge
        add_half_edge_sphere(voronoi_diagram, cells, half_edges, vertex, left, right);
        
        //finish old edges
        int left_id = event->arc->left_edge_idx;
        int right_id = event->arc->right_edge_idx;
        if (left_id != -1)
        {
            finish_half_edge_sphere(voronoi_diagram, half_edges, left_id, vertex);
        }
        if (right_id != -1)
        {
            finish_half_edge_sphere(voronoi_diagram, half_edges, right_id, vertex);
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
        
        remove_arc_sphere(event->arc, beach_head);
        
        //check for new circle events
        check_circle_event(left, circle_event_queue_ptr, cells);
        check_circle_event(right, circle_event_queue_ptr, cells);
    }

    void check_circle_event(ArcSphere * arc, priority_queue<CircleEventSphere *, vector<CircleEventSphere *>, PriorityQueueCompare> * circle_event_queue_ptr, vector<VoronoiCellSphere> * cells)
    {
        if (arc == NULL || arc->prev[0] == NULL || arc->next[0] == NULL || arc->prev[0] == arc->next[0] || arc == arc->next[0] || arc->prev[0] == arc)
        {
            //cout << "Invalid circle event.\n";
            return;
        }
        
        PointSphere circumcenter;
        Real lowest_theta;
        
        PointSphere cur_site = (*cells)[arc->cell_idx].site;
        PointSphere prev_site = (*cells)[arc->prev[0]->cell_idx].site;
        PointSphere next_site = (*cells)[arc->next[0]->cell_idx].site;
        
        make_circle(prev_site, cur_site, next_site, circumcenter, lowest_theta);
        
        //This if statement causes errors for some reason...
        //if (lowest_theta > sweep_line)
        {
            arc->event = new CircleEventSphere(arc, circumcenter, lowest_theta);
            circle_event_queue_ptr->push(arc->event);
        }
    }

    void make_circle(PointSphere a, PointSphere b, PointSphere c, PointSphere & circumcenter, Real & lowest_theta)
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

    bool parabolic_intersection(PointSphere left, PointSphere right, Real & phi_intersection, ArcSphere * & beach_head, Real sweep_line, Real sin_sweep_line, Real cos_sweep_line)
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
        else if (left.theta == sweep_line /*&& right.theta < sweep_line*/)
        {
            /*
             *  The left site is on our sweep line so it contains our intersection phi.
             */
            phi_intersection = left.phi;
            //cout << "Left site is on sweep line.\n";
            return true;
        }
        else if (right.theta == sweep_line /*&& left.theta < sweep_line*/)
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
        
        if (abs(e) > sqrt_a_b)
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
            //cout << "Arc is close to sweep line. " << setprecision(50) << abs(e) << " > " << sqrt_a_b << endl << left << right << sweep_line << endl;
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

    PointSphere phi_to_point(PointSphere arc, Real phi, Real sweep_line, Real sin_sweep_line, Real cos_sweep_line)
    {
        Real a = cos(arc.theta) - cos_sweep_line;
        Real b = sin_sweep_line - sin(arc.theta) * cos(phi - arc.phi);
        return PointSphere(atan2(a, b), phi);
    }

    void add_initial_arc_sphere(int cell_id, ArcSphere * & beach_head)
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

    void add_arc_sphere(int cell_id, ArcSphere * left, ArcSphere * right, ArcSphere * & beach_head)
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

    void remove_arc_sphere(ArcSphere * arc, ArcSphere * & beach_head)
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

    void add_half_edge_sphere(VoronoiDiagramSphere * voronoi_diagram, vector<VoronoiCellSphere> * cells, vector<HalfEdgeSphere> * half_edges, PointCartesian start, ArcSphere * left, ArcSphere *right)
    {
        int edge_id = (int)half_edges->size();
        int voronoi_vertex_id = (int)voronoi_diagram->voronoi_verticies.size();
        
        voronoi_diagram->delaunay_edges.push_back(Edge(left->cell_idx, right->cell_idx));
        
        voronoi_diagram->voronoi_verticies.push_back(start);
        
        HalfEdgeSphere half_edge(voronoi_vertex_id);
        half_edges->push_back(half_edge);
        
        (*cells)[left->cell_idx].edge_ids.push_back(edge_id);
        (*cells)[right->cell_idx].edge_ids.push_back(edge_id);
        
        left->right_edge_idx = right->left_edge_idx = edge_id;
    }

    void finish_half_edge_sphere(VoronoiDiagramSphere * voronoi_diagram, vector<HalfEdgeSphere> * half_edges, int edge_idx, PointCartesian end)
    {
        if (!(*half_edges)[edge_idx].is_finished)
        {
            (*half_edges)[edge_idx].end_idx = (int)voronoi_diagram->voronoi_verticies.size();
            (*half_edges)[edge_idx].is_finished = true;
            
            voronoi_diagram->voronoi_verticies.push_back(end);
            
            voronoi_diagram->voronoi_edges.push_back(Edge((*half_edges)[edge_idx].start_idx, (*half_edges)[edge_idx].end_idx));
        }
    }

    ArcSphere * traverse_skiplist_to_site(ArcSphere * arc, Real phi, vector<VoronoiCellSphere> * cells)
    {
        int level = arc->height - 1;
        
        bool move_right = true;
        
        while (level >= 0)
        {
            Real delta = abs(phi - (*cells)[arc->cell_idx].site.phi);
            if (delta > M_PI) {delta = 2.f * M_PI - delta;}
            
            Real next_delta = abs(phi - (*cells)[arc->next[level]->cell_idx].site.phi);
            if (next_delta > M_PI) {delta = 2.f * M_PI - next_delta;}

            Real prev_delta = abs(phi - (*cells)[arc->prev[level]->cell_idx].site.phi);
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

    int random_height()
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
}
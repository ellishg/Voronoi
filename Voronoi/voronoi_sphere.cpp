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
    
    VoronoiDiagramSphere generate_voronoi(std::vector<std::tuple<Real, Real, Real>> * verts, THREAD_NUMBER num_threads, void (*render)(VoronoiDiagramSphere, ArcSphere *, vector<VoronoiCellSphere> *, Real), bool (*is_sleeping)())
    {
        switch (num_threads) {
            case ONE_THREAD:
                return generate_voronoi_one_thread(verts, render, is_sleeping);
                break;
            case TWO_THREADS:
                return generate_voronoi_two_threads(verts);
                break;
            case FOUR_THREADS:
                return generate_voronoi_four_threads(verts);
                break;
        }
    }

    VoronoiDiagramSphere generate_voronoi_four_threads(vector<tuple<Real, Real, Real>> * verts)
    {
        /*
         *  The voronoi diagram is computed from sites on the sphere corresponding
         *  to the point on a regular tetrahedron. The sweepline needs to continue
         *  until it has a radius of ARCTAN_2_ROOT_2. This value was obtained
         *  from the radius of the circumcircle from three points of the tetrahedron.
         *  A: (0, 0)
         *  B: (arcsin(1/3) + PI/2, 0)
         *  C: (arcsin(1/3) + PI/2, 2PI/3)
         *  D: (arcsin(1/3) + PI/2, 4PI/3)
         */
        
        VoronoiDiagramSphere diagram_a, diagram_b, diagram_c, diagram_d;
        
        vector<tuple<Real, Real, Real>> b_verts, c_verts, d_verts;
        
        for (auto point : *verts)
        {
            // We transform vertices to corresponding points on the tetrahedron
            auto b = rotate_y<Real>(point, SIN_NEG_ARCSIN_ONE_THIRD_PLUS_PI_2, COS_NEG_ARCSIN_ONE_THIRD_PLUS_PI_2);
            
            auto c = rotate_z<Real>(point, SIN_TWO_PI_3, COS_TWO_PI_3);
            c = rotate_y<Real>(c, SIN_NEG_ARCSIN_ONE_THIRD_PLUS_PI_2, COS_NEG_ARCSIN_ONE_THIRD_PLUS_PI_2);
            
            auto d = rotate_z<Real>(point, SIN_FOUR_PI_3, COS_FOUR_PI_3);
            d = rotate_y<Real>(d, SIN_NEG_ARCSIN_ONE_THIRD_PLUS_PI_2, COS_NEG_ARCSIN_ONE_THIRD_PLUS_PI_2);
            
            b_verts.push_back(b);
            c_verts.push_back(c);
            d_verts.push_back(d);
        }
        
        thread b_thread(compute_priority_queues, &diagram_b, &b_verts, ARCTAN_2_ROOT_2);
        thread c_thread(compute_priority_queues, &diagram_c, &c_verts, ARCTAN_2_ROOT_2);
        thread d_thread(compute_priority_queues, &diagram_d, &d_verts, ARCTAN_2_ROOT_2);
        
        compute_priority_queues(&diagram_a, verts, ARCTAN_2_ROOT_2);
        
        b_thread.join();
        c_thread.join();
        d_thread.join();
        
        unsigned long a_voronoi_vertex_length = diagram_a.voronoi_vertices.size();
        unsigned long b_voronoi_vertex_length = diagram_b.voronoi_vertices.size();
        unsigned long c_voronoi_vertex_length = diagram_c.voronoi_vertices.size();

        // Merge the voronoi vertices
        for (int i = 0; i < diagram_b.voronoi_vertices.size(); i++)
        {
            auto vert = diagram_b.voronoi_vertices[i];
            auto transformed_vert = rotate_y<Real>(make_tuple(vert.x, vert.y, vert.z), SIN_ARCSIN_ONE_THIRD_PLUS_PI_2, COS_ARCSIN_ONE_THIRD_PLUS_PI_2);
            
            diagram_a.voronoi_vertices.push_back(PointCartesian(get<0>(transformed_vert), get<1>(transformed_vert), get<2>(transformed_vert)));
        }
        for (int i = 0; i < diagram_c.voronoi_vertices.size(); i++)
        {
            auto vert = diagram_c.voronoi_vertices[i];
            auto transformed_vert = rotate_y<Real>(make_tuple(vert.x, vert.y, vert.z), SIN_ARCSIN_ONE_THIRD_PLUS_PI_2, COS_ARCSIN_ONE_THIRD_PLUS_PI_2);
            transformed_vert = rotate_z<Real>(transformed_vert, SIN_FOUR_PI_3, COS_FOUR_PI_3);
            
            diagram_a.voronoi_vertices.push_back(PointCartesian(get<0>(transformed_vert), get<1>(transformed_vert), get<2>(transformed_vert)));
        }
        for (int i = 0; i < diagram_d.voronoi_vertices.size(); i++)
        {
            auto vert = diagram_d.voronoi_vertices[i];
            auto transformed_vert = rotate_y<Real>(make_tuple(vert.x, vert.y, vert.z), SIN_ARCSIN_ONE_THIRD_PLUS_PI_2, COS_ARCSIN_ONE_THIRD_PLUS_PI_2);
            transformed_vert = rotate_z<Real>(transformed_vert, SIN_TWO_PI_3, COS_TWO_PI_3);
            
            diagram_a.voronoi_vertices.push_back(PointCartesian(get<0>(transformed_vert), get<1>(transformed_vert), get<2>(transformed_vert)));
        }

        // Merge the voronoi edges
        for (auto voronoi_edge : diagram_b.voronoi_edges)
        {
            voronoi_edge.vidx[0] += a_voronoi_vertex_length;
            voronoi_edge.vidx[1] += a_voronoi_vertex_length;

            diagram_a.voronoi_edges.push_back(voronoi_edge);
        }
        for (auto voronoi_edge : diagram_c.voronoi_edges)
        {
            voronoi_edge.vidx[0] += a_voronoi_vertex_length + b_voronoi_vertex_length;
            voronoi_edge.vidx[1] += a_voronoi_vertex_length + b_voronoi_vertex_length;
            
            diagram_a.voronoi_edges.push_back(voronoi_edge);
        }
        for (auto voronoi_edge : diagram_d.voronoi_edges)
        {
            voronoi_edge.vidx[0] += a_voronoi_vertex_length + b_voronoi_vertex_length + c_voronoi_vertex_length;
            voronoi_edge.vidx[1] += a_voronoi_vertex_length + b_voronoi_vertex_length + c_voronoi_vertex_length;
            
            diagram_a.voronoi_edges.push_back(voronoi_edge);
        }
        
        // Merge the delaunay edges
        for (auto delaunay_edge : diagram_b.delaunay_edges)
        {
            diagram_a.delaunay_edges.push_back(delaunay_edge);
        }
        for (auto delaunay_edge : diagram_c.delaunay_edges)
        {
            diagram_a.delaunay_edges.push_back(delaunay_edge);
        }
        for (auto delaunay_edge : diagram_d.delaunay_edges)
        {
            diagram_a.delaunay_edges.push_back(delaunay_edge);
        }

        return diagram_a;
    }

    VoronoiDiagramSphere generate_voronoi_two_threads(vector<tuple<Real, Real, Real>> * verts)
    {
        VoronoiDiagramSphere diagram_top_down, diagram_bottom_up;
        
        vector<tuple<Real, Real, Real>> bottom_up_verts;
        
        for (auto point : *verts)
        {
            bottom_up_verts.push_back(make_tuple(get<0>(point), get<1>(point), -get<2>(point)));
        }
     
        // Create a new thread to process the southern hemisphere
        thread diagram_bottom_up_thread(compute_priority_queues, &diagram_bottom_up, &bottom_up_verts, M_PI_2);
        
        compute_priority_queues(&diagram_top_down, verts, M_PI_2);

        diagram_bottom_up_thread.join();
        
        // Merge the two voronoi diagrams
        unsigned int vertices_length = (unsigned int)diagram_top_down.voronoi_vertices.size();
        
        for (auto vertex : diagram_bottom_up.voronoi_vertices)
        {
            diagram_top_down.voronoi_vertices.push_back(PointCartesian(vertex.x, vertex.y, -vertex.z));
        }
        
        for (auto voronoi_edge : diagram_bottom_up.voronoi_edges)
        {
            voronoi_edge.vidx[0] += vertices_length;
            voronoi_edge.vidx[1] += vertices_length;
            
            diagram_top_down.voronoi_edges.push_back(voronoi_edge);
        }
        
        for (auto delaunay_edge : diagram_bottom_up.delaunay_edges)
        {
            diagram_top_down.delaunay_edges.push_back(delaunay_edge);
        }
        
        return diagram_top_down;
    }
    
    void compute_priority_queues(VoronoiDiagramSphere * voronoi_diagram, vector<tuple<Real, Real, Real>> * verts, Real bound_theta)
    {        
        Real sweep_line = 0;
        
        ArcSphere * beach_head = NULL;
                
        vector<HalfEdgeSphere> half_edges;
        
        vector<VoronoiCellSphere> cells;

        priority_queue<VoronoiCellSphere, vector<VoronoiCellSphere>, PriorityQueueCompare> site_event_queue;

        priority_queue<CircleEventSphere *, vector<CircleEventSphere *>, PriorityQueueCompare> circle_event_queue;
        
        for (auto point : *verts)
        {
            PointCartesian point_cartesian = PointCartesian(get<0>(point), get<1>(point), get<2>(point));
   
            VoronoiCellSphere cell = VoronoiCellSphere(point_cartesian, (unsigned int)cells.size());
            
            site_event_queue.push(cell);
            
            cells.push_back(cell);
            
            voronoi_diagram->sites.push_back(point_cartesian);
        }
        
        while (!site_event_queue.empty() || !circle_event_queue.empty())
        {
            if (sweep_line > bound_theta)
            {
                // We need to check if this hemisphere is finished
                bool is_finished = true;
                ArcSphere * cur = beach_head;
                do
                {
                    if (cells[cur->cell_idx].site.theta < bound_theta)
                    {
                        // There is still a site in the beachline we need to process
                        is_finished = false;
                        break;
                    }
                } while ((cur = cur->next[0]) != beach_head);

                if (is_finished) {break;} // We have finished this hemisphere
            }
            
            if (!circle_event_queue.empty() && !circle_event_queue.top()->is_valid)
            {
                // Remove invalid circle events
                delete circle_event_queue.top();
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
        
        // Clean up
        while (beach_head->next[0] != beach_head && beach_head->prev[0] != beach_head)
        {
            // Deallocate the beachline
            remove_arc_sphere(beach_head, beach_head);
        }
        remove_arc_sphere(beach_head, beach_head);
        
        while (!circle_event_queue.empty())
        {
            // Deallocate all the remaining circle events
            delete circle_event_queue.top();
            circle_event_queue.pop();
        }
    }

    VoronoiDiagramSphere generate_voronoi_one_thread(vector<tuple<Real, Real, Real>> * verts, void (*render)(VoronoiDiagramSphere, ArcSphere *, vector<VoronoiCellSphere> *, Real), bool (*is_sleeping)())
    {
        bool should_render = (false) && render != NULL && is_sleeping != NULL;

        VoronoiDiagramSphere voronoi_diagram;
        
        ArcSphere * beach_head = NULL;
        
        Real sweep_line = 0;
        
        vector<HalfEdgeSphere> half_edges;
        
        vector<VoronoiCellSphere> cells;

        priority_queue<VoronoiCellSphere, vector<VoronoiCellSphere>, PriorityQueueCompare> site_event_queue;
        
        priority_queue<CircleEventSphere *, vector<CircleEventSphere *>, PriorityQueueCompare> circle_event_queue;
        
        for (auto point : *verts)
        {
            PointCartesian point_cartesian = PointCartesian(get<0>(point), get<1>(point), get<2>(point));
            
            VoronoiCellSphere cell = VoronoiCellSphere(point_cartesian, (unsigned int)cells.size());
            
            site_event_queue.push(cell);
            
            cells.push_back(cell);
            
            voronoi_diagram.sites.push_back(point_cartesian);
        }
        
        while (!site_event_queue.empty() || !circle_event_queue.empty())
        {
            //if (sweep_line >= 1.00736) {should_render = true;}
            //cout << sweep_line << endl;
            
            //cout << "[" << site_event_queue.size() << ", " << circle_event_queue.size() << "]\n";
            
            if (!circle_event_queue.empty() && !circle_event_queue.top()->is_valid)
            {
                // Remove invalid circle events
                delete circle_event_queue.top();
                circle_event_queue.pop();
            }
            else if (site_event_queue.empty() || (!circle_event_queue.empty() && circle_event_queue.top()->lowest_theta < site_event_queue.top().site.theta))
            {
                CircleEventSphere * circle = circle_event_queue.top();
                sweep_line = circle->lowest_theta;
                handle_circle_event(circle, &voronoi_diagram, &cells, &half_edges, &circle_event_queue, beach_head);
                delete circle;
                circle_event_queue.pop();
                
                while (should_render && is_sleeping())
                {
                    render(voronoi_diagram, beach_head, &cells, (Real)sweep_line);
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
                    render(voronoi_diagram, beach_head, &cells, (Real)sweep_line);
                }
            }
        }
        
        // Clean up
        while (beach_head != NULL && beach_head->next[0] != beach_head && beach_head->prev[0] != beach_head)
        {
            // Deallocate the beachline
            remove_arc_sphere(beach_head, beach_head);
        }
        remove_arc_sphere(beach_head, beach_head);
        
        // The final two edges need to be connected.
        for (int i = 0; i < half_edges.size(); i++)
        {
            if (!half_edges[i].is_finished)
            {
                for (int k = i + 1; k < half_edges.size(); k++)
                {
                    if (!half_edges[k].is_finished)
                    {
                        finish_half_edge_sphere(&voronoi_diagram, &half_edges, i, voronoi_diagram.voronoi_vertices[half_edges[k].start_idx]);
                        finish_half_edge_sphere(&voronoi_diagram, &half_edges, k, voronoi_diagram.voronoi_vertices[half_edges[i].start_idx]);
                        
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
        //cout << "Handle site event " << cell.site;
        
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

            Real phi_start = 0;
            Real phi_end = 0;
            
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
                assert(0);
                return;
            }
            else if (valid_arc && ((phi_start < phi_end && phi_start <= cell.site.phi && cell.site.phi <= phi_end) || (phi_start > phi_end && (phi_start <= cell.site.phi || cell.site.phi <= phi_end))))
            {
                // The arc is found!
                /*cout << setprecision(16) << hexfloat;
                cout << "Prev arc = " << (*cells)[arc->prev[0]->cell_idx].site;
                cout << "Cur arc = " << (*cells)[arc->cell_idx].site;
                cout << "Next arc = " << (*cells)[arc->next[0]->cell_idx].site;
                cout << "Sweepline = " << sweep_line << endl;
                cout << "Phi Start = " << phi_start << endl;
                cout << "Phi End = " << phi_end << endl << endl;
                */
                
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
        
        //This if statement breaks my code for some reason...
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
            cout << "Intersection is at north pole.\n";
            return false;
        }
        else if (left.theta == sweep_line)
        {
            /*
             *  The left site is on our sweep line so it contains our intersection phi.
             */
            phi_intersection = left.phi;
            //cout << "Left site is on sweep line.\n";
            return true;
        }
        else if (right.theta == sweep_line)
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
             *  Theoretically this is impossible because as the sweep line
             *  approaches either the left or right theta e / sqrt_a_b approaches
             *  1. In reality, e can be greater than sqrt_a_b due to rounding.
             *  This excecutes when a site is close enought to the sweep line
             *  to cause problems. So one of the parabolas is actually a line.
             *  The lower site contains the phi of our intersection.
             */
            phi_intersection = (left.theta > right.theta) ? left.phi : right.phi;
            /*cout << "cos(left.theta) = " << cos_left_theta << endl;
            cout << "a = " << a << "\nb = " << b << endl;
            cout << "Arc is close to sweep line.\n" << abs(e) << " > " << sqrt_a_b << endl << left << right << sweep_line << endl << endl;*/
            return true;
            //return false;
        }
        
        Real sin_phi_int_plus_gamma = e / sqrt_a_b;
        
        Real gamma = atan2(a, b);
        
        phi_intersection = asin(sin_phi_int_plus_gamma) - gamma;
        
        if (phi_intersection > M_PI)
        {
            phi_intersection -= 2 * M_PI;
        }
        else if (phi_intersection <= -M_PI)
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
        if (beach_head == NULL) {return;}
        
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
        int voronoi_vertex_id = (int)voronoi_diagram->voronoi_vertices.size();
        
        voronoi_diagram->delaunay_edges.push_back(Edge(left->cell_idx, right->cell_idx));
        
        voronoi_diagram->voronoi_vertices.push_back(start);
        
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
            (*half_edges)[edge_idx].end_idx = (int)voronoi_diagram->voronoi_vertices.size();
            (*half_edges)[edge_idx].is_finished = true;
            
            voronoi_diagram->voronoi_vertices.push_back(end);
            
            voronoi_diagram->voronoi_edges.push_back(Edge((*half_edges)[edge_idx].start_idx, (*half_edges)[edge_idx].end_idx));
        }
    }

    ArcSphere * traverse_skiplist_to_site(ArcSphere * arc, Real phi, vector<VoronoiCellSphere> * cells)
    {
        /*
         *  Remember that arcs on the beachline are duplicated when an arc is inserted.
         *  So that the beachline is not necessarily sorted by phi.
         *  This function does NOT work as you might expect because of this!
         *  But it will most likely move closer to the destination.
         */
        
        int level = arc->height - 1;
        
        while (level >= 0)
        {
            int prev_index = arc->prev[level]->cell_idx;
            int cur_index = arc->cell_idx;
            int next_index = arc->next[level]->cell_idx;
        
            Real delta = abs(phi - (*cells)[cur_index].site.phi);
            if (delta > M_PI) {delta = 2.f * M_PI - delta;}
            
            Real next_delta = abs(phi - (*cells)[next_index].site.phi);
            if (next_delta > M_PI) {next_delta = 2.f * M_PI - next_delta;}

            Real prev_delta = abs(phi - (*cells)[prev_index].site.phi);
            if (prev_delta > M_PI) {prev_delta = 2.f * M_PI - prev_delta;}
            
            if (next_delta < delta && next_delta < prev_delta) // if next_delta is the min
            {
                arc = arc->next[level];
            }
            else if (prev_delta < delta) // if prev_delta is the min
            {
                arc = arc->prev[level];
            }
            else // if delta is the min
            {
                level--;
            }
        }
        return arc;
    }
    
    template <typename T>
    tuple<T, T, T> rotate_y(tuple<T, T, T> point, T sin_theta, T cos_theta)
    {
        T x = get<0>(point);
        T y = get<1>(point);
        T z = get<2>(point);
        
        return make_tuple(cos_theta * x - sin_theta * z, y, sin_theta * x + cos_theta * z);
    }
    
    template <typename T>
    tuple<T, T, T> rotate_z(tuple<T, T, T> point, T sin_theta, T cos_theta)
    {
        T x = get<0>(point);
        T y = get<1>(point);
        T z = get<2>(point);
        
        return make_tuple(cos_theta * x + sin_theta * y, cos_theta * y - sin_theta * x, z);
    }
    
    int random_height()
    {
        int height = MAX_SKIPLIST_HEIGHT;
        int r = (rand() % (0b1 << MAX_SKIPLIST_HEIGHT)) + 1;
        while (--height > 0 && r < (0b1 << height));
        return MAX_SKIPLIST_HEIGHT - height;
    }
}
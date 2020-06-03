#include "Dijkstra_Path.h"


Dijkstra_Path::~Dijkstra_Path()
{
}

void Dijkstra_Path::ComputePath(VH & start, VH & end)
{
	std::priority_queue<Node> Q;

	mesh_.add_property(forward_he);
	mesh_.add_property(version);
	mesh_.add_property(dist);

	for (VH v_h:mesh_.vertices())
	{
		mesh_.property(dist, v_h) = DBL_MAX;
	}

	// 1
	for (HEH voh : mesh_.voh_range(start))
	{
		VH v_h = mesh_.to_vertex_handle(voh);

		// property
		mesh_.property(forward_he, v_h) = voh;
		mesh_.property(version, v_h) += 1;
		if (mesh_.is_boundary(mesh_.edge_handle(voh))) {
			mesh_.property(dist, v_h) = DBL_EPSILON;
		}
		else {
			mesh_.property(dist, v_h) = mesh_.calc_edge_length(voh);
		}

		// Node 
		Node temp_node;
		temp_node.v_h = v_h;
		temp_node.version = mesh_.property(version, v_h);
		temp_node.dist_ = mesh_.property(dist, v_h);
		Q.push(temp_node);
	}

	// 2
	while (!Q.empty()) {
		// 2.1
		Node temp_P = Q.top();
		Q.pop();

		// 2.2
		VH p = temp_P.v_h;
		if (mesh_.property(version, p) != temp_P.version) continue;

		// 2.3
		if (p == end) break;

		// 2.4
		double temp_dist = mesh_.property(dist, p);
		VH vv_h;
		double new_dist;
		for (HEH voh : mesh_.voh_range(p)) {
			vv_h = mesh_.to_vertex_handle(voh);
			new_dist = temp_dist + mesh_.calc_edge_length(voh);
			if (mesh_.property(dist, vv_h) > new_dist)
			{
				// property
				mesh_.property(forward_he, vv_h) = voh;
				mesh_.property(version, vv_h) += 1;
				mesh_.property(dist, vv_h) = new_dist;
				
				Node temp_node;
				temp_node.v_h = vv_h;
				temp_node.version = mesh_.property(version, vv_h);
				temp_node.dist_ = mesh_.property(dist, vv_h);
				Q.push(temp_node);

			}

	
		}
	}

	Path.clear();
	VH v_h = end;
	HEH he_h;
	do
	{
		he_h = mesh_.property(forward_he, v_h);
		Path.insert(Path.begin(), he_h);
		v_h = mesh_.from_vertex_handle(he_h);
	} while (v_h != start);

	mesh_.remove_property(forward_he);
	mesh_.remove_property(version);
	mesh_.remove_property(dist);
}

std::vector<HEH> Dijkstra_Path::return_path()
{
	return Path;
}

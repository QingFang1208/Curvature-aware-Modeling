#include "Ope.h"


Ope::~Ope()
{
}

void Ope::Split()
{
	for (auto e : mesh_.edges())
	{
		if (!mesh_.property(is_remesh_e_, e)) continue;

		if (mesh_.calc_edge_length(e) > high_)
		{
			auto mid_point = mesh_.calc_edge_midpoint(e);
			auto v_h = mesh_.split(e, mid_point);

			mesh_.property(is_remesh_v_, v_h) = true;
			for (auto e_h:mesh_.ve_range(v_h))
			{
				mesh_.property(is_remesh_e_, e_h) = true;
			}
		}
	}
}

void Ope::Collapse()
{
	// add checked property
	OpenMesh::EPropHandleT<bool> checked;
	mesh_.add_property(checked);

	// init property
	for (auto e:mesh_.edges()) {
		if (!mesh_.property(is_remesh_e_, e)) mesh_.property(checked, e) = true;
		else mesh_.property(checked, e) = false;
	}

	bool finished = false;
	while (!finished) {

		finished = true;

		for (auto e_h:mesh_.edges()) {

			if (mesh_.property(checked, e_h)) continue;
			
			mesh_.property(checked, e_h) = true;
			auto he_h = mesh_.halfedge_handle(e_h, 0);

			auto from_v = mesh_.from_vertex_handle(he_h);
			auto to_v = mesh_.to_vertex_handle(he_h);

			double edgeLength = mesh_.calc_edge_length(e_h);
			if ((edgeLength < low_) && (edgeLength > DBL_EPSILON)) {
				auto to_point = mesh_.point(to_v);
				bool collapse_ok = true;

				for (auto he_h : mesh_.voh_range(from_v)) {
					if (((to_point - mesh_.point(mesh_.to_vertex_handle(he_h))).norm() > high_)
						|| (mesh_.status(mesh_.edge_handle(he_h)).feature())
						|| (mesh_.is_boundary(mesh_.edge_handle(he_h)))) {
						collapse_ok = false;
						break;
					}
				}
					
				if (collapse_ok && mesh_.is_collapse_ok(he_h)) {
					mesh_.collapse(he_h);

					finished = false;
				}
			}
		}
	}
	mesh_.remove_property(checked);
	// mesh_.garbage_collection();
}

void Ope::Equalize_valences()
{
	// valences
	OpenMesh::VPropHandleT<int> reg_val;
	mesh_.add_property(reg_val);
	for (auto v:mesh_.vertices())
	{
		if (mesh_.is_boundary(v)) mesh_.property(reg_val, v) = 4;
		else mesh_.property(reg_val, v) = 6;
	}

	for (auto e_h:mesh_.edges()) {

		if (!mesh_.property(is_remesh_e_, e_h)) continue;
		if (!mesh_.is_flip_ok(e_h)) continue;
		if (mesh_.status(e_h).feature()) continue;

		const auto& h0 = mesh_.halfedge_handle(e_h, 0);
		const auto& h1 = mesh_.halfedge_handle(e_h, 1);

		if (h0.is_valid() && h1.is_valid())
			if (mesh_.face_handle(h0).is_valid() && mesh_.face_handle(h1).is_valid()) {
				//get vertices of corresponding faces
				const auto& a = mesh_.to_vertex_handle(h0);
				const auto& b = mesh_.to_vertex_handle(h1);
				const auto& c = mesh_.to_vertex_handle(mesh_.next_halfedge_handle(h0));
				const auto& d = mesh_.to_vertex_handle(mesh_.next_halfedge_handle(h1));

				// not flip for boundary point
				if (is_bound_flip(a, b, c, d)) continue;

				const int deviation_pre = abs((int)(mesh_.valence(a) - mesh_.property(reg_val, a)))
					+ abs((int)(mesh_.valence(b) - mesh_.property(reg_val, b)))
					+ abs((int)(mesh_.valence(c) - mesh_.property(reg_val, c)))
					+ abs((int)(mesh_.valence(d) - mesh_.property(reg_val, d)));
				mesh_.flip(e_h);

				const int deviation_post = abs((int)(mesh_.valence(a) - mesh_.property(reg_val, a)))
					+ abs((int)(mesh_.valence(b) - mesh_.property(reg_val, b)))
					+ abs((int)(mesh_.valence(c) - mesh_.property(reg_val, c)))
					+ abs((int)(mesh_.valence(d) - mesh_.property(reg_val, d)));

				if (deviation_pre <= deviation_post)
					mesh_.flip(e_h);
			}
	}
	mesh_.remove_property(reg_val);
}

void Ope::Tangential_relaxation()
{
	mesh_.update_normals();

	//add center property
	OpenMesh::VPropHandleT<Mesh::Point> barycenter;
	mesh_.add_property(barycenter);

	//first compute barycenters
	for (auto v_h:mesh_.vertices()) {

		Mesh::Point tmp(0.0, 0.0, 0.0);
		for (auto vv_h:mesh_.vv_range(v_h)) {
			tmp += mesh_.point(vv_h);
		}

		mesh_.property(barycenter, v_h) = tmp / mesh_.valence(v_h);
	}

	//move to new position
	for (auto v_h:mesh_.vertices()) {
		if (!mesh_.is_boundary(v_h) && mesh_.property(is_remesh_v_, v_h))
			mesh_.set_point(v_h, mesh_.property(barycenter, v_h) + 
			(mesh_.normal(v_h) | (mesh_.point(v_h) - mesh_.property(barycenter, v_h))) * mesh_.normal(v_h));
	}

	mesh_.remove_property(barycenter);

}

void Ope::compute_average_edge_length()
{

	double avelen = 0;
	int cont = 0;
	for (auto e : mesh_.edges())
	{
		if (mesh_.property(is_remesh_e_, e)) continue;

		avelen += mesh_.calc_edge_length(e);
		cont++;
	}

	avelen /= cont;

	high_ = 4.0 / 3 * avelen;
	low_ = 4.0 / 5 * avelen;
}

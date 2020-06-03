#include "Mesh_Ope.h"
#include "Ope.h"

Mesh_Ope::~Mesh_Ope()
{
}

Mesh Mesh_Ope::Cut_Seg(int i, std::unordered_map<VH, VH>& old_to_new, std::unordered_map<VH, VH>& new_to_old)
{
	Mesh seg_mesh = mesh_;

	OpenMesh::VPropHandleT<OpenMesh::VertexHandle> old_handle;
	seg_mesh.add_property(old_handle);
	for (VH v : seg_mesh.vertices())
	{
		seg_mesh.property(old_handle, v) = v;
	}


	// delete
	for (FH f_h : seg_mesh.faces())
	{
		if (seg_mesh.property(seg_, f_h) != i) {
			seg_mesh.delete_face(f_h);
		}
	}
	seg_mesh.garbage_collection();

	// construct relation between old and new
	new_to_old.clear();
	old_to_new.clear();
	for (VH v : seg_mesh.vertices())
	{
		old_to_new.insert(std::make_pair(seg_mesh.property(old_handle, v), v));
		new_to_old.insert(std::make_pair(v, seg_mesh.property(old_handle, v)));
	}


	seg_mesh.remove_property(old_handle);

	return seg_mesh;
}

FSet Mesh_Ope::Close_draw_bound(HESet& bound, int & seg_num)
{
	FSet change_face(0);

	FH f_h = mesh_.face_handle(bound[0]);
	int cur_seg = mesh_.property(seg_, f_h);


	FSet stackV(0);
	for (HEH he_h : bound)
	{
		// property
		FH f_h = mesh_.face_handle(he_h);
		mesh_.property(seg_, f_h) = seg_num;
		change_face.push_back(f_h);

		// init stackV
		if (stackV.size() != 0) continue;

		HEH	next_oppo_next = mesh_.next_halfedge_handle(mesh_.opposite_halfedge_handle(mesh_.next_halfedge_handle(he_h)));
		if (find(bound.begin(), bound.end(), next_oppo_next) == bound.end())
		{
			FH init_f = mesh_.face_handle(next_oppo_next);
			stackV.push_back(init_f);
		}
	}

	while (stackV.size() != 0)
	{
		FH cur_f_h = stackV.back();
		stackV.pop_back();

		for (FH f_h : mesh_.ff_range(cur_f_h))
		{
			if (mesh_.property(seg_, f_h) != cur_seg) continue;

			// property
			mesh_.property(seg_, f_h) = seg_num;
			change_face.push_back(f_h);

			stackV.push_back(f_h);
		}
	}

	seg_num++;
	return change_face;
}



#ifndef INTERACTIVE_VIEWER_WIDGET
#define INTERACTIVE_VIEWER_WIDGET

#include "MeshViewerWidget.h"
#include "ANN\ANN.h"
#include "Dijkstra_Path.h"
#include<QKeyEvent>

extern const double max_dists;

class InteractiveViewerWidget : public MeshViewerWidget
{
	Q_OBJECT
public:
	InteractiveViewerWidget(QWidget* parent = 0);
	InteractiveViewerWidget(QGLFormat& _fmt, QWidget* _parent);
	~InteractiveViewerWidget();

	virtual void clearMesh();

	// Draw Seg
	void Draw_Undo();
	void Finish_Seg_Draw();
	void Undo_Finished();
	void Undo_Seg_Change();
	
public slots:
	void Fix_Vertices() {
		fixedVertices.clear();
		fixedPositions.clear();
		switch (edit_mode_)
		{
		case BOUND_PICK:
		{
			for (size_t i = 0; i < segment.n_contours(); i++)
			{
				if (!contour_pick_[i]) continue;
				auto c = segment.contour(i);
				for (auto he : c)
				{
					fixedVertices.push_back(segment.curMesh().to_vertex_handle(he));
					fixedPositions.push_back(segment.curMesh().point(segment.curMesh().to_vertex_handle(he)));
				}
			}
			break;
		}
		case VERT_SELECT: {
			for (auto v_h : selectedVertices)
			{
				fixedVertices.push_back(v_h);
				fixedPositions.push_back(segment.curMesh().point(v_h));
			}
			selectedVertices.clear();
			break;
		}
		default:
			break;
		}

		// arap
		/*arap_set_handle(fixedVertices);*/

	}
	
signals:
	void mouse_press_signal(Mesh::Point P);
	void mouse_move_signal(OpenMesh::Vec3d xy);
	void mouse_release_signal(Mesh::Point P);
	void draw_from_out_signal();
	void setEditMode_signal(int);

public:
	enum { TRANS, SEG_DRAW, BOUND_PICK, SEG_PICK, SEG_CHANGE,  VERT_SELECT, FACE_SELECT, VERT_DRAG, INIT_SPHERE, INIT_CYLINDER};
	void setEditMode(int mm);
	//int editMode() const { return edit_mode_; }

protected:
	virtual void mousePressEvent(QMouseEvent *_event);
	virtual void mouseReleaseEvent(QMouseEvent *_event);
	virtual void mouseMoveEvent(QMouseEvent *_event);
	virtual void wheelEvent(QWheelEvent* _event);
	void keyPressEvent(QKeyEvent *event)
	{
		switch (event->key())
		{
		case Qt::Key_Shift:
			addPart = true;
			break;
		default:
			break;
		}
	}
	void keyReleaseEvent(QKeyEvent *event)
	{
		switch (event->key())
		{
		case Qt::Key_Shift:
			addPart = false;
			break;
		default:
			break;
		}
	}

	int edit_mode_;

protected:
	void pick_vertex(int x,int y);
	void pick_face(int x,int y);
	void pick_edge(int x,int y);
	void pick_point(int x,int y);
	void move_point_based_lastVertex(int x,int y);

	int find_vertex_using_selected_point();
	int find_face_using_selected_point();
	int find_edge_using_selected_point();

	void buildIndex();
	ANNkd_tree* kdTree;

	void draw_picked_seg();
	void draw_picked_bound();
	void draw_selected_faces();
	void draw_selected_vertices();
	void draw_seg_lines();
	void draw_fixed_vertices();
	void draw_init_sphere();
	void draw_init_cylinder();
	virtual void draw_scene(int drawmode);


protected:
	double selectedPoint[3];
	int lastestVertex;
	int lastestFace;
	int lastestEdge;
	bool addPart = false;

	// Path Draw
	Dijkstra_Path* dijkstra_;
	std::vector<HEH> cur_path_;
	std::vector<std::vector<HEH>> path_;
	VH init_point;
	VH start_point;
	VH end_point;

	// Selected Draw
	QPoint select_start_, select_end_;
	std::vector<OpenMesh::FaceHandle> selectedFaces;
	void compute_selected_faces();

	// Add Seg
	std::vector<OpenMesh::FaceHandle> last_faces_vector;
	std::vector<int> last_faces_seg_vector;

	// Change Seg
	std::vector<OpenMesh::FaceHandle> change_face_vector;
	std::vector<int> change_face_seg_vector;

	// Edit
	std::vector<OpenMesh::VertexHandle> selectedVertices;

	void compute_selected_vertices();


	

protected:
	void dragEnterEvent(QDragEnterEvent *event);
	void dropEvent(QDropEvent *event);
	//------------------------------------------------------------------
	// ************  add for current mesh selected point  **************
	ANNkd_tree *curMesh_kdTree = NULL;
	void build_curMesh_tree();
	int pick_curMesh_idx();


	//------------------------------------------------------------------
	// edge pick
	ANNkd_tree* contour_kdTree = NULL;
	std::unordered_map<int, int> treeidx_to_bound;
	void build_bound_tree();
	int pick_contour();

	//ARAP* arap;
	//Mesh temp_mesh;
	//std::unordered_map<VH, VH> old_to_new;
	//std::unordered_map<VH, VH> new_to_old;
	//void init_arap() {
	//	temp_mesh = ope->Cut_Seg(pick_seg_, old_to_new, new_to_old);
	//	arap = new ARAP(temp_mesh);
	//}
	//void arap_set_handle(std::vector<OpenMesh::VertexHandle> fix_handle) {
	//	std::vector<OpenMesh::VertexHandle> fix_point_idx;
	//	for (auto v_h : fix_handle)
	//	{
	//		fix_point_idx.push_back(old_to_new[v_h]);
	//	}
	//	arap->Set_Fix_Idx(fix_point_idx);
	//}

	//void arap_move_handle(OpenMesh::VertexHandle move_handle) {
	//	std::vector<OpenMesh::VertexHandle> move_point_idx = { old_to_new[move_handle] };
	//	arap->Set_Move_Idx(move_point_idx);
	//}

	//void arap_move_position(OpenMesh::VertexHandle move_handle, OpenMesh::Vec3d v_position) {
	//	std::vector<OpenMesh::VertexHandle> move_point_idx = { old_to_new[move_handle] };
	//	std::vector<OpenMesh::Vec3d> move_point_position = { v_position };
	//	arap->Set_Move_Position(move_point_idx, move_point_position);
	//	arap->Compute_ARAP();

	//	for (auto v_h : temp_mesh.vertices())
	//	{
	//		mesh.set_point(new_to_old[v_h], temp_mesh.point(v_h));
	//	}
	//}


private:
	inline bool is_in_rect(GLdouble x, GLdouble y) {
		return (x - select_start_.x()) * (x - select_end_.x()) < 0 &&
			(height() - y - select_start_.y()) * (height() - y - select_end_.y()) < 0;
	}

	bool pickFaceValid() {
		if (lastestFace >= 0 && lastestFace < mesh.n_faces())
			return true;
		return false;
	}

	bool pickEdgeValid() {
		if (lastestEdge >= 0 && lastestEdge < mesh.n_edges())
			return true;
		return false;
	}

	bool pickVertexValid() {
		if (lastestVertex >= 0 && lastestVertex < mesh.n_vertices())
			return true;
		return false;
	}
};

#endif
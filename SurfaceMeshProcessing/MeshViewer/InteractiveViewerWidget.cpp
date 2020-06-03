#include <QMouseEvent>
#include <QLineEdit>
#include <QDragEnterEvent>
#include <QDropEvent>
#include <QtCore>
#include <QUrl>
#include <ctime>
#include "InteractiveViewerWidget.h"

const double max_dists = 20;

InteractiveViewerWidget::InteractiveViewerWidget(QWidget* parent /* = 0 */)
	:MeshViewerWidget(parent)
{
	kdTree = NULL;
	//arap = NULL;
	contour_kdTree = NULL;

	// Draw
	dijkstra_ = new Dijkstra_Path(mesh);
	cur_path_.resize(0);
	path_.resize(0);

}

InteractiveViewerWidget::InteractiveViewerWidget(QGLFormat& _fmt, QWidget* _parent)
:MeshViewerWidget(_fmt, _parent)
{
	kdTree = NULL;
	//arap = NULL;
	contour_kdTree = NULL;

	// Draw
	dijkstra_ = new Dijkstra_Path(mesh);
	cur_path_.resize(0);
	path_.resize(0);
}

InteractiveViewerWidget::~InteractiveViewerWidget()
{
	/*if (arap) delete arap;*/
	if (kdTree) delete kdTree;
	if (dijkstra_) delete dijkstra_;
	if (contour_kdTree) delete contour_kdTree;
}

void InteractiveViewerWidget::clearMesh()
{
	MeshViewerWidget::clearMesh();
}

void InteractiveViewerWidget::Draw_Undo()
{
	if (path_.size() == 0) return;

	path_.pop_back();
	cur_path_.clear();
	if (path_.size() > 1)
	{
		start_point = mesh.to_vertex_handle(path_.back().back());
		end_point = start_point;
	}
	if (path_.size() == 1)
	{
		start_point = init_point;
		end_point = init_point;
	}
	updateGL();
}

void InteractiveViewerWidget::Finish_Seg_Draw()
{
	last_faces_vector.clear();
	last_faces_seg_vector.clear();
	switch (edit_mode_)
	{
	case SEG_DRAW: {
		if (init_point == end_point && path_.size() != 0)
		{
			std::vector<HEH> bound(0);
			for (size_t i = 0; i < path_.size(); i++) {
				for (HEH he_h : path_[i]) {
					bound.push_back(he_h);
				}
			}

			// last
			int ori_seg_idx = mesh.property(seg, mesh.face_handle(bound[0]));
			last_faces_vector = ope->Close_draw_bound(bound, seg_num);
			last_faces_seg_vector.resize(last_faces_vector.size(), ori_seg_idx);

			updateGL();

			path_.clear();
			cur_path_.clear();
		}
		else {
			std::cerr << "It is not a close loop!" << std::endl;
		}
		break;
	}
	case FACE_SELECT: {
		for (auto f_h : selectedFaces)
		{
			// last
			last_faces_vector.push_back(f_h);
			last_faces_seg_vector.push_back(mesh.property(seg, f_h));

			mesh.property(seg, f_h) = seg_num;
		}

		seg_num++;
		break;
	}
	default:
		break;
	}
	std::cout << "Seg Num:" << seg_num << std::endl;
	updateGL();
}

void InteractiveViewerWidget::Undo_Finished()
{
	if (seg_num > 1)
	{
		int cont = 0;
		for (auto f_h : last_faces_vector)
		{
			mesh.property(seg, f_h) = last_faces_seg_vector[cont];
			cont++;
		}
		seg_num--;
		last_faces_vector.clear();
		last_faces_seg_vector.clear();
		std::cout << "Seg Num:" << seg_num << std::endl;
		updateGL();
	}
}

void InteractiveViewerWidget::Undo_Seg_Change()
{
	if (change_face_vector.size() > 0)
	{
		auto f_h = change_face_vector.back();
		change_face_vector.pop_back();
		int temp_seg = change_face_seg_vector.back();
		change_face_seg_vector.pop_back();
		mesh.property(seg, f_h) = temp_seg;
		std::cout << "Undo" << std::endl;
		updateGL();
	}
}

void InteractiveViewerWidget::setEditMode(int mm)
{
	edit_mode_ = mm;
	if (TRANS != edit_mode_)
	{
		buildIndex();
	}
	if (FACE_SELECT == edit_mode_)
	{
		selectedFaces.clear();
	}
	//pick_bound_ = -1; 
	//pick_seg_ = -1;
	emit setEditMode_signal(mm);
}

void InteractiveViewerWidget::mousePressEvent(QMouseEvent *_event)
{

	MeshViewerWidget::mousePressEvent(_event);
	
	pick_point(_event->x(), _event->y());

	if (edit_mode_ == SEG_DRAW)
	{
		pick_vertex(_event->x(), _event->y());
	}
	else if (edit_mode_ == FACE_SELECT)
	{
		select_start_ = _event->pos();
	}
	else if (edit_mode_ == SEG_CHANGE)
	{
		pick_face(_event->x(), _event->y());
		if (pickFaceValid())
		{
			/*auto f_h = mesh.face_handle(lastestFace);
			if (std::find(pick_segs.begin(), pick_segs.end(), mesh.property(seg, f_h)) != pick_segs.end())
			{
				change_face_vector.push_back(f_h);
				change_face_seg_vector.push_back(mesh.property(seg, f_h));
				mesh.property(seg, f_h) = pick_segs[0];

			}*/
		}
	}
	else if (edit_mode_ == SEG_PICK)
	{
		pick_face(_event->x(), _event->y());
		if (pickFaceValid() && seg.is_valid())
		{
			set_cur_seg(mesh.property(seg, mesh.face_handle(lastestFace)), addPart);
			build_curMesh_tree();
			build_bound_tree();

			//init_arap();

			/*segment.select(std::vector<int>(1, pick_seg_));
			contour_pick_ = std::vector<bool>(segment.n_contours(), false);*/
		}
	}
	else if (edit_mode_ == BOUND_PICK)
	{
		if (treeidx_to_bound.size() != 0) {
			int v_idx = pick_contour();
			if (v_idx != -1)
			{
				pick_bound_ = treeidx_to_bound[v_idx];
				contour_pick_[pick_bound_] = !contour_pick_[pick_bound_];

				std::cout << "Pick Bound: " << pick_bound_ << std::endl;
			}
		}
	}
	else if (edit_mode_ == VERT_SELECT)
	{
		select_start_ = _event->pos();
	}
	else if (edit_mode_ == VERT_DRAG)
	{
		pick_point(_event->x(), _event->y());
		ANNpoint tp = annAllocPt(3); tp[0] = selectedPoint[0]; tp[1] = selectedPoint[1]; tp[2] = selectedPoint[2];
		int neighNum = 1;
		ANNidxArray nnIdx = new ANNidx[neighNum]; ANNdistArray dists = new ANNdist[neighNum];
		curMesh_kdTree->annkSearch(tp, neighNum, nnIdx, dists);

		for (int i = 0; i < neighNum; i++)
		{
			if (dists[i] > max_dists)
			{
				nnIdx[i] = -1;
			}
		}

		lastestVertex = nnIdx[0];

		moveVertices.clear();
		for (int i = 0; i < neighNum; i++)
		{
			if (nnIdx[i] < 0) continue;
			moveVertices.push_back(segment.curMesh().vertex_handle(nnIdx[i]));
			/*for (auto vv : segment.curMesh().vv_range(moveVertices[0]))
			{
				moveVertices.push_back(vv);
			}*/
		}


		if (pickVertexValid())
		{
			/*auto v_h = mesh.vertex_handle(lastestVertex);
			auto temp_vector = fixedVertices;
			temp_vector.push_back(v_h);*/
			/*auto vInner = segmesh->inner();
			
			OpenMesh::Vec3d dis(1.0, 0.0, 0.0);

			deform->domain.build(vInner, temp_vector);
			std::cout << "Ori Position:" << mesh.point(v_h).data()[0]
				<< "  " << mesh.point(v_h).data()[1]
				<< "  " << mesh.point(v_h).data()[2] << std::endl;

			deform->shapeOpt1(std::vector<OpenMesh::VertexHandle>(1, v_h), std::vector<OpenMesh::Vec3d>(1, mesh.point(v_h) + dis), op);*/
		}
		

		//// arap
		//if (pickVertexValid())
		//{
		//	auto v_h = mesh.vertex_handle(lastestVertex);
		//	//arap_move_handle(v_h);
		//	std::cout << "Ori Position:" << mesh.point(v_h).data()[0]
		//		<< "  " << mesh.point(v_h).data()[1] 
		//		<< "  " << mesh.point(v_h).data()[2] << std::endl;
		//}
	}
	else if (edit_mode_ == INIT_SPHERE)
	{
		pick_face(_event->x(), _event->y());
		if (pickFaceValid() && seg.is_valid())
		{
			is_init_sphere[mesh.property(seg, mesh.face_handle(lastestFace))]
				= !is_init_sphere[mesh.property(seg, mesh.face_handle(lastestFace))];
		}
	}
	else if (edit_mode_ == INIT_CYLINDER)
	{
		pick_face(_event->x(), _event->y());
		if (pickFaceValid() && seg.is_valid())
		{
			is_init_cylinder[mesh.property(seg, mesh.face_handle(lastestFace))]
				= !is_init_cylinder[mesh.property(seg, mesh.face_handle(lastestFace))];
		}
	}
	updateGL();
}

void InteractiveViewerWidget::mouseMoveEvent(QMouseEvent *_event)
{
	if (edit_mode_ == TRANS)
	{
		MeshViewerWidget::mouseMoveEvent(_event);
	}
	else if (edit_mode_ == SEG_DRAW && mouse_mode == Qt::LeftButton)
	{
		pick_point(_event->x(), _event->y());
		pick_vertex(_event->x(), _event->y());
		if (pickVertexValid())
		{
			end_point = mesh.vertex_handle(lastestVertex);
			if (path_.size() != 0)
			{
				dijkstra_->ComputePath(start_point, end_point);
				cur_path_ = dijkstra_->return_path();
			}
			else
			{
				init_point = end_point;
				start_point = end_point;
			}
		}
	}
	else if (edit_mode_ == FACE_SELECT)
	{
		select_end_ = _event->pos();
		compute_selected_faces();
	}
	else if (edit_mode_ == VERT_SELECT)
	{
		select_end_ = _event->pos();
		compute_selected_vertices();
	}
	else if (edit_mode_ == VERT_DRAG)
	{
		if (pickVertexValid())
		{
			//move_point_based_lastVertex(_event->x(), _event->y());
			//auto v_h = mesh.vertex_handle(lastestVertex);
			//OpenMesh::Vec3d ori_position = mesh.point(v_h);
			/*OpenMesh::Vec3d new_position(selectedPoint[0], selectedPoint[1], selectedPoint[2]);*/
			//OpenMesh::Vec3d displacement = new_position - ori_position;

			


			/*std::cout << "New Position:" << new_position.data()[0]
				<< "  " << new_position.data()[1]
				<< "  " << new_position.data()[2] << std::endl;
			arap_move_position(v_h, new_position);*/

			/*std::vector<OpenMesh::VertexHandle> v_handle_vector = fixedVertices;
			std::vector<OpenMesh::Vec3d> new_position_vector = fixedPositions;


			v_handle_vector.push_back(v_h);
			new_position_vector.push_back(new_position);
			for (OpenMesh::VertexHandle vv_h:mesh.vv_range(v_h))
			{
				v_handle_vector.push_back(vv_h);
				new_position_vector.push_back(displacement + mesh.point(vv_h));
			}
			

			std::cout << "Pick Vert:" << v_h.idx() << std::endl;

			back = mesh;
			MeshViewerWidget::arap_move(v_handle_vector, new_position_vector);*/
		}
	}
	updateGL();
}

void InteractiveViewerWidget::mouseReleaseEvent(QMouseEvent *_event)
{
	if (edit_mode_ == TRANS)
	{
		MeshViewerWidget::mouseReleaseEvent(_event);
	}
	else if (edit_mode_ == FACE_SELECT)
	{
		select_start_ = select_end_;
	}
	else if (edit_mode_ == SEG_DRAW)
	{
		if (pickVertexValid())
		{
			start_point = mesh.vertex_handle(lastestVertex);
			path_.push_back(cur_path_);
			updateGL();
		}
	}
	else if (edit_mode_ == VERT_SELECT)
	{
		select_start_ = select_end_;
	}
	else if (edit_mode_ == VERT_DRAG)
	{
		if (pickVertexValid())
		{
			move_point_based_lastVertex(_event->x(), _event->y());
			auto v_h = segment.curMesh().vertex_handle(lastestVertex);
			OpenMesh::Vec3d ori_position = segment.curMesh().point(v_h);
			OpenMesh::Vec3d new_position(selectedPoint[0], selectedPoint[1], selectedPoint[2]);
			OpenMesh::Vec3d displacement = new_position - ori_position;

			
			std::cout << "Displacement:" << displacement.data()[0]
				<< "  " << displacement.data()[1]
				<< "  " << displacement.data()[2] << std::endl;
			
			//moveVertices.clear();
			//moveVertices.push_back(v_h);
			
			moveVec = displacement;

			/*std::vector<OpenMesh::VertexHandle> v_handle_vector = fixedVertices;
			std::vector<OpenMesh::Vec3d> new_position_vector = fixedPositions;


			v_handle_vector.push_back(v_h);
			new_position_vector.push_back(new_position);
			for (OpenMesh::VertexHandle vv_h:mesh.vv_range(v_h))
			{
				v_handle_vector.push_back(vv_h);
				new_position_vector.push_back(displacement + mesh.point(vv_h));
			}


			std::cout << "Pick Vert:" << v_h.idx() << std::endl;

			back = mesh;
			MeshViewerWidget::arap_move(v_handle_vector, new_position_vector);*/
		}
	}
}

void InteractiveViewerWidget::wheelEvent(QWheelEvent* _event)
{
	MeshViewerWidget::wheelEvent(_event);
}

void InteractiveViewerWidget::dragEnterEvent(QDragEnterEvent* event)
{
	if( event->mimeData()->hasFormat("text/uri-list") )
	{
		event->acceptProposedAction();
	}
}

void InteractiveViewerWidget::dropEvent(QDropEvent* event)
{
	QList<QUrl> urls = event->mimeData()->urls();
	if( urls.isEmpty() )
		return;
	QString fileName = urls.first().toLocalFile();
	if (fileName.isEmpty())
		return;

	if( fileName.endsWith(".off") || fileName.endsWith(".obj") || fileName.endsWith(".stl") || fileName.endsWith(".ply"))
	{
		if( openMesh(fileName.toLocal8Bit()))
		{
			emit(loadMeshOK(true,fileName));
			setDrawMode(FLAT_POINTS);
			setEditMode(TRANS);

			std::string filename_s(fileName.toLocal8Bit());
			int middle = filename_s.find_last_of('/');
			std::string name = filename_s.substr(middle);
			std::string folder = filename_s.substr(0, middle);
			filename_s = folder.substr(0, folder.find_last_of('/')) + "/newseg" + name.substr(0, name.find_last_of('.')) + ".seg";
			if (Load_Seg(filename_s.data()))
			{
				setDrawMode(InteractiveViewerWidget::SEGMENT);
			}
		}
		else
		{
			emit(loadMeshOK(false,"No Mesh"));
		}
	}
}

void InteractiveViewerWidget::pick_vertex(int x,int y)
{
	int r = find_vertex_using_selected_point();
	lastestVertex = r;
	//printf("Select Vertex : %d\n", r);
	/*updateGL();*/
}
void InteractiveViewerWidget::pick_face(int x,int y)
{
	int desiredFace = find_face_using_selected_point();
	// if(desiredFace < 0) return;
	lastestFace = desiredFace;
	updateGL();
}
void InteractiveViewerWidget::pick_edge(int x,int y)
{
	int desiredEdge = find_edge_using_selected_point();
	if(desiredEdge < 0) return;
	lastestEdge = desiredEdge;
	printf("Select Edge : %d\n", desiredEdge);
	updateGL();
}
void InteractiveViewerWidget::pick_point(int x,int y)
{
	GLint viewport[4];
	glGetIntegerv(GL_VIEWPORT, viewport);
	GLdouble winX = double(x);
	GLdouble winY = double( height() - y );
	GLfloat winZ = 0.0;
	glReadPixels((int)winX, (int)winY, 1, 1, GL_DEPTH_COMPONENT, GL_FLOAT, &winZ);
	gluUnProject(winX, winY, (GLdouble)winZ, &ModelViewMatrix[0], &ProjectionMatrix[0], viewport, &selectedPoint[0], &selectedPoint[1], &selectedPoint[2]);
}

void InteractiveViewerWidget::move_point_based_lastVertex(int x,int y)
{
	if(lastestVertex<0 || lastestVertex>=mesh.n_vertices())
	{
		return;
	}
	GLint viewport[4];
	glGetIntegerv(GL_VIEWPORT, viewport);
	GLdouble winX = 0.0;
	GLdouble winY = 0.0;
	GLdouble winZ = 0.0;
	OpenMesh::Vec3d p = segment.curMesh().point(segment.curMesh().vertex_handle(lastestVertex));
	gluProject(p[0], p[1], p[2],  &ModelViewMatrix[0], &ProjectionMatrix[0], viewport, &winX, &winY, &winZ);
	
	gluUnProject((GLdouble)(x), (GLdouble)( height() - y ), winZ,  &ModelViewMatrix[0], &ProjectionMatrix[0], viewport, &selectedPoint[0], &selectedPoint[1], &selectedPoint[2]);
}

int InteractiveViewerWidget::find_vertex_using_selected_point()
{
	ANNpoint tp = annAllocPt(3); tp[0] = selectedPoint[0]; tp[1] = selectedPoint[1]; tp[2] = selectedPoint[2];
	ANNidxArray nnIdx = new ANNidx[1]; ANNdistArray dists = new ANNdist[1];
	kdTree->annkSearch(tp, 1, nnIdx, dists);
	if (dists[0] > max_dists)
	{
		nnIdx[0] = -1;
	}

	return nnIdx[0];
}

int InteractiveViewerWidget::find_face_using_selected_point()
{
	int desiredFace = -1;
	int rv = find_vertex_using_selected_point();
	if (rv>= 0 && rv < mesh.n_vertices())
	{
		Mesh::VertexFaceIter vf_it = mesh.vf_iter(mesh.vertex_handle(rv));
		//double minLen = 10*radius();
		std::vector<OpenMesh::Vec3d> tri_p(3); int tri_count = 0;
		Mesh::Point resultP(selectedPoint[0], selectedPoint[1], selectedPoint[2]);
		for (vf_it; vf_it.is_valid(); ++vf_it)
		{
			tri_count = 0;
			for (Mesh::FaceVertexIter fv_it = mesh.fv_iter(*vf_it); fv_it.is_valid(); ++fv_it)
			{
				tri_p[tri_count] = mesh.point(*fv_it); ++tri_count;
			}
			if (check_in_triangle_face(tri_p, resultP))
			{
				desiredFace = vf_it->idx(); break;
			}
		}
		if (desiredFace < 0)
		{
			for (Mesh::FaceIter f_it = mesh.faces_begin(); f_it != mesh.faces_end(); ++f_it)
			{
				tri_count = 0;
				for (Mesh::FaceVertexIter fv_it = mesh.fv_iter(*f_it); fv_it.is_valid(); ++fv_it)
				{
					tri_p[tri_count] = mesh.point(*fv_it); ++tri_count;
				}
				if (check_in_triangle_face(tri_p, resultP))
				{
					desiredFace = f_it->idx(); break;
				}
			}
		}
	}
	
	return  desiredFace;
}

int InteractiveViewerWidget::find_edge_using_selected_point()
{
	int desiredFace = find_face_using_selected_point(); if(desiredFace < 0) return -1;
	Mesh::FaceHandle fh = mesh.face_handle(desiredFace);
	double min_len= 1e30; int desiredEdge = -1;
	Mesh::Point resultP(selectedPoint[0],selectedPoint[1],selectedPoint[2]);
	for(Mesh::FaceHalfedgeIter fhe_it = mesh.fh_iter(fh); fhe_it.is_valid(); ++fhe_it)
	{
		OpenMesh::Vec3d s = mesh.point( mesh.from_vertex_handle(*fhe_it) );
		OpenMesh::Vec3d e = mesh.point( mesh.to_vertex_handle(*fhe_it) );
		double dis = ((resultP - s) % (resultP - e)).norm() / (s - e).norm();
		if(dis < min_len){ min_len = dis; desiredEdge = mesh.edge_handle(*fhe_it).idx(); }
	}
	
	return desiredEdge;
}

void InteractiveViewerWidget::buildIndex()
{
	if(mesh.n_vertices() == 0)
		return;

	Mesh::VertexIter v_it(mesh.vertices_begin());
	Mesh::VertexIter v_end(mesh.vertices_end());
	Mesh::Point p;
	unsigned nv = mesh.n_vertices();
	ANNpointArray dataPts = annAllocPts(nv, 3);
	int count = 0;
	for(; v_it != v_end; ++v_it)
	{
		p = mesh.point(*v_it);
		dataPts[count][0] = p[0]; dataPts[count][1] = p[1]; dataPts[count][2] = p[2];
		++count;
	}

	if(kdTree) delete kdTree;
	kdTree = new ANNkd_tree(dataPts, nv, 3);
}

void InteractiveViewerWidget::draw_picked_seg()
{
	const Mesh &tmp = segment.curMesh();
	if (tmp.faces_empty()) return;
	glColor3f(1.0, 0.5, 1.0);
	glLineWidth(1);
	glBegin(GL_TRIANGLES);
	for (auto f : tmp.faces())
	{
		for (OpenMesh::VertexHandle fv : tmp.fv_range(f))
		{
			glVertex3dv(tmp.point(fv).data());
		}
	}
	glEnd();
}

void InteractiveViewerWidget::draw_picked_bound()
{
	//glEnable(GL_LIGHTING);
	float color[2][3] = { 1,0,0, 0,0,1 };
	int count = 0;
	glLineWidth(20);
	
	for (size_t i = 0; i < segment.n_contours(); i++)
	{
		if (!contour_pick_[i]) continue;
		glColor3f(color[count][0], color[count][1], color[count][2]);
		glBegin(GL_LINES);
		for (auto he : segment.contour(i))
		{
			auto from_v = segment.curMesh().from_vertex_handle(he);
			auto to_v = segment.curMesh().to_vertex_handle(he);
			glVertex3dv(segment.curMesh().point(from_v).data());
			glVertex3dv(segment.curMesh().point(to_v).data());
		}
		glEnd();
		count++;
	}
}

void InteractiveViewerWidget::draw_selected_vertices()
{
	glPointSize(5);
	glColor3f(1.0, 0.0, 0.0);
	glBegin(GL_POINTS);
	for (auto v_h : selectedVertices)
	{
		glVertex3dv(segment.curMesh().point(v_h).data());
	}
	glEnd();
}

void InteractiveViewerWidget::draw_selected_faces()
{
	glColor3f(1.0, 0.0, 0.0);
	glBegin(GL_TRIANGLES);
	for (auto f_h : selectedFaces)
	{
		for (auto v_h :mesh.fv_range(f_h))
		{
			glVertex3dv(mesh.point(v_h).data());
		}
	}
	glEnd();
}

void InteractiveViewerWidget::draw_seg_lines()
{
	glPointSize(5);
	glColor3f(0.0, 0.0, 0.0);
	glBegin(GL_POINTS);
	if (path_.size() >= 1)
	{
		glVertex3dv(mesh.point(start_point).data());
		glVertex3dv(mesh.point(end_point).data());
	}

	for (size_t i = 1; i < path_.size(); i++)
	{
		glVertex3dv(mesh.point(mesh.from_vertex_handle(path_[i].front())).data());
		glVertex3dv(mesh.point(mesh.to_vertex_handle(path_[i].back())).data());
	}
	glEnd();

	glLineWidth(3);
	glColor3f(1.0, 0.0, 0.0);
	glBegin(GL_LINES);
	for (size_t i = 1; i < path_.size(); i++)
	{
		for (auto he_h : path_[i])
		{
			auto to_v = mesh.to_vertex_handle(he_h);
			auto from_v = mesh.from_vertex_handle(he_h);
			glVertex3dv(mesh.point(to_v).data());
			glVertex3dv(mesh.point(from_v).data());
		}
	}
	glEnd();

	glLineWidth(3);
	glColor3f(0.0, 0.0, 1.0);
	glBegin(GL_LINES);
	for (auto he_h : cur_path_)
	{
		auto to_v = mesh.to_vertex_handle(he_h);
		auto from_v = mesh.from_vertex_handle(he_h);
		glVertex3dv(mesh.point(to_v).data());
		glVertex3dv(mesh.point(from_v).data());
	}
	glEnd();
}

void InteractiveViewerWidget::draw_fixed_vertices()
{
	glPointSize(5);
	glColor3f(0.0, 0.0, 1.0);
	glBegin(GL_POINTS);
	for (auto v_h : fixedVertices)
	{
		glVertex3dv(segment.curMesh().point(v_h).data());
	}
	glEnd();
}

void InteractiveViewerWidget::draw_scene(int drawmode)
{
	if (!mesh.n_vertices()) { return; }
	glViewport(0, 0, width(), height());
	glMatrixMode(GL_PROJECTION);
	glLoadMatrixd(&ProjectionMatrix[0]);
	glMatrixMode(GL_MODELVIEW);
	glLoadMatrixd(&ModelViewMatrix[0]);

	emit draw_from_out_signal();
	//draw select vertex, face, edge.
	glDisable(GL_LIGHTING);
	glDisable(GL_TEXTURE_2D);

	switch (edit_mode_)
	{
	case TRANS:
	case SEG_PICK:
	case BOUND_PICK:
		glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
		draw_picked_seg();
		glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
		draw_picked_bound();
		break;
	case SEG_CHANGE:
		break;
	case SEG_DRAW:
		draw_seg_lines();
		break;
	case VERT_SELECT:
		glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
		draw_picked_seg();
		glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
		draw_selected_vertices();
		draw_fixed_vertices();
		break;
	case FACE_SELECT:
		draw_selected_faces();
		break;
	case VERT_DRAG:
		glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
		draw_picked_seg();
		glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
		draw_selected_vertices();
		draw_fixed_vertices();
		break;
	case INIT_SPHERE:
		draw_init_sphere();
		break;
	case INIT_CYLINDER:
		draw_init_cylinder();
		break;
	default:
		break;
	}

	draw_scene_mesh(drawmode);
}

void InteractiveViewerWidget::draw_init_sphere() {
	glColor3f(0.0, 0.0, 1.0);
	glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
	glBegin(GL_TRIANGLES);
	for (auto f_h : mesh.faces())
	{
		if (!is_init_sphere[mesh.property(seg, f_h)]) continue;
		
		for (auto v_h : mesh.fv_range(f_h))
		{
			glVertex3dv(mesh.point(v_h).data());
		}
	}
	glEnd();
	glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
}

void InteractiveViewerWidget::draw_init_cylinder() {
	glColor3f(0.0, 1.0, 0.0);
	glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
	glBegin(GL_TRIANGLES);
	for (auto f_h : mesh.faces())
	{
		if (!is_init_cylinder[mesh.property(seg, f_h)]) continue;

		for (auto v_h : mesh.fv_range(f_h))
		{
			glVertex3dv(mesh.point(v_h).data());
		}
	}
	glEnd();
	glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
}

void InteractiveViewerWidget::compute_selected_faces()
{
	selectedFaces.clear();

	GLint viewport[4];
	glGetIntegerv(GL_VIEWPORT, viewport);
	GLdouble winX = 0.0;
	GLdouble winY = 0.0;
	GLdouble winZ = 0.0;
	for (auto f_h : mesh.faces())
	{
		int cont = 0;
		for (auto v_h : mesh.fv_range(f_h))
		{
			OpenMesh::Vec3d p = mesh.point(v_h);
			gluProject(p[0], p[1], p[2], &ModelViewMatrix[0], &ProjectionMatrix[0], viewport, &winX, &winY, &winZ);
			if (is_in_rect(winX, winY)) {
				cont++;
			}
		}
		if (cont == 3)
		{
			selectedFaces.push_back(f_h);
		}
	}
}

//---------------------------------------------------------------

void InteractiveViewerWidget::build_curMesh_tree()
{
	if (segment.curMesh().faces_empty()) return;

	int cont;

	Mesh::Point p;
	ANNpointArray dataPts = annAllocPts(segment.curMesh().n_vertices(), 3);

	cont = 0;
	for (auto v : segment.curMesh().vertices())
	{
		p = segment.curMesh().point(v);
		dataPts[cont][0] = p[0]; dataPts[cont][1] = p[1]; dataPts[cont][2] = p[2];
		++cont;
	}

	if (curMesh_kdTree) delete curMesh_kdTree;
	curMesh_kdTree = new ANNkd_tree(dataPts, segment.curMesh().n_vertices(), 3);
}

int InteractiveViewerWidget::pick_curMesh_idx()
{
	ANNpoint tp = annAllocPt(3); tp[0] = selectedPoint[0]; tp[1] = selectedPoint[1]; tp[2] = selectedPoint[2];
	ANNidxArray nnIdx = new ANNidx[1]; ANNdistArray dists = new ANNdist[1];
	curMesh_kdTree->annkSearch(tp, 1, nnIdx, dists);
	if (dists[0] > max_dists)
	{
		nnIdx[0] = -1;
	}

	return nnIdx[0];
}

//---------------------------------------------------------------

void InteractiveViewerWidget::build_bound_tree()
{
	treeidx_to_bound.clear();
	if (segment.n_contours() == 0) return;

	int nv = 0;
	int cont;

	// idx_to_bound
	cont = 0;
	std::unordered_map<VH, int> idx_to_bound;
	for (int i=0; i<segment.n_contours(); i++)
	{
		auto c = segment.contour(i);
		nv += c.size();
		for (auto he:c)
		{
			idx_to_bound.insert(std::make_pair(segment.curMesh().to_vertex_handle(he), cont));
		}
		cont++;
	}

	Mesh::Point p;
	ANNpointArray dataPts = annAllocPts(nv, 3);

	cont = 0;
	for (int i = 0; i < segment.n_contours(); i++)
	{
		auto c = segment.contour(i);
		for (auto he : c)
		{
			p = segment.curMesh().point(segment.curMesh().to_vertex_handle(he));
			dataPts[cont][0] = p[0]; dataPts[cont][1] = p[1]; dataPts[cont][2] = p[2];

			treeidx_to_bound.insert(std::make_pair(cont, idx_to_bound[segment.curMesh().to_vertex_handle(he)]));

			++cont;
		}
	}

	if (contour_kdTree) delete contour_kdTree;
	contour_kdTree = new ANNkd_tree(dataPts, nv, 3);
}

int InteractiveViewerWidget::pick_contour()
{
	ANNpoint tp = annAllocPt(3); tp[0] = selectedPoint[0]; tp[1] = selectedPoint[1]; tp[2] = selectedPoint[2];
	ANNidxArray nnIdx = new ANNidx[1]; ANNdistArray dists = new ANNdist[1];
	contour_kdTree->annkSearch(tp, 1, nnIdx, dists);
	if (dists[0] > max_dists)
	{
		nnIdx[0] = -1;
	}

	return nnIdx[0];
}

void InteractiveViewerWidget::compute_selected_vertices()
{
	selectedVertices.clear();

	GLint viewport[4];
	glGetIntegerv(GL_VIEWPORT, viewport);
	GLdouble winX = 0.0;
	GLdouble winY = 0.0;
	GLdouble winZ = 0.0;
	for (auto v_h : segment.curMesh().vertices())
	{
		OpenMesh::Vec3d p = segment.curMesh().point(v_h);
		gluProject(p[0], p[1], p[2], &ModelViewMatrix[0], &ProjectionMatrix[0], viewport, &winX, &winY, &winZ);
		if (is_in_rect(winX, winY)) {
			selectedVertices.push_back(v_h);
		}
	}
}

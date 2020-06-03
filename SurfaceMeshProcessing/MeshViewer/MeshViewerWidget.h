#ifndef MESHPROCESSING_MESHVIEWERWIDGET_H
#define MESHPROCESSING_MESHVIEWERWIDGET_H

#include <iostream>
#include <QString>
#include <QMessageBox>
#include <QFileDialog>
#include <Eigen/Dense>
#include "QGLViewerWidget.h"
#include "MeshDefinition.h"
#include <algorithm/SegmentMesh.hpp>
#include <algorithm/Deformation.hpp>
#include "Mesh_Ope.h"

enum DeformKind {
	spherizeOpen,
	cylinderize,
	planePolygon,
	planeFree,
	spherePolygon,
	smoothControl,
	featureControl,
	dragVertices
};

class MeshViewerWidget : public QGLViewerWidget 
{
	Q_OBJECT
public:
	MeshViewerWidget(QWidget* parent = 0);
	MeshViewerWidget(QGLFormat& _fmt, QWidget* _parent);
	~MeshViewerWidget();
public:
	bool openMesh(const char* filename);
	void initMesh();
	bool saveMesh(const char* filename);
	bool saveScreen(const char* filePath);
	void updateMesh();
	virtual void clearMesh();

	void printBasicMeshInfo();


	bool Load_Seg(const char* filename);
	bool Output_Segment(const char * filename);
	
	void ModeMaterial(int& mode);

	void Load_Tex(const char* filename);

signals:
	void loadMeshOK(bool,QString);

public slots:
	void shine(int i);
	void specular(int i);
	void deformKind(int i);
	void sliderPara(int i);
	void run_algorithm();
	void refresh_mesh();
	void assemble_mesh();

protected:
	void updateMeshCenter(); // used by update_mesh().
	void updateMeshNormals(); // used by update_mesh().


protected:
	virtual void draw_scene(int drawmode);
	void draw_scene_mesh(int drawmode);
	
private:
	void draw_mesh_wireframe();
	void draw_mesh_solidflat();
	void draw_seg();
	void draw_tex();


protected:
	bool first_init;
	OpenMesh::Vec3d bbMin;
	OpenMesh::Vec3d bbMax;

	// Material
	int shine_num_;
	int specular_num_;

	Mesh mesh;
	QImage tex;
	std::vector<bool> is_init_sphere;
	std::vector<bool> is_init_cylinder;
	// seg
	OpenMesh::FPropHandleT<int> seg;
	int seg_num;
	void initSeg();

	
	Mesh_Ope* ope;


	SegmentMesh segment;
	std::vector<bool> contour_pick_;
	
	//opParas op;
	DeformKind kind = spherizeOpen;
	CurvatureDeformation *deform;

	std::vector<OpenMesh::VertexHandle> fixedVertices;
	std::vector<OpenMesh::Vec3d> fixedPositions;
	std::vector<OpenMesh::VertexHandle> moveVertices;
	OpenMesh::Vec3d moveVec;

	// Pick_idx
	std::vector<int> pick_segs;
	int pick_bound_;
	int cur_slider_ = 20;


	void set_cur_seg(int i, bool add);


	/*void arap_move(std::vector<OpenMesh::VertexHandle> v_handle, std::vector<OpenMesh::Vec3d> v_position) {
		mesh.assign(back);
		deform->ARAP(v_handle, v_position, op);
	}*/
};

#endif // MESHPROCESSING_MESHVIEWERWIDGET_H

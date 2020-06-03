#include <iostream>
#include <algorithm>
#include <cmath>
#include <fstream>
#include <qapplication.h>

#include <OpenMesh/Core/Utils/vector_cast.hh>
#include <OpenMesh/Core/IO/MeshIO.hh>

#include "MeshViewerWidget.h"
#include "Common/CommonDefinitions.h"

using namespace Qt;

const int Color_N = 11;
const GLfloat COLOR[Color_N][3] = {
	{0.9020,    0.8941,    0.7529},
	{0.9059,    0.7922,    0.4980},
	{1.0000,    0.8431,    0.7333},
	{1.0000,    0.7020,    0.6588},
	{0.5098,    0.2235,    0.2078},
	{0.6235,    0.4902,    0.3137},
	{0.6745,    0.3176,    0.0941},
	{0.8706,    0.6118,    0.3255},
	{0.7882,    0.7294,    0.5137},
	{0.6314,    0.0902,    0.0824},
	{0.8706,    0.4902,    0.1725},
};

const std::string savefilename = "E:/untexture_mesh/savefile/editOrder.txt";
const std::string saveDrag = "E:/untexture_mesh/savefile/Drag/";
int dragcount = 0;
std::ofstream savefstream;
std::string enumStrings[8] = {
	"spherizeOpen",
	"cylinderize",
	"planePolygon",
	"planeFree",
	"spherePolygon",
	"smoothControl",
	"featureControl",
	"dragVertices"
};

MeshViewerWidget::MeshViewerWidget(QWidget* parent)
 : QGLViewerWidget(parent), segment(mesh, seg)
{
	mesh.request_vertex_status();
	mesh.request_edge_status();
	mesh.request_face_status();

	mesh.request_face_normals();
	mesh.request_vertex_normals();
	first_init = true;

	shine_num_ = 20;
	specular_num_ = 0;

	ope = NULL;
	deform = NULL;

	savefstream.open(savefilename);
}

MeshViewerWidget::MeshViewerWidget(QGLFormat& _fmt, QWidget* _parent)
: QGLViewerWidget(_fmt, _parent), segment(mesh, seg)
{
	mesh.request_vertex_status();
	mesh.request_edge_status();
	mesh.request_face_status();

	mesh.request_face_normals();
	mesh.request_vertex_normals();
	first_init = true;

	shine_num_ = 20;
	specular_num_ = 0;

	ope = NULL;
	deform = NULL;

	savefstream.open(savefilename);
}

MeshViewerWidget::~MeshViewerWidget()
{
	if (ope) delete ope;
	if (deform) delete deform;
	savefstream.close();
}

void MeshViewerWidget::updateMeshCenter()
{
	typedef Mesh::Point Point;
	Mesh::VertexIter vIt = mesh.vertices_begin();
	Mesh::VertexIter vEnd = mesh.vertices_end();
	bbMin = bbMax = mesh.point(*vIt);

	size_t count = 0;
	for (; vIt != vEnd; ++vIt, ++count)
	{
		bbMin = bbMin.minimize(mesh.point(*vIt));
		bbMax = bbMax.maximize(mesh.point(*vIt));
	}

	Mesh::EdgeIter e_it = mesh.edges_begin();
	Mesh::EdgeIter e_end = mesh.edges_end();
	double aveLen = 0.0; double maxLen = 0.0; double minLen = mesh.calc_edge_length(*e_it);
	double e_len = 0.0;
	for (; e_it != e_end; ++e_it)
	{
		double e_len = mesh.calc_edge_length(*e_it);
		if (e_len > maxLen)
		{
			maxLen = e_len;
		}
		else if (e_len < minLen)
		{
			minLen = e_len;
		}
		aveLen += e_len;
	}

	if (first_init)
	{
		set_scene_pos((bbMin + bbMax)*0.5, (bbMin - bbMax).norm()*0.5);
		first_init = false;
	}
	else
	{
		set_scene_pos((bbMin + bbMax)*0.5, (bbMin - bbMax).norm()*0.5);
	}

	printf("BoundingBox:\nX : [ %f , %f ]\n", bbMin[0], bbMax[0]);
	printf("Y : [ %f , %f ]\n", bbMin[1], bbMax[1]);
	printf("Z : [ %f , %f ]\n", bbMin[2], bbMax[2]);
	printf("Diag length of BBox : %f\n", (bbMax - bbMin).norm());
	printf("Edge Length : Max : %f; Min : %f; AVG : %f\n", maxLen, minLen, aveLen / mesh.n_edges());
}

void MeshViewerWidget::updateMeshNormals()
{
	mesh.update_face_normals();
	mesh.update_vertex_normals();
}

bool MeshViewerWidget::openMesh(const char* filename)
{
	clearMesh();
	mesh.request_halfedge_texcoords2D();
	OpenMesh::IO::Options opt_tex = OpenMesh::IO::Options::FaceTexCoord;
	bool read_OK = OpenMesh::IO::read_mesh( mesh, filename , opt_tex);

	printf("%s\n", filename);
	if ( read_OK )
	{
		initMesh();
		initSeg();
		return true;
	}
	return false;
}

void MeshViewerWidget::initMesh()
{
	mesh.request_vertex_status();
	mesh.request_edge_status();
	mesh.request_face_status();

	mesh.request_face_normals();
	mesh.request_vertex_normals();
	printBasicMeshInfo();
	updateMesh();
}

void MeshViewerWidget::printBasicMeshInfo()
{
	if (mesh.n_vertices() == 0)
		printf("No Mesh\n");

	printf("Information of the input mesh:\nVertex : %d;\nFace : %d;\nEdge : %d, HalfEdge : %d\n",
		mesh.n_vertices(),mesh.n_faces(),mesh.n_edges(),mesh.n_halfedges());;
}

bool MeshViewerWidget::Load_Seg(const char * filename)
{
	// Mesh
	if (mesh.n_faces() == 0)
	{
		std::cerr << "There is no face!" << std::endl;
		return false;
	}

	// Open File
	std::vector<int> seg_vector;
	std::ifstream srcFile(filename, std::ios::in);
	if (!srcFile) {
		std::cerr << "Error opening seg file." << std::endl;
		return false;
	}

	// Read Segmentation
	int x;
	while (srcFile >> x)
		seg_vector.push_back(x);
	srcFile.close();

	// Add seg Property
	if (seg_vector.size() == mesh.n_faces())
	{
		// Seg_num
		seg_num = *max_element(seg_vector.begin(), seg_vector.end()) + 1;
		std::cout << "Seg Num: " << seg_num << std::endl;

		// Add face_seg property
		for (auto f : mesh.faces())
		{
			mesh.property(seg, f) = seg_vector[f.idx()];
		}
		

		is_init_sphere.resize(seg_num, false);
		is_init_cylinder.resize(seg_num, false);
	}
	else
	{
		std::cerr << "Error load segmentation." << std::endl;
	}

	segment.update();
	return true;
}

void MeshViewerWidget::Load_Tex(const char * filename)
{
	QImage buf;
	buf.load(filename);
	tex = QGLWidget::convertToGLFormat(buf);

}

bool MeshViewerWidget::Output_Segment(const char * filename)
{
	std::ofstream ofile;
	ofile.open(filename);
	for (FH f : mesh.faces())
	{
		ofile << mesh.property(seg, f) << std::endl;
	}
	
	return true;
}

void MeshViewerWidget::ModeMaterial(int& mode)
{
	GLfloat mat_a[] = { COLOR[mode][0],  COLOR[mode][1],  COLOR[mode][2], 1.0f };
	GLfloat mat_d[] = { COLOR[mode][0],  COLOR[mode][1],  COLOR[mode][2], 1.0f };
	GLfloat mat_s[] = { 0.01f*specular_num_, 0.01f*specular_num_, 0.01f*specular_num_, 1.0f };
	GLfloat shine[] = { shine_num_*1.0f };

	glMaterialfv(GL_FRONT_AND_BACK, GL_AMBIENT, mat_a);
	glMaterialfv(GL_FRONT_AND_BACK, GL_DIFFUSE, mat_d);
	glMaterialfv(GL_FRONT_AND_BACK, GL_SPECULAR, mat_s);
	glMaterialfv(GL_FRONT_AND_BACK, GL_SHININESS, shine);
}


bool MeshViewerWidget::saveMesh(const char* filename)
{
	mesh.request_halfedge_texcoords2D();
	OpenMesh::IO::Options opt_tex = OpenMesh::IO::Options::FaceTexCoord;
	return OpenMesh::IO::write_mesh(mesh, filename, opt_tex);
}
bool MeshViewerWidget::saveScreen(const char* filePath)
{
	QImage image = grabFrameBuffer();
	image.save(filePath);
	return true;
}

void MeshViewerWidget::updateMesh()
{
	updateMeshCenter();
	updateMeshNormals();
}

void MeshViewerWidget::clearMesh()
{
	mesh.clear();
	updateGL();
}

void MeshViewerWidget::draw_scene(int drawmode)
{
	QFont Text_Font("Courier", 12);
	glViewport ( 0,0, width(),height());
	glMatrixMode( GL_PROJECTION );
	glLoadMatrixd( &ProjectionMatrix[0] );
	glMatrixMode( GL_MODELVIEW );
	glLoadMatrixd( &ModelViewMatrix[0] );

	draw_scene_mesh(drawmode);
}

void MeshViewerWidget::draw_scene_mesh(int drawmode)
{
	if(mesh.n_vertices() == 0) { return; }

	switch (drawmode)
	{
	case FLAT_POINTS:
		glEnable(GL_POLYGON_OFFSET_FILL);
		glPolygonOffset(1.5f,2.0f);
		glEnable(GL_LIGHTING);
		glShadeModel(GL_FLAT);
		draw_mesh_solidflat();
		glDisable(GL_POLYGON_OFFSET_FILL);
		//draw_meshpointset();
		glDisable(GL_LIGHTING);
		glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
		draw_mesh_wireframe();
		glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
		break;
	case SEGMENT:
		if (seg.is_valid() && seg_num > 0)
		{
			glEnable(GL_LIGHTING);
			draw_seg();
			glDisable(GL_LIGHTING);
		}
		break;
	case TEXTURE:
		draw_tex();
		break;
	default:
		break;
	}
}


void MeshViewerWidget::draw_mesh_wireframe()
{
	glLineWidth(1);
	glColor3f(0.0, 0.0, 0.25);
	Mesh::ConstFaceIter fIt(mesh.faces_begin()),
		fEnd(mesh.faces_end());
	Mesh::ConstFaceVertexIter fvIt;
	for (; fIt != fEnd; ++fIt)
	{
		fvIt = mesh.cfv_iter(*fIt);
		glBegin(GL_POLYGON);
		for (fvIt; fvIt.is_valid(); ++fvIt)
		{
			glVertex3dv(mesh.point(*fvIt).data());
		}
		glEnd();
	}

}

void MeshViewerWidget::draw_mesh_solidflat()
{
	Mesh::ConstFaceIter fIt(mesh.faces_begin()),
		fEnd(mesh.faces_end());
	Mesh::ConstFaceVertexIter fvIt;

	int mode = 0;
	ModeMaterial(mode);

	glBegin(GL_TRIANGLES);
	for (fIt; fIt != fEnd; ++fIt)
	{
		glNormal3dv(mesh.normal(*fIt).data());
		fvIt = mesh.cfv_iter(*fIt);
		for (fvIt; fvIt.is_valid(); ++fvIt)
		{
			glVertex3dv(mesh.point(*fvIt).data());
		}
	}
	glEnd();
}

void MeshViewerWidget::draw_seg()
{
	for (auto f : mesh.faces())
	{
		glNormal3dv(mesh.normal(f).data());
		int temp_seg = mesh.property(seg, f);

		if (std::find(pick_segs.begin(), pick_segs.end(), temp_seg) != pick_segs.end()) continue;

		int mode = temp_seg % Color_N;
		ModeMaterial(mode);

		glBegin(GL_TRIANGLES);
		for (OpenMesh::VertexHandle fv : mesh.fv_range(f))
		{
			glVertex3dv(mesh.point(fv).data());
		}
		glEnd();
	}
}

void MeshViewerWidget::draw_tex()
{
	GLuint tex_id;
	
	glGenTextures(1, &tex_id);
	glBindTexture(GL_TEXTURE_2D, tex_id);

	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);

	glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA, tex.width(), tex.height(), 0,
		GL_RGBA, GL_UNSIGNED_BYTE, tex.bits());

	//std::cout << tex.width() << "  " << tex.height() << "  " << tex.bits()[0] << "  " << tex.bits()[1] << std::endl;

	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
	glEnable(GL_TEXTURE_2D);
	glTexEnvf(GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_DECAL);
	glBindTexture(GL_TEXTURE_2D, tex_id);
	glTexEnvf(GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_DECAL);
	glBegin(GL_TRIANGLES);
	for (Mesh::FaceIter f_it = mesh.faces_begin(); f_it != mesh.faces_end(); ++f_it) {
		for (HEH he_h:mesh.fh_range(*f_it))
		{
			glTexCoord2dv(mesh.texcoord2D(he_h).data());
			glVertex3dv(mesh.point(mesh.to_vertex_handle(he_h)).data());
		}
	}
	glEnd();
	glDisable(GL_TEXTURE_2D);
}

void MeshViewerWidget::initSeg()
{
	mesh.add_property(seg);
	for (auto f : mesh.faces())
	{
		mesh.property(seg, f) = 0;
	}
	seg_num = 1;


	pick_segs.clear();
	pick_bound_ = -1;

	if (ope) delete ope;
	ope = new Mesh_Ope(mesh, seg);
	segment.update();

	is_init_sphere.resize(seg_num, false);
	is_init_cylinder.resize(seg_num, false);
}

void MeshViewerWidget::deformKind(int i)
{
	kind = DeformKind(i);
	if (pick_segs.empty()) return;
	contour_pick_ = std::vector<bool>(segment.n_contours(), false);
	if (deform)
	{
		delete deform;
		deform = NULL;
	}
	switch (kind)
	{
	case spherizeOpen:
	{
		std::cout << "Spherize Init...\n";
		deform = new Spherelize(segment.curMesh(), segment.hContours());
		deform->setCurvatureParameters(cur_slider_);
		break;
	}
	case cylinderize:
	{
		if (segment.n_contours() >= 2)
		{
			std::cout << "Cylinderize Init...\n";
			deform = new Cylinderize(segment.curMesh(), segment.hContours());
			deform->setCurvatureParameters(cur_slider_);
			contour_pick_[0] = contour_pick_[1] = true;
		}
		break;
	}
	case planePolygon:
	{
		std::cout << "PlanePolygon Init...\n";
		deform = new PlanePolygon(segment.curMesh(), segment.hContours());
		deform->setCurvatureParameters(cur_slider_);
		break;
	}
	case planeFree:
	{
		std::cout << "PlaneFree Init...\n";
		deform = new PlaneFree(segment.curMesh(), segment.hContours());
		deform->setCurvatureParameters(cur_slider_);
		break;
	}
	case spherePolygon:
	{
		std::cout << "SpherePolygon Init...\n";
		deform = new SpherePolygon(segment.curMesh(), segment.hContours());
		deform->setCurvatureParameters(cur_slider_);
		break;
	}
	case smoothControl:
	{
		std::cout << "SmoothControl Init...\n";
		deform = new SmoothControl(segment.curMesh(), segment.hContours());
		deform->setCurvatureParameters(cur_slider_);
		break;
	}
	case featureControl:
	{
		std::cout << "FeatureControl Init...\n";
		deform = new FeatureControl(segment.curMesh(), segment.hContours());
		deform->setCurvatureParameters(cur_slider_);
		break;
	}
	case dragVertices:
	{
		std::cout << "DragVertices Init...\n";
		deform = new DragVertices(segment.curMesh(), segment.hContours());
		deform->setCurvatureParameters(cur_slider_);
		break;
	}
	default:
		break;
	}

	updateGL();
	savefstream << "\n" << enumStrings[kind] << '\t' << "pick_segs : ";
	for (int i = 0; i < pick_segs.size(); i++) savefstream << pick_segs[i] << '\t';
}

void MeshViewerWidget::sliderPara(int i)
{
	std::cout << "Curvature: " << i << std::endl;
	std::cout << "Slider: " << i << std::endl;
	cur_slider_ = i;
	if (deform) deform->setCurvatureParameters(cur_slider_);
}

void MeshViewerWidget::shine(int i)
{
	std::cout << "shine: " << i << std::endl;
	shine_num_ = i;
}

void MeshViewerWidget::specular(int i)
{
	std::cout << "specular: " << i << std::endl;
	specular_num_ = i;
}

void MeshViewerWidget::run_algorithm()
{
	std::cout << "curvature deformation : ... \n";

	if (deform)
	{
		if (kind == DeformKind::cylinderize)
		{
			dynamic_cast<Cylinderize *>(deform)->selectContours(contour_pick_);
			int contourIds[2];
			int count = 0;
			for (int i = 0; i < contour_pick_.size(); i++)
			{
				if (contour_pick_[i])
				{
					if (count < 2) contourIds[count] = i;
					count++;
				}
			}
			if (count == 2)
			{
				savefstream << "selectContours : " << contourIds[0] << ' ' << contourIds[1] << '\t';
			}
		}
		if (kind == DeformKind::dragVertices)
		{
			Vectori achs;
			Vector3d dpos;
			Mesh &tmp = segment.curMesh();
			for (auto v : fixedVertices)
			{
				achs.push_back(v.idx());
				dpos.push_back(Vec3d(0, 0, 0));
			}
			/*for (auto v : tmp.vertices())
			{
				bool nearbound = tmp.is_boundary(v);
				for (auto vv : tmp.vv_range(v))
				{
					nearbound = nearbound || tmp.is_boundary(vv);
					if (nearbound) break;
				}
				if (nearbound)
				{
					achs.push_back(v.idx());
					dpos.push_back(Vec3d(0, 0, 0));
				}
			}*/
			for (auto v : moveVertices)
			{
				achs.push_back(v.idx());
				dpos.push_back(Vec3d(moveVec[0], moveVec[1], moveVec[2]));
			}

			/*std::ofstream fout1(saveDrag + "fixedHandle" + std::to_string(dragcount) + ".txt");
			if (fout1)
			{
				for (auto v : fixedVertices)
				{
					fout1 << v.idx() << std::endl;
				}
			}
			fout1.close();*/

			std::ofstream fout2(saveDrag + "movedHandle" + std::to_string(dragcount) + ".txt");
			if (fout2)
			{
				for (auto v : moveVertices)
				{
					fout2 << v.idx() << std::endl;
				}
			}
			fout2.close();

			std::ofstream fout3(saveDrag + "movedVector" + std::to_string(dragcount) + ".txt");
			if (fout3)
			{
				for (auto v : moveVertices)
				{
					fout3 << moveVec[0] << ' ' << moveVec[1] << ' ' << moveVec[2] << std::endl;
				}
			}
			fout3.close();

			savefstream << "dragNumber : " << dragcount << '\t';
			dragcount++;

			/*int idx;
			double vx, vy, vz;
			std::ifstream fin1("G:/final_backup/107_models/fixedHandle.txt");
			if (fin1)
			{
				while (fin1 >> idx)
				{
					achs.push_back(idx);
					dpos.push_back(Vec3d(0, 0, 0));
				}
			}
			fin1.close();

			std::ifstream fin2("G:/final_backup/107_models/movedHandle.txt");
			if (fin2)
			{
				while (fin2 >> idx)
				{
					achs.push_back(idx);
				}
			}
			fin2.close();
			
			std::ifstream fin3("G:/final_backup/107_models/movedVector.txt");
			if (fin3)
			{
				while (fin3 >> vx && fin3 >> vy && fin3 >> vz)
				{
					dpos.push_back(Vec3d(vx, vy, vz));
				}
			}
			fin3.close();*/

			deform->anchorPos(achs, dpos);
		}
		savefstream << "sliderValue : " << cur_slider_ << '\t' << "Time : " << deform->run() << '\t';
		segment.assembleMesh(deform->deformedPos());
		updateGL();
	}
}

void MeshViewerWidget::refresh_mesh()
{
	if (deform)
	{
		segment.refreshMesh(pick_segs);
		contour_pick_ = std::vector<bool>(segment.n_contours(), false);
		if (deform)
		{
			delete deform;
			deform = NULL;
		}
		switch (kind)
		{
		case spherizeOpen:
		{
			std::cout << "Spherize Init...\n";
			deform = new Spherelize(segment.curMesh(), segment.hContours());
			deform->setCurvatureParameters(cur_slider_);
			break;
		}
		case cylinderize:
		{
			if (segment.n_contours() >= 2)
			{
				std::cout << "Cylinderize Init...\n";
				deform = new Cylinderize(segment.curMesh(), segment.hContours());
				deform->setCurvatureParameters(cur_slider_);
				contour_pick_[0] = contour_pick_[1] = true;
			}
			break;
		}
		case planePolygon:
		{
			std::cout << "PlanePolygon Init...\n";
			deform = new PlanePolygon(segment.curMesh(), segment.hContours());
			deform->setCurvatureParameters(cur_slider_);
			break;
		}
		case planeFree:
		{
			std::cout << "PlaneFree Init...\n";
			deform = new PlaneFree(segment.curMesh(), segment.hContours());
			deform->setCurvatureParameters(cur_slider_);
			break;
		}
		case spherePolygon:
		{
			std::cout << "SpherePolygon Init...\n";
			deform = new SpherePolygon(segment.curMesh(), segment.hContours());
			deform->setCurvatureParameters(cur_slider_);
			break;
		}
		case smoothControl:
		{
			std::cout << "SmoothControl Init...\n";
			deform = new SmoothControl(segment.curMesh(), segment.hContours());
			deform->setCurvatureParameters(cur_slider_);
			break;
		}
		case featureControl:
		{
			std::cout << "FeatureControl Init...\n";
			deform = new FeatureControl(segment.curMesh(), segment.hContours());
			deform->setCurvatureParameters(cur_slider_);
			break;
		}
		case dragVertices:
		{
			std::cout << "DragVertices Init...\n";
			deform = new DragVertices(segment.curMesh(), segment.hContours());
			deform->setCurvatureParameters(cur_slider_);
			break;
		}
		default:
			break;
		}
		updateGL();
		std::cout << "refreshMesh\n";
		savefstream << "refreshMesh\t";
	}
}

void MeshViewerWidget::assemble_mesh()
{
	if (deform)
	{
		segment.updateMesh();
		updateGL();
		std::cout << "assembleMesh\n";
		savefstream << "assembleMesh\t";
	}
}

void MeshViewerWidget::set_cur_seg(int i, bool add) {
	std::cout << "Pick Seg: " << i << std::endl;
	std::vector<int>::iterator iteri = std::find(pick_segs.begin(), pick_segs.end(), i);
	if (iteri == pick_segs.end())
	{
		if (add) pick_segs.push_back(i);
		else pick_segs = std::vector<int>(1, i);
	}
	else
	{
		if(add) pick_segs.erase(iteri);
		else pick_segs = std::vector<int>(1, i);
	}
	if (deform)
	{
		delete deform;
		deform = NULL;
	}
	segment.select(pick_segs);
	contour_pick_ = std::vector<bool>(segment.n_contours(), false);
	if (pick_segs.empty()) return;
	switch (kind)
	{
	case spherizeOpen:
	{
		std::cout << "Spherize Init...\n";
		deform = new Spherelize(segment.curMesh(), segment.hContours());
		deform->setCurvatureParameters(cur_slider_);
		break;
	}
	case cylinderize:
	{
		if (segment.n_contours() >= 2)
		{
			std::cout << "Cylinderize Init...\n";
			deform = new Cylinderize(segment.curMesh(), segment.hContours());
			deform->setCurvatureParameters(cur_slider_);
			contour_pick_[0] = contour_pick_[1] = true;
		}
		break;
	}
	case planePolygon:
	{
		std::cout << "PlanePolygon Init...\n";
		deform = new PlanePolygon(segment.curMesh(), segment.hContours());
		deform->setCurvatureParameters(cur_slider_);
		break;
	}
	case planeFree:
	{
		std::cout << "PlaneFree Init...\n";
		deform = new PlaneFree(segment.curMesh(), segment.hContours());
		deform->setCurvatureParameters(cur_slider_);
		break;
	}
	case spherePolygon:
	{
		std::cout << "SpherePolygon Init...\n";
		deform = new SpherePolygon(segment.curMesh(), segment.hContours());
		deform->setCurvatureParameters(cur_slider_);
		break;
	}
	case smoothControl:
	{
		std::cout << "SmoothControl Init...\n";
		deform = new SmoothControl(segment.curMesh(), segment.hContours());
		deform->setCurvatureParameters(cur_slider_);
		break;
	}
	case featureControl:
	{
		std::cout << "FeatureControl Init...\n";
		deform = new FeatureControl(segment.curMesh(), segment.hContours());
		deform->setCurvatureParameters(cur_slider_);
		break;
	}
	case dragVertices:
	{
		std::cout << "DragVertices Init...\n";
		deform = new DragVertices(segment.curMesh(), segment.hContours());
		deform->setCurvatureParameters(cur_slider_);
		break;
	}
	default:
		break;
	}

	savefstream << "\n" << enumStrings[kind] << '\t' << "pick_segs : ";
	for (int i = 0; i < pick_segs.size(); i++) savefstream << pick_segs[i] << '\t';
}
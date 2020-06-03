#ifndef MESHPROCESSING_MAIN_VIEWWE_WIDGET_H
#define MESHPROCESSING_MAIN_VIEWWE_WIDGET_H

#include <QtGui>
#include <QString>
#include <QMessageBox>
#include <QFileDialog>
//main widget
#include "InteractiveViewerWidget.h"
#include "MeshParamDialog.h"

class MainViewerWidget : public QDialog
{
	Q_OBJECT

public:
	MainViewerWidget(QWidget* _parent = 0);
	~MainViewerWidget();

	void setDrawMode(int dm);
	void setEditMode(int dm);


public slots:
	void open_mesh_query();
	void save_mesh_query();
	void saveOpenGLScreen();

	virtual void clear_all_mesh();
	void LoadMeshFromInner(bool OK, QString fname);

	void print_info();

	void load_seg();

	void test() {
		open_mesh_gui("models/cat/12221_Cat_v1_l3_tri.obj");
		MeshViewer->Load_Seg("models/cat/test.seg");
		MeshViewer->Load_Tex("models/cat/Cat_diffuse.jpg");
		//MeshViewer->run_algorithm();
	}

	// Draw Seg
	void draw_undo() {
		MeshViewer->Draw_Undo();
	}
	void undo_finished() {
		MeshViewer->Undo_Finished();
	}
	void undo_seg_change() {
		MeshViewer->Undo_Seg_Change();
	}
	void finish_seg_draw() {
		MeshViewer->Finish_Seg_Draw();
	}

	void undo() {

	}
	
	void output_seg();
	void load_tex();
	void fix_vertices() {
		MeshViewer->Fix_Vertices();
	}

signals:
	void haveLoadMesh(QString filePath);
	void setEditMode_signal_main(int);
	void setDrawMode_signal_main(int);
	void setTexShow_main();
	void setSegShow_main();

protected:
	virtual void initViewerWindow();
	virtual void createParamIDialog();
	virtual void createViewerDialog();
	virtual void save_mesh_gui(QString fname);
	virtual void open_mesh_gui(QString fname);
	virtual void save_screen_gui(QString fname);

protected:
	bool LoadMeshSuccess;

private:
	InteractiveViewerWidget* MeshViewer;
	MeshParamDialog* MeshParam;
};

#endif
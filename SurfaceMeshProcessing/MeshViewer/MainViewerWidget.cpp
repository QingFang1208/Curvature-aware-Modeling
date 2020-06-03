#include "MainViewerWidget.h"

MainViewerWidget::MainViewerWidget(QWidget* _parent/* =0 */)
{
	initViewerWindow();
	LoadMeshSuccess = false;
}
MainViewerWidget::~MainViewerWidget()
{
}
void MainViewerWidget::setDrawMode(int dm)
{
	MeshViewer->setDrawMode(dm);
}

void MainViewerWidget::setEditMode(int dm)
{
	MeshViewer->setEditMode(dm);
}

void MainViewerWidget::open_mesh_query() {
	QString fileName = QFileDialog::getOpenFileName(this,
		tr("Open mesh file"),
		tr("../models/"),
		tr("OFF Files (*.off);;"
			"OBJ Files (*.obj);;"
			"PLY Files (*.ply);;"
			"STL Files (*.stl);;"
			"All Files (*)"));
	if (!fileName.isEmpty())
	{
		open_mesh_gui(fileName);
	}
}

void MainViewerWidget::save_mesh_query() {
	QString fileName = QFileDialog::getSaveFileName(this,
		tr("Save mesh file"),
		tr("../models/untitled.off"),
		tr("OFF Files (*.off);;"
			"OBJ Files (*.obj);;"
			"PLY Files (*.ply);;"
			"STL Files (*.stl);;"
			"All Files (*)"));
	if (!fileName.isEmpty())
	{
		save_mesh_gui(fileName);
	}
}

void MainViewerWidget::saveOpenGLScreen()
{
	QString fileName = QFileDialog::getSaveFileName(this,
		("Save screen as image file"),
		("../Results/untitled.png"),
		("PNG Files (*.png);;BMP Files (*.bmp);;JPG Files (*.jpg);;"
			"All Files (*)"));
	if (!fileName.isEmpty())
	{
		save_screen_gui(fileName);
	}
}

void MainViewerWidget::clear_all_mesh()
{
	if (LoadMeshSuccess)
	{
		LoadMeshSuccess = false;
		MeshViewer->clearMesh();
	}
}

void MainViewerWidget::LoadMeshFromInner(bool OK, QString fname)
{
	LoadMeshSuccess = OK;
	std::string filename_s(fname.toLocal8Bit());
	int middle = filename_s.find_last_of('/');
	std::string name = filename_s.substr(middle);
	std::string folder = filename_s.substr(0, middle);
	filename_s = folder.substr(0, folder.find_last_of('/')) + "/newseg" + name.substr(0, name.find_last_of('.')) + ".seg";
	if (MeshViewer->Load_Seg(filename_s.data()))
	{
		MeshViewer->setDrawMode(InteractiveViewerWidget::SEGMENT);
	}
	emit(haveLoadMesh(fname));
}

void MainViewerWidget::initViewerWindow()
{
	createParamIDialog();
	createViewerDialog();

	QHBoxLayout* main_layout = new QHBoxLayout();
	main_layout->addWidget(MeshParam, 2);
	main_layout->addWidget(MeshViewer, 7);
	this->setLayout(main_layout);

	connect(MeshViewer, SIGNAL(setEditMode_signal(int)), SIGNAL(setEditMode_signal_main(int)));
	connect(MeshViewer, SIGNAL(setDrawMode_signal(int)), SIGNAL(setDrawMode_signal_main(int)));

	connect(MeshParam, SIGNAL(PrintInfoSignal()), SLOT(print_info()));

	connect(MeshParam, SIGNAL(ShineSignal(int)), MeshViewer, SLOT(shine(int)));
	connect(MeshParam, SIGNAL(SpecularSignal(int)), MeshViewer, SLOT(specular(int)));
	connect(MeshParam, SIGNAL(runAlgorithm()), MeshViewer, SLOT(run_algorithm()));
	connect(MeshParam, SIGNAL(refreshMesh()), MeshViewer, SLOT(refresh_mesh()));
	connect(MeshParam, SIGNAL(assembleMesh()), MeshViewer, SLOT(assemble_mesh()));
	connect(MeshParam, SIGNAL(customizeSignal(int)), MeshViewer, SLOT(deformKind(int)));
	connect(MeshParam, SIGNAL(sliderSignal(int)), MeshViewer, SLOT(sliderPara(int)));
}

void MainViewerWidget::createParamIDialog()
{
	MeshParam = new MeshParamDialog();
}

void MainViewerWidget::createViewerDialog()
{
	QGLFormat glFormat;
	glFormat.setSampleBuffers(true);
	glFormat.setSamples(16);

	MeshViewer = new InteractiveViewerWidget(glFormat, NULL);
	MeshViewer->setAcceptDrops(true);
	connect(MeshViewer, SIGNAL(loadMeshOK(bool, QString)), this, SLOT(LoadMeshFromInner(bool, QString)));
}

void MainViewerWidget::open_mesh_gui(QString fname)
{
	if (fname.isEmpty() || !MeshViewer->openMesh(fname.toLocal8Bit()))
	{
		QString msg = "Cannot read mesh from file:\n '";
		msg += fname;
		msg += "'";
		QMessageBox::critical(NULL, windowTitle(), msg);
	}
	else
	{
		LoadMeshSuccess = true;
		MeshViewer->setDrawMode(InteractiveViewerWidget::FLAT_POINTS);
		MeshViewer->setEditMode(InteractiveViewerWidget::TRANS);
		emit(haveLoadMesh(fname));
	}
}

void MainViewerWidget::save_mesh_gui(QString fname)
{
	if (fname.isEmpty() || !MeshViewer->saveMesh(fname.toLocal8Bit()))
	{
		QString msg = "Cannot read mesh from file:\n '";
		msg += fname;
		msg += "'";
		QMessageBox::critical(NULL, windowTitle(), msg);
	}
}

void MainViewerWidget::save_screen_gui(QString fname)
{
	if (fname.isEmpty() || !MeshViewer->saveScreen(fname.toLocal8Bit()))
	{
		QString msg = "Cannot save image to file:\n '";
		msg += fname;
		msg += "'";
		QMessageBox::critical(NULL, windowTitle(), msg);
	}
}

void MainViewerWidget::print_info()
{
	MeshViewer->printBasicMeshInfo();
}

void MainViewerWidget::load_tex()
{
	QString fileName = QFileDialog::getOpenFileName(this,
		tr("Open texture file"),
		tr(""),
		tr("JPG Files (*.jpg);;"
		));
	if (!fileName.isEmpty())
	{
		MeshViewer->Load_Tex(fileName.toLocal8Bit());
		emit setTexShow_main();
	}
}

void MainViewerWidget::load_seg()
{
	QString fileName = QFileDialog::getOpenFileName(this,
		tr("Open Segmentation"),
		tr("../models/"),
		tr("seg Files (*.seg);;"));


	if (!fileName.isEmpty())
	{
		MeshViewer->Load_Seg(fileName.toLocal8Bit());
		emit setSegShow_main();
	}
}

void MainViewerWidget::output_seg()
{
	QString fileName = QFileDialog::getSaveFileName(this,
		tr("Save mesh file"),
		tr(""),
		tr("SEG Files (*.seg);;"));
	if (!fileName.isEmpty())
	{
		if (fileName.isEmpty() || !MeshViewer->Output_Segment(fileName.toLocal8Bit()))
		{
			QString msg = "Cannot read mesh from file:\n '";
			msg += fileName;
			msg += "'";
			QMessageBox::critical(NULL, windowTitle(), msg);
		}
	}
}

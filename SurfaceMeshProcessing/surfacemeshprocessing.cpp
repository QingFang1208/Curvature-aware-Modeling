#include "surfacemeshprocessing.h"
#include "MeshViewer/MainViewerWidget.h"

SurfaceMeshProcessing::SurfaceMeshProcessing(QWidget *parent)
	: QMainWindow(parent)
{
	viewer = new MainViewerWidget(this);
	setCentralWidget(viewer);
	initWindow();
	createActions();
	createMenus();
	createToolBars();
	createStatusBar();
}

SurfaceMeshProcessing::~SurfaceMeshProcessing()
{
}

void SurfaceMeshProcessing::initWindow()
{
	connect(viewer, SIGNAL(setEditMode_signal_main(int)), this, SLOT(setEditMode_slot(int)));
	connect(viewer, SIGNAL(setDrawMode_signal_main(int)), this, SLOT(setDrawMode_slot(int)));
	connect(viewer, SIGNAL(setTexShow_main()), this, SLOT(TexShow()));
	connect(viewer, SIGNAL(setSegShow_main()), this, SLOT(SegShow()));
}

void SurfaceMeshProcessing::createActions()
{
	// File
	openAction = new QAction(tr("&Open"), this);
	openAction->setIcon(QIcon(":/SurfaceMeshProcessing/Images/Open.png"));
	openAction->setShortcut(QKeySequence::Open);
	openAction->setStatusTip(tr("Open a mesh file"));
	connect(openAction, SIGNAL(triggered()), viewer, SLOT(open_mesh_query()));

	saveAction = new QAction(tr("&Save"), this);
	saveAction->setIcon(QIcon(":/SurfaceMeshProcessing/Images/Save.png"));
	saveAction->setShortcut(QKeySequence::Save);
	saveAction->setStatusTip(tr("Save the mesh to file"));
	connect(saveAction, SIGNAL(triggered()), viewer, SLOT(save_mesh_query()));

	saveAsAction = new QAction(tr("Save &As..."), this);
	saveAsAction->setStatusTip(tr("Save the mesh under a new name"));
	connect(saveAsAction, SIGNAL(triggered()), this, SLOT(saveAs()));

	exitAction = new QAction(tr("E&xit"), this);
	exitAction->setShortcut(tr("Ctrl+Q"));
	exitAction->setStatusTip(tr("Exit the application"));
	connect(exitAction, SIGNAL(triggered()), this, SLOT(close()));

	clearAction = new QAction(("Clear Mesh"), this);
	clearAction->setIcon(QIcon(":/SurfaceMeshProcessing/Images/ClearMesh.png"));
	clearAction->setStatusTip(tr("Clear the Current Mesh"));
	connect(clearAction, SIGNAL(triggered()), viewer, SLOT(clear_all_mesh()));

	texLoadAction = new QAction("Load Tex", this);
	connect(texLoadAction, SIGNAL(triggered()), viewer, SLOT(load_tex()));

	// Load Segmentation
	segLoadAction = new QAction("Load Segmentation", this);
	connect(segLoadAction, SIGNAL(triggered()), viewer, SLOT(load_seg()));

	// Output Segment Mesh
	segOutputAction = new QAction("OutPut Seg", this);
	connect(segOutputAction, SIGNAL(triggered()), viewer, SLOT(output_seg()));


	// View
	flatPointsAction = new QAction(tr("Flat&Points"), this);
	flatPointsAction->setIcon(QIcon(":/SurfaceMeshProcessing/Images/flatlines.png"));
	flatPointsAction->setStatusTip(tr("Using flatpoints showing method"));
	flatPointsAction->setCheckable(true);
	flatPointsAction->setChecked(false);
	connect(flatPointsAction, SIGNAL(triggered()), this, SLOT(FlatPointsShow()));

	segShowAction = new QAction("Segmentation", this);
	segShowAction->setIcon(QIcon(":/SurfaceMeshProcessing/Images/curvature.png"));
	segShowAction->setStatusTip(tr("Segmentation"));
	segShowAction->setCheckable(true);
	segShowAction->setChecked(false);
	connect(segShowAction, SIGNAL(triggered()), this, SLOT(SegShow()));

	texShowAction = new QAction("Segmentation Boundary", this);
	texShowAction->setIcon(QIcon(":/SurfaceMeshProcessing/Images/bbox.png"));
	texShowAction->setStatusTip(tr("Segmentation Boundary"));
	texShowAction->setCheckable(true);
	texShowAction->setChecked(false);
	connect(texShowAction, SIGNAL(triggered()), this, SLOT(TexShow()));

	// Edit
	segPickAction = new QAction(tr("&Select Face"), this);
	segPickAction->setIcon(QIcon(":/SurfaceMeshProcessing/Images/selectFaces.png"));
	segPickAction->setStatusTip(tr("Select Face of the mesh"));
	segPickAction->setCheckable(true);
	segPickAction->setChecked(false);
	connect(segPickAction, SIGNAL(triggered()), this, SLOT(SegPick()));

	boundPickAction = new QAction(tr("&Select Edge"), this);
	boundPickAction->setIcon(QIcon(":/SurfaceMeshProcessing/Images/selectEdges.png"));
	boundPickAction->setStatusTip(tr("Select Edge of the mesh"));
	boundPickAction->setCheckable(true);
	boundPickAction->setChecked(false);
	connect(boundPickAction, SIGNAL(triggered()), this, SLOT(BoundPick()));

	// Select Vertices
	vertSelectAction = new QAction("Select Vertices", this);
	vertSelectAction->setStatusTip(tr("Select Vertices"));
	vertSelectAction->setIcon(QIcon(":/SurfaceMeshProcessing/Images/selected_vert.png"));
	vertSelectAction->setCheckable(true);
	vertSelectAction->setChecked(false);
	connect(vertSelectAction, SIGNAL(triggered()), this, SLOT(VerticesSelect()));

	vertFixAction = new QAction("Fix", this);
	vertFixAction->setStatusTip(tr("Fix"));
	vertFixAction->setIcon(QIcon(":/SurfaceMeshProcessing/Images/layer_edit_locked.png"));
	connect(vertFixAction, SIGNAL(triggered()), viewer, SLOT(fix_vertices()));

	// Undo
	undoAction = new QAction("Undo", this);
	undoAction->setStatusTip(tr("Undo"));
	undoAction->setIcon(QIcon(":/SurfaceMeshProcessing/Images/edit-undo.png"));
	connect(undoAction, SIGNAL(triggered()), viewer, SLOT(undo()));


	// Drag
	vertDragAction = new QAction("Drag Vertex", this);
	vertDragAction->setStatusTip(tr("Drag Vertex"));
	vertDragAction->setIcon(QIcon(":/SurfaceMeshProcessing/Images/cursor_move.png"));
	vertDragAction->setCheckable(true);
	vertDragAction->setChecked(false);
	connect(vertDragAction, SIGNAL(triggered()), this, SLOT(VertexDrag()));


	// Draw Seg
	segDrawAction = new QAction("Draw Seg Boundary", this);
	segDrawAction->setStatusTip(tr("Draw Seg Boundary"));
	segDrawAction->setIcon(QIcon(":/SurfaceMeshProcessing/Images/selectVertices.png"));
	segDrawAction->setCheckable(true);
	segDrawAction->setChecked(false);
	connect(segDrawAction, SIGNAL(triggered()), this, SLOT(SegDraw()));

	drawUndoAction = new QAction("Undo", this);
	drawUndoAction->setStatusTip(tr("Undo Seg"));
	drawUndoAction->setIcon(QIcon(":/SurfaceMeshProcessing/Images/edit-undo.png"));
	connect(drawUndoAction, SIGNAL(triggered()), viewer, SLOT(draw_undo()));

	faceSelectAction = new QAction("Select Faces", this);
	faceSelectAction->setStatusTip(tr("Select Faces"));
	faceSelectAction->setIcon(QIcon(":/SurfaceMeshProcessing/Images/selected_face.png"));
	faceSelectAction->setCheckable(true);
	faceSelectAction->setChecked(false);
	connect(faceSelectAction, SIGNAL(triggered()), this, SLOT(FacesSelect()));

	drawFinishAction = new QAction("Close Picked Seg Vertex", this);
	drawFinishAction->setStatusTip(tr("Close Picked Seg Vertex"));
	drawFinishAction->setIcon(QIcon(":/SurfaceMeshProcessing/Images/Boundary.png"));
	connect(drawFinishAction, SIGNAL(triggered()), viewer, SLOT(finish_seg_draw()));

	finishedUndoAction = new QAction("Close Picked Seg Vertex", this);
	finishedUndoAction->setStatusTip(tr("Close Picked Seg Vertex"));
	finishedUndoAction->setIcon(QIcon(":/SurfaceMeshProcessing/Images/edit-undo.png"));
	connect(finishedUndoAction, SIGNAL(triggered()), viewer, SLOT(undo_finished()));

	segChangeAction = new QAction("Change Seg", this);
	segChangeAction->setStatusTip(tr("Change Segmentation"));
	segChangeAction->setIcon(QIcon(":/SurfaceMeshProcessing/Images/topology-edgeFlip.png"));
	segChangeAction->setCheckable(true);
	segChangeAction->setChecked(false);
	connect(segChangeAction, SIGNAL(triggered()), this, SLOT(SegChange()));

	changedUndoAction = new QAction("Undo", this);
	changedUndoAction->setStatusTip(tr("Undo Seg"));
	changedUndoAction->setIcon(QIcon(":/SurfaceMeshProcessing/Images/edit-undo.png"));
	connect(changedUndoAction, SIGNAL(triggered()), viewer, SLOT(undo_seg_change()));

	test_ = new QAction("Test", this);
	connect(test_, SIGNAL(triggered()), viewer, SLOT(test()));

	planelize = new QAction("Planar", this);
	connect(planelize, SIGNAL(triggered()), viewer, SLOT(planar()));
}

void SurfaceMeshProcessing::createMenus()
{
	fileMenu = menuBar()->addMenu(tr("&File"));
	fileMenu->addAction(openAction);
	fileMenu->addAction(saveAction);
	fileMenu->addAction(saveAsAction);
	fileMenu->addAction(clearAction);
	fileMenu->addAction(exitAction);

	viewMenu = menuBar()->addMenu(tr("&View"));
	viewMenu->addAction(flatPointsAction);
	viewMenu->addAction(segShowAction);
	viewMenu->addAction(texShowAction);
	viewMenu->addSeparator();

	editMenu = menuBar()->addMenu(tr("&Edit"));
	editMenu->addAction(segPickAction);
	editMenu->addAction(boundPickAction);
}

void SurfaceMeshProcessing::createToolBars()
{
	fileToolBar = addToolBar(tr("&File"));
	fileToolBar->addAction(openAction);
	fileToolBar->addAction(saveAction);
	fileToolBar->addAction(clearAction);
	fileToolBar->addAction(texLoadAction);
	fileToolBar->addAction(segLoadAction);
	fileToolBar->addAction(segOutputAction);

	viewToolBar = addToolBar(tr("&View"));
	viewToolBar->addAction(flatPointsAction);
	viewToolBar->addAction(segShowAction);
	viewToolBar->addAction(texShowAction);

	drawToolBar = addToolBar(tr("Draw"));
	drawToolBar->addAction(segDrawAction);
	drawToolBar->addAction(drawUndoAction);
	drawToolBar->addSeparator();
	drawToolBar->addAction(faceSelectAction);
	drawToolBar->addSeparator();
	drawToolBar->addAction(drawFinishAction);
	drawToolBar->addAction(finishedUndoAction);
	drawToolBar->addSeparator();
	drawToolBar->addAction(segChangeAction);
	drawToolBar->addAction(changedUndoAction);


	editToolBar = addToolBar(tr("Edit"));
	editToolBar->addAction(segPickAction);
	editToolBar->addAction(boundPickAction);
	editToolBar->addAction(vertSelectAction);
	editToolBar->addAction(vertFixAction);
	editToolBar->addAction(vertDragAction);
	editToolBar->addAction(undoAction);
	
	editToolBar->addAction(test_);
	editToolBar->addAction(planelize);
}

void SurfaceMeshProcessing::createStatusBar()
{
	statusLabel = new QLabel(tr("No mesh"));
	statusLabel->setAlignment(Qt::AlignHCenter);
	connect(viewer, SIGNAL(haveLoadMesh(QString)), statusLabel, SLOT(setText(QString)));
	statusBar()->addWidget(statusLabel);
}


bool SurfaceMeshProcessing::save()
{
	return true;
}

bool SurfaceMeshProcessing::saveAs()
{
	return true;
}

void SurfaceMeshProcessing::FlatPointsShow()
{
	viewer->setDrawMode(InteractiveViewerWidget::FLAT_POINTS);
	viewer->setEditMode(InteractiveViewerWidget::TRANS);
}

void SurfaceMeshProcessing::SegShow()
{
	viewer->setDrawMode(InteractiveViewerWidget::SEGMENT);
	viewer->setEditMode(InteractiveViewerWidget::TRANS);
}

void SurfaceMeshProcessing::TexShow()
{
	viewer->setDrawMode(InteractiveViewerWidget::TEXTURE);
	viewer->setEditMode(InteractiveViewerWidget::TRANS);
}

void SurfaceMeshProcessing::setAllViewActionChecked(bool b)
{
	flatPointsAction->setChecked(b);
	segShowAction->setChecked(b);
	texShowAction->setChecked(b);
}


void SurfaceMeshProcessing::BoundPick()
{
	setAllViewActionChecked(false);
	setAllEditActionChecked(false);
	boundPickAction->setChecked(true);
	viewer->setEditMode(InteractiveViewerWidget::BOUND_PICK);
}
void SurfaceMeshProcessing::SegPick()
{
	setAllViewActionChecked(false);
	setAllEditActionChecked(false);
	segPickAction->setChecked(true);
	viewer->setEditMode(InteractiveViewerWidget::SEG_PICK);
}

void SurfaceMeshProcessing::SegChange() {
	setAllViewActionChecked(false);
	setAllEditActionChecked(false);
	segChangeAction->setChecked(true);
	viewer->setEditMode(InteractiveViewerWidget::SEG_CHANGE);
}

void SurfaceMeshProcessing::SegDraw(){
	setAllViewActionChecked(false);
	setAllEditActionChecked(false);
	segDrawAction->setChecked(true);
	viewer->setEditMode(InteractiveViewerWidget::SEG_DRAW);
}

void SurfaceMeshProcessing::VerticesSelect()
{
	setAllViewActionChecked(false);
	setAllEditActionChecked(false);
	vertSelectAction->setChecked(true);
	viewer->setEditMode(InteractiveViewerWidget::VERT_SELECT);
}

void SurfaceMeshProcessing::FacesSelect()
{
	setAllViewActionChecked(false);
	setAllEditActionChecked(false);
	faceSelectAction->setChecked(true);
	viewer->setEditMode(InteractiveViewerWidget::FACE_SELECT);
}

void SurfaceMeshProcessing::VertexDrag()
{
	setAllViewActionChecked(false);
	setAllEditActionChecked(false);
	vertDragAction->setChecked(true);
	viewer->setEditMode(InteractiveViewerWidget::VERT_DRAG);
}

void SurfaceMeshProcessing::setAllEditActionChecked(bool b)
{
	segPickAction->setChecked(b);
	boundPickAction->setChecked(b);
	segChangeAction->setChecked(b);
	segDrawAction->setChecked(b);
	vertSelectAction->setChecked(b);
	faceSelectAction->setChecked(b);
	vertDragAction->setChecked(b);
}

void SurfaceMeshProcessing::setEditMode_slot(int mm)
{
	if (mm == InteractiveViewerWidget::TRANS)
	{
		setAllViewActionChecked(false);
		setAllEditActionChecked(false);
		switch (latestDrawMode)
		{
		case InteractiveViewerWidget::FLAT_POINTS:
			flatPointsAction->setChecked(true);
			break;
		case InteractiveViewerWidget::SEGMENT:
			segShowAction->setChecked(true);
			break;
		case InteractiveViewerWidget::TEXTURE:
			texShowAction->setChecked(true);
			break;
		}
	}
	else
	{
		setAllViewActionChecked(false);
	}
}

void SurfaceMeshProcessing::setDrawMode_slot(int dm)
{
	setAllEditActionChecked(false);
	latestDrawMode = dm;
}
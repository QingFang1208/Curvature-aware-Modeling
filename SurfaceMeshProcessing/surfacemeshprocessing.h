#ifndef SURFACEMESHPROCESSING_H
#define SURFACEMESHPROCESSING_H

#include <QMainWindow>
#include <QtGui>
#include <QtWidgets>

class QAction;
class MainViewerWidget;

class SurfaceMeshProcessing : public QMainWindow
{
	Q_OBJECT
public:
	SurfaceMeshProcessing(QWidget *parent = 0);
	~SurfaceMeshProcessing();

private slots:
	bool save();
	bool saveAs();

	void setEditMode_slot(int mm);
	void setDrawMode_slot(int dm);

private:
	void initWindow();
	void createActions();
	void createMenus();
	void createToolBars();
	void createStatusBar();

private slots:
	void FlatPointsShow();
	void SegShow();
	void TexShow();

	void BoundPick();
	void SegPick();
	void SegChange();
	void SegDraw();
	void VerticesSelect();
	void FacesSelect();
	void VertexDrag();

private:
	// File Actions.
	QAction* openAction;
	QAction* saveAction;
	QAction* saveAsAction;
	QAction* exitAction;
	QAction* clearAction;

	QAction* texLoadAction;
	QAction* segLoadAction;
	QAction* segOutputAction;

	// View Actions.
	QAction* flatPointsAction;
	QAction* segShowAction;
	QAction* texShowAction;

	// Menus.
	QMenu* fileMenu;
	QMenu* viewMenu;
	QMenu* editMenu;

	// ToolBars
	QToolBar* fileToolBar;
	QToolBar* viewToolBar;
	QToolBar* editToolBar;
	QToolBar* drawToolBar;

	// Label
	QLabel* statusLabel;

	// Edit Actions
	QAction* segPickAction;
	QAction* boundPickAction;
	QAction* vertSelectAction;
	QAction* vertFixAction;
	QAction* vertDragAction;
	QAction* undoAction;

	// Draw
	QAction* segDrawAction;
	QAction* drawUndoAction;
	QAction* faceSelectAction;
	QAction* drawFinishAction;
	QAction* finishedUndoAction;
	QAction* segChangeAction;
	QAction* changedUndoAction;

	// test
	QAction* test_;
	QAction* planelize;


private:
	void setAllEditActionChecked(bool b);
	void setAllViewActionChecked(bool b);

private:
	MainViewerWidget* viewer;
	int latestDrawMode;

};

#endif // SURFACEMESHPROCESSING_H

#pragma once

#include <QWidget>
#include <QtGui>
#include <QtWidgets>

class MeshParamDialog : public QWidget
{
	Q_OBJECT

public:
	MeshParamDialog(QWidget *parent = 0);
	~MeshParamDialog(void);

private:
	void CreateTabWidget(void);
	void CreateLayout(void);
	void CreatePrintWidget(QFont& ft, QPalette& pe);
	void CreateMaterialWidget(QFont& ft, QPalette& pe);
	void CreateEditWidget(QFont& ft, QPalette& pe);

signals:
	// Print
	void PrintInfoSignal();

	// Init Draw
	void InitSphereSignal();
	void InitCylinderSignal();
	void InitDrawSignal();

	// Customize
	void customizeSignal(int);
	void sliderSignal(int);
	void runAlgorithm();
	void refreshMesh();
	void assembleMesh();

	// Material
	void ShineSignal(int);
	void SpecularSignal(int);

private:
	QTabWidget *twParam;
	QWidget *wParam;
	QScrollArea *saParam;

	// Print
	QLabel *PrintLab;
	QPushButton *PrintInfoButton;

	// Material
	QLabel *MaterialLab;
	QLabel *ShineLab;
	QLabel *SpecularLab;
	QSlider *ShineSlider;
	QSlider *SpecularSlider;

	// Draw
	QLabel *editLab;
	QButtonGroup *customize;
	QPushButton *spherizeOpen;
	QPushButton *cylinderize;
	QPushButton *planePolygon;
	QPushButton *planeFree;
	QPushButton *SpherePolygon;
	QPushButton *SmoothControl;
	QPushButton *FeatureControl;
	QPushButton *DragVertices;
	QPushButton *algorithmRun;
	QPushButton *meshRefresh;
	QPushButton *meshAssemble;
	QSlider *paraSlider;
	QLabel *paraSliderLable0;
	QLabel *paraSliderLable1;

};

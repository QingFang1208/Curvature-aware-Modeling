#include "MeshParamDialog.h"

MeshParamDialog::MeshParamDialog(QWidget *parent)
	: QWidget(parent)
{
	CreateTabWidget();
	CreateLayout();
}

MeshParamDialog::~MeshParamDialog()
{
}

void MeshParamDialog::CreateTabWidget(void)
{
	// Font
	QFont ft;
	ft.setPointSize(10);
	ft.setBold(true);

	// Palette
	QPalette pe;
	pe.setColor(QPalette::WindowText, Qt::darkBlue);

	// Print
	CreatePrintWidget(ft, pe);

	// Material
	CreateMaterialWidget(ft, pe);

	// Draw
	CreateEditWidget(ft, pe);

	// Layout
	QGridLayout *layout = new QGridLayout();

	layout->addWidget(PrintLab, 0, 0, Qt::AlignLeft);
	layout->addWidget(PrintInfoButton, 0, 3);

	layout->addWidget(MaterialLab, 1, 0);
	layout->addWidget(ShineLab, 2, 0);
	layout->addWidget(ShineSlider, 2, 1, 1, 4);
	layout->addWidget(SpecularLab, 3, 0);
	layout->addWidget(SpecularSlider, 3, 1, 1, 4);

	layout->addWidget(editLab, 4, 0);
	layout->addWidget(spherizeOpen, 5, 0);
	layout->addWidget(cylinderize, 5, 1);
	//layout->addWidget(planePolygon, 5, 2);
	//layout->addWidget(planeFree, 5, 3);
	//layout->addWidget(SpherePolygon, 5, 4);
	layout->addWidget(SmoothControl, 5, 2);
	layout->addWidget(FeatureControl, 6, 0);
	layout->addWidget(DragVertices, 6, 1);
	layout->addWidget(algorithmRun, 7, 0);
	layout->addWidget(meshAssemble, 7, 1);
	layout->addWidget(meshRefresh, 7, 2);
	layout->addWidget(paraSlider, 8, 0, 2, 5);
	layout->addWidget(paraSliderLable0, 9, 0, Qt::AlignLeft);
	layout->addWidget(paraSliderLable1, 9, 4, Qt::AlignRight);

	wParam = new QWidget();
	wParam->setLayout(layout);
	saParam = new QScrollArea();
	saParam->setFocusPolicy(Qt::NoFocus);
	saParam->setFrameStyle(QFrame::NoFrame);
	saParam->setWidget(wParam);
	saParam->setWidgetResizable(true);
}

void MeshParamDialog::CreateLayout(void)
{
	twParam = new QTabWidget();
	twParam->addTab(saParam, "Tab");
	QGridLayout *layout = new QGridLayout();
	layout->addWidget(twParam, 0, 0, 1, 1);
	this->setLayout(layout);
}

void MeshParamDialog::CreatePrintWidget(QFont& ft, QPalette& pe)
{
	// Label
	PrintLab = new QLabel(tr("Print Info: "));
	PrintLab->setFont(ft);
	PrintLab->setPalette(pe);

	// Pushbutton
	PrintInfoButton = new QPushButton(tr("Print"));
	PrintInfoButton->setFixedSize(60, 20);
	connect(PrintInfoButton, SIGNAL(clicked()), SIGNAL(PrintInfoSignal()));
}

void MeshParamDialog::CreateMaterialWidget(QFont & ft, QPalette & pe)
{
	MaterialLab = new QLabel(tr("Material:"));
	MaterialLab->setFont(ft);
	MaterialLab->setPalette(pe);

	ShineLab = new QLabel(tr("Shine:"));

	ShineSlider = new QSlider(this);
	ShineSlider->setOrientation(Qt::Horizontal);
	ShineSlider->setTickPosition(QSlider::TicksAbove);
	ShineSlider->setMinimum(0);
	ShineSlider->setMaximum(128);
	ShineSlider->setValue(20);
	ShineSlider->setSingleStep(10);
	connect(ShineSlider, SIGNAL(valueChanged(int)), SIGNAL(ShineSignal(int)));

	SpecularLab = new QLabel(tr("Specular:"));
	SpecularSlider = new QSlider(this);
	SpecularSlider->setOrientation(Qt::Horizontal);
	SpecularSlider->setTickPosition(QSlider::TicksAbove);
	SpecularSlider->setMinimum(0);
	SpecularSlider->setMaximum(100);
	SpecularSlider->setValue(0);
	SpecularSlider->setSingleStep(10);
	connect(SpecularSlider, SIGNAL(valueChanged(int)), SIGNAL(SpecularSignal(int)));

}

void MeshParamDialog::CreateEditWidget(QFont & ft, QPalette & pe)
{
	editLab = new QLabel(tr("Editing:"));
	editLab->setFont(ft);
	editLab->setPalette(pe);

	spherizeOpen = new QPushButton();
	spherizeOpen->setIcon(QIcon(":/SurfaceMeshProcessing/Images/sphere1.png"));
	spherizeOpen->setIconSize(QSize(40, 40));
	spherizeOpen->setCheckable(true);
	spherizeOpen->setChecked(true);

	cylinderize = new QPushButton();
	cylinderize->setIcon(QIcon(":/SurfaceMeshProcessing/Images/cylinder1.png"));
	cylinderize->setIconSize(QSize(40, 40));
	cylinderize->setCheckable(true);

	planePolygon = new QPushButton();
	planePolygon->setIcon(QIcon(":/SurfaceMeshProcessing/Images/planePolygon1.png"));
	planePolygon->setIconSize(QSize(40, 40));
	planePolygon->setCheckable(true);

	planeFree = new QPushButton();
	planeFree->setIcon(QIcon(":/SurfaceMeshProcessing/Images/planeFree1.png"));
	planeFree->setIconSize(QSize(40, 40));
	planeFree->setCheckable(true);

	SpherePolygon = new QPushButton();
	SpherePolygon->setIcon(QIcon(":/SurfaceMeshProcessing/Images/spherePolygon1.png"));
	SpherePolygon->setIconSize(QSize(40, 40));
	SpherePolygon->setCheckable(true);

	SmoothControl = new QPushButton();
	SmoothControl->setIcon(QIcon(":/SurfaceMeshProcessing/Images/smoothControl1.png"));
	SmoothControl->setIconSize(QSize(40, 40));
	SmoothControl->setCheckable(true);

	FeatureControl = new QPushButton();
	FeatureControl->setIcon(QIcon(":/SurfaceMeshProcessing/Images/featureControl1.png"));
	FeatureControl->setIconSize(QSize(40, 40));
	FeatureControl->setCheckable(true);

	DragVertices = new QPushButton();
	DragVertices->setIcon(QIcon(":/SurfaceMeshProcessing/Images/dragVertices1.png"));
	DragVertices->setIconSize(QSize(40, 40));
	DragVertices->setCheckable(true);

	customize = new QButtonGroup();
	customize->addButton(spherizeOpen, 0);
	customize->addButton(cylinderize, 1);
	customize->addButton(planePolygon, 2);
	customize->addButton(planeFree, 3);
	customize->addButton(SpherePolygon, 4);
	customize->addButton(SmoothControl, 5);
	customize->addButton(FeatureControl, 6);
	customize->addButton(DragVertices, 7);
	customize->setExclusive(true);
	connect(customize, SIGNAL(buttonClicked(int)), this, SIGNAL(customizeSignal(int)));

	algorithmRun = new QPushButton(tr("run"));
	connect(algorithmRun, SIGNAL(clicked()), SIGNAL(runAlgorithm()));

	meshRefresh = new QPushButton(tr("refresh"));
	connect(meshRefresh, SIGNAL(clicked()), SIGNAL(refreshMesh()));

	meshAssemble = new QPushButton(tr("assemble"));
	connect(meshAssemble, SIGNAL(clicked()), SIGNAL(assembleMesh()));

	paraSlider = new QSlider(this);
	paraSlider->setOrientation(Qt::Horizontal);
	paraSlider->setTickPosition(QSlider::TicksAbove);
	paraSlider->setMinimum(0);
	paraSlider->setMaximum(40);
	paraSlider->setValue(20);
	paraSlider->setSingleStep(4);
	connect(paraSlider, SIGNAL(valueChanged(int)), SIGNAL(sliderSignal(int)));

	paraSliderLable0 = new QLabel(tr("0"));
	paraSliderLable1 = new QLabel(tr("40"));
}

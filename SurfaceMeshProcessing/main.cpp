#include "surfacemeshprocessing.h"
#include <QtWidgets/QApplication>

int main(int argc, char *argv[])
{
	QApplication app(argc, argv);

	SurfaceMeshProcessing mainWin;
	mainWin.showMaximized();

	return app.exec();
}

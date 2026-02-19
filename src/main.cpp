#include <QApplication>
#include <QMainWindow>
#include <QSurfaceFormat>

#include "view/gl_widget.h"
#include "view/interface.h"

int main(int argc, char* argv[]) {
	QSurfaceFormat fmt;
	fmt.setVersion(3, 3);
	fmt.setProfile(QSurfaceFormat::CoreProfile);
	fmt.setDepthBufferSize(24);
	QSurfaceFormat::setDefaultFormat(fmt);

	QApplication app(argc, argv);

	QMainWindow w;
	auto* gl = new s21::GlWidget(&w);

	auto* ui = new s21::Interface(gl);

	// w.setCentralWidget(gl);
	w.setCentralWidget(ui);
	w.resize(900, 700);
	w.show();

	return app.exec();
}

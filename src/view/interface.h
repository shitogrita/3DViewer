#pragma once

#include <QLabel>
#include <QWidget>

class QButtonGroup;
class QFrame;
class QStackedWidget;

namespace s21 {

	class GlWidget;

	class Interface final : public QWidget {
		Q_OBJECT

	   public:
		explicit Interface(GlWidget* viewport, QWidget* parent = nullptr);
		~Interface() override;

	protected:
		bool eventFilter(QObject* obj, QEvent* ev) override;
		void resizeEvent(QResizeEvent* e) override;

	private:
		void BuildUi_();
		void ApplyStyle_();
		void MakeControlsNoFocus_(QWidget* root);

		void UpdatePanelHeight_();

		void OnTabClicked(int id);
		void OnClosePanel();

		static QByteArray ImageToBytes_(const QImage& img,
								 const char* format,
								 int quality);

		void SaveRenderedImage_();

		QWidget* MakePageProjection_();
		QWidget* MakePageEdges_();
		QWidget* MakePageVertices_();
		QWidget* MakePageBackground_();
		QWidget* MakePageFile_();

		GlWidget* viewport_ = nullptr;
		QWidget* topBar_ = nullptr;
		QButtonGroup* tabs_ = nullptr;
		QFrame* panel_ = nullptr;
		QFrame* panelBox_ = nullptr;
		QStackedWidget* pages_ = nullptr;

		// interface.h (внутри class Interface)
	private:
		QWidget* infoBar_ = nullptr;
		QLabel*  fileValue_ = nullptr;
		QLabel*  vtxValue_  = nullptr;
		QLabel*  edgeValue_ = nullptr;

		void SetInfoBarVisible_(bool on);
	};
}

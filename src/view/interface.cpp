#include "view/interface.h"
#include "view/gl_widget.h"

#include <QApplication>
#include <QButtonGroup>
#include <QCheckBox>
#include <QColorDialog>
#include <QComboBox>
#include <QCoreApplication>
#include <QDoubleSpinBox>
#include <QEvent>
#include <QFrame>
#include <QHBoxLayout>
#include <QLabel>
#include <QPushButton>
#include <QResizeEvent>
#include <QStackedWidget>
#include <QToolButton>
#include <QVBoxLayout>
#include <QOpenGLWidget>
#include <QFileDialog>
#include <QFileInfo>
#include <QBuffer>
#include <QMessageBox>

namespace s21 {
	QToolButton *MakeTabButton(const QString &text) {
		auto *b = new QToolButton();
		b->setText(text);
		b->setCheckable(true);
		b->setAutoRaise(true);
		b->setToolButtonStyle(Qt::ToolButtonTextOnly);
		b->setSizePolicy(QSizePolicy::Expanding, QSizePolicy::Fixed);
		b->setFixedHeight(36);
		b->setFocusPolicy(Qt::NoFocus);
		return b;
	}

	bool IsTextInputWidget(QWidget *w) {
		if (!w) return false;
		return w->inherits("QLineEdit") || w->inherits("QTextEdit") || w->inherits("QPlainTextEdit") ||
				w->inherits("QAbstractSpinBox") || w->inherits("QComboBox");
	}

	Interface::Interface(GlWidget *viewport, QWidget *parent)
		: QWidget(parent), viewport_(viewport) {
		BuildUi_();
		ApplyStyle_();

		if (viewport_) {
			if (fileValue_) fileValue_->setText(viewport_->GetModelFileName());
			if (vtxValue_)  vtxValue_->setText(QString::number(viewport_->GetVertexCount()));

			connect(viewport_, &s21::GlWidget::ModelInfoChanged, this,
					[this](const QString& name, int vtx, int edges) {
					  if (fileValue_) fileValue_->setText(name.isEmpty() ? "—" : name);
					  if (vtxValue_)  vtxValue_->setText(QString::number(vtx));
					  if (edgeValue_) edgeValue_->setText(QString::number(edges));
					});
		}
		setFocusPolicy(Qt::StrongFocus);
		if (viewport_) {
			setFocusProxy(viewport_);
			viewport_->setFocusPolicy(Qt::StrongFocus);
			viewport_->setFocus();
		}

		MakeControlsNoFocus_(this);
		qApp->installEventFilter(this);

		UpdatePanelHeight_();
	}

	void Interface::SetInfoBarVisible_(bool on) {
		if (infoBar_) infoBar_->setVisible(on);
	}

	void Interface::OnTabClicked(int id) {
		panel_->setVisible(true);
		pages_->setCurrentIndex(id);
		UpdatePanelHeight_();
		SetInfoBarVisible_(false);
		if (viewport_) viewport_->setFocus();
	}

	void Interface::OnClosePanel() {
		panel_->setVisible(false);
		if (auto* checked = tabs_->checkedButton()) checked->setChecked(false);
		SetInfoBarVisible_(true);
		if (viewport_) viewport_->setFocus();
	}

	Interface::~Interface() {
		qApp->removeEventFilter(this);
	}

	void Interface::resizeEvent(QResizeEvent *e) {
		QWidget::resizeEvent(e);
		UpdatePanelHeight_();
	}

	void Interface::UpdatePanelHeight_() {
		if (!panel_ || !topBar_) return;

		if (!panel_->isVisible()) return;

		const int topH = topBar_->height();
		const int avail = std::max(0, height() - topH);
		const int half = std::max(180, avail / 2);
		panel_->setFixedHeight(half);

		if (panelBox_) panelBox_->setSizePolicy(QSizePolicy::Expanding, QSizePolicy::Expanding);
		if (pages_) pages_->setSizePolicy(QSizePolicy::Expanding, QSizePolicy::Expanding);
	}

	bool Interface::eventFilter(QObject *obj, QEvent *ev) {
		if (!viewport_) return QWidget::eventFilter(obj, ev);

		if (ev->type() != QEvent::KeyPress && ev->type() != QEvent::KeyRelease) {
			return QWidget::eventFilter(obj, ev);
		}

		if (QApplication::activeModalWidget() != nullptr) {
			return QWidget::eventFilter(obj, ev);
		}

		if (obj == viewport_) {
			return QWidget::eventFilter(obj, ev);
		}

		QWidget *fw = QApplication::focusWidget();
		if (IsTextInputWidget(fw)) {
			return QWidget::eventFilter(obj, ev);
		}

		thread_local bool forwarding = false;
		if (forwarding) return QWidget::eventFilter(obj, ev);

		forwarding = true;
		QCoreApplication::sendEvent(viewport_, ev);
		forwarding = false;
		return true;
	}

	void Interface::BuildUi_() {
		auto *root = new QVBoxLayout(this);
		root->setContentsMargins(0, 0, 0, 0);
		root->setSpacing(0);

		topBar_ = new QWidget(this);
		topBar_->setObjectName("TopBar");
		topBar_->setFixedHeight(36);
		topBar_->setSizePolicy(QSizePolicy::Expanding, QSizePolicy::Fixed);

		auto *topL = new QHBoxLayout(topBar_);
		topL->setContentsMargins(0, 0, 0, 0);
		topL->setSpacing(0);

		tabs_ = new QButtonGroup(this);
		tabs_->setExclusive(true);

		const struct {
			const char *text;
			int id;
		} items[] = {
					{"Проекция", 0},
					{"Рёбра", 1},
					{"Вершины", 2},
					{"Покраска", 3},
					{"Файл", 4},
				};

		for (const auto &it: items) {
			auto *b = MakeTabButton(QString::fromUtf8(it.text));
			tabs_->addButton(b, it.id);
			topL->addWidget(b, 1);
		}

		connect(tabs_, &QButtonGroup::idClicked, this, &Interface::OnTabClicked);

		infoBar_ = new QFrame(this);
		infoBar_->setObjectName("InfoBar");
		infoBar_->setFixedHeight(28);
		infoBar_->setSizePolicy(QSizePolicy::Expanding, QSizePolicy::Fixed);

		auto* infoL = new QHBoxLayout(infoBar_);
		infoL->setContentsMargins(12, 0, 12, 0);
		infoL->setSpacing(10);

		auto* fileLab = new QLabel("Файл:", infoBar_);
		fileLab->setMinimumWidth(50);

		fileValue_ = new QLabel("—", infoBar_);
		fileValue_->setSizePolicy(QSizePolicy::Expanding, QSizePolicy::Preferred);

		auto* vtxLab = new QLabel("Вершин:", infoBar_);
		vtxLab->setMinimumWidth(70);

		vtxValue_ = new QLabel("0", infoBar_);
		vtxValue_->setMinimumWidth(80);
		vtxValue_->setAlignment(Qt::AlignRight | Qt::AlignVCenter);


		auto* edgeLab = new QLabel("Ребер: ", infoBar_);
		edgeLab->setMinimumWidth(80);

		edgeValue_ = new QLabel("0", infoBar_);
		edgeValue_->setMinimumWidth(80);
		edgeValue_->setAlignment(Qt::AlignRight | Qt::AlignVCenter);
		infoL->addWidget(fileLab, 0);
		infoL->addWidget(fileValue_, 1);
		infoL->addSpacing(10);
		infoL->addWidget(vtxLab, 0);
		infoL->addWidget(vtxValue_, 0);
		infoL->addSpacing(10);
		infoL->addWidget(edgeLab, 0);
		infoL->addWidget(edgeValue_, 0);

		panel_ = new QFrame(this);
		panel_->setObjectName("Panel");
		panel_->setVisible(false);
		panel_->setFrameShape(QFrame::NoFrame);
		panel_->setSizePolicy(QSizePolicy::Expanding, QSizePolicy::Fixed);

		auto *panelL = new QVBoxLayout(panel_);
		panelL->setContentsMargins(10, 10, 10, 10);
		panelL->setSpacing(0);

		panelBox_ = new QFrame(panel_);
		panelBox_->setObjectName("PanelBox");
		panelBox_->setFrameShape(QFrame::NoFrame);
		panelBox_->setSizePolicy(QSizePolicy::Expanding, QSizePolicy::Expanding);

		auto *boxL = new QVBoxLayout(panelBox_);
		boxL->setContentsMargins(12, 10, 12, 12);
		boxL->setSpacing(10);

		auto *header = new QWidget(panelBox_);
		header->setSizePolicy(QSizePolicy::Expanding, QSizePolicy::Fixed);

		auto *headerL = new QHBoxLayout(header);
		headerL->setContentsMargins(0, 0, 0, 0);
		headerL->setSpacing(8);

		auto *title = new QLabel("Настройки", header);
		title->setObjectName("PanelTitle");
		title->setSizePolicy(QSizePolicy::Expanding, QSizePolicy::Preferred);

		auto *closeBtn = new QPushButton("×", header);
		closeBtn->setObjectName("CloseBtn");
		closeBtn->setFixedSize(28, 28);
		closeBtn->setFocusPolicy(Qt::NoFocus);
		connect(closeBtn, &QPushButton::clicked, this, &Interface::OnClosePanel);

		headerL->addWidget(title, 1);
		headerL->addWidget(closeBtn, 0);

		pages_ = new QStackedWidget(panelBox_);
		pages_->setSizePolicy(QSizePolicy::Expanding, QSizePolicy::Expanding);

		QWidget *p0 = MakePageProjection_();
		QWidget *p1 = MakePageEdges_();
		QWidget *p2 = MakePageVertices_();
		QWidget *p3 = MakePageBackground_();
		QWidget *p4 = MakePageFile_();

		for (QWidget *p: {p0, p1, p2, p3, p4}) {
			p->setSizePolicy(QSizePolicy::Expanding, QSizePolicy::Expanding);
		}

		pages_->addWidget(p0);
		pages_->addWidget(p1);
		pages_->addWidget(p2);
		pages_->addWidget(p3);
		pages_->addWidget(p4);

		boxL->addWidget(header, 0);
		boxL->addWidget(pages_, 1);

		panelL->addWidget(panelBox_, 1);

		if (viewport_) {
			viewport_->setParent(this);
			viewport_->setSizePolicy(QSizePolicy::Expanding, QSizePolicy::Expanding);
		}

		root->addWidget(topBar_, 0);
		root->addWidget(infoBar_, 0);
		root->addWidget(panel_, 0);
		if (viewport_) root->addWidget(viewport_, 1);
	}

	void Interface::ApplyStyle_() {
		setStyleSheet(R"(
      QWidget { background: #d6d6d6; color: #000000; }
      QWidget#TopBar { background: #c7c7c7; }

      QToolButton {
        border: 0px;
        margin: 0px;
        padding: 0px 10px;
        background: #c7c7c7;
        color: #000000;
        font-weight: 600;
      }
      QToolButton:hover { background: #bbbbbb; }
      QToolButton:checked { background: #b0b0b0; font-weight: 700; }

      QFrame#Panel { background: transparent; }
      QFrame#PanelBox {
        background: #e1e1e1;
        border: 1px solid #7a7a7a;
        border-radius: 8px;
      }

		QFrame#InfoBar {
		  background: #dcdcdc;
		  border-top: 1px solid #b0b0b0;
		  border-bottom: 1px solid #b0b0b0;
		}

      QLabel#PanelTitle { background: transparent; color: #000000; font-weight: 800; }
      QLabel { background: transparent; color: #000000; }

      QPushButton {
        background: #e3e3e3;
        color: #000000;
        border: 1px solid #8f8f8f;
        padding: 6px 10px;
        min-height: 34px;
      }
      QPushButton:hover { background: #d9d9d9; }
      QPushButton:pressed { background: #cfcfcf; }

      QPushButton#CloseBtn {
        min-height: 28px; min-width: 28px;
        max-height: 28px; max-width: 28px;
        padding: 0px; font-weight: 900;
      }

      QCheckBox { background: transparent; color: #000000; spacing: 8px; }

      QComboBox, QDoubleSpinBox, QSpinBox {
        background: #e3e3e3;
        color: #000000;
        border: 1px solid #8f8f8f;
        padding: 4px 8px;
        min-height: 32px;
      }

      QComboBox QAbstractItemView {
        background: #e3e3e3;
        color: #000000;
        selection-background-color: #bdbdbd;
        selection-color: #000000;
        outline: 0;
      }
    )");
	}

	void Interface::MakeControlsNoFocus_(QWidget *root) {
		const auto children = root->findChildren<QWidget *>();
		for (auto *w: children) {
			if (IsTextInputWidget(w)) continue;
			w->setFocusPolicy(Qt::NoFocus);
		}
	}

	QWidget *Interface::MakePageProjection_() {
		auto *w = new QWidget();
		auto *l = new QVBoxLayout(w);
		l->setContentsMargins(0, 0, 0, 0);
		l->setSpacing(10);

		auto *row = new QWidget(w);
		auto *rowL = new QHBoxLayout(row);
		rowL->setContentsMargins(0, 0, 0, 0);
		rowL->setSpacing(8);

		auto *lab = new QLabel("Тип проекции:", row);
		lab->setMinimumWidth(130);

		auto *bOrtho = new QPushButton("Параллельная", row);
		auto *bPersp = new QPushButton("Центральная", row);
		bOrtho->setFocusPolicy(Qt::NoFocus);
		bPersp->setFocusPolicy(Qt::NoFocus);
		bOrtho->setCheckable(true);
		bPersp->setCheckable(true);

		auto *grp = new QButtonGroup(row);
		grp->setExclusive(true);
		grp->addButton(bOrtho, 0);
		grp->addButton(bPersp, 1);
		bOrtho->setChecked(true);

		connect(grp, &QButtonGroup::idClicked, this, [this](int pid) {
			if (viewport_) viewport_->SetProjectionType(pid);
			if (viewport_) viewport_->setFocus();
		});

		rowL->addWidget(lab, 0);
		rowL->addWidget(bOrtho, 1);
		rowL->addWidget(bPersp, 1);

		l->addWidget(row);

		l->addWidget(new QLabel("Точные повороты (градусы):", w));

		auto makeAngle = []() {
			auto *s = new QDoubleSpinBox();
			s->setRange(-360000.0, 360000.0);
			s->setDecimals(1);
			s->setSingleStep(1.0);
			s->setKeyboardTracking(false);
			s->setFocusPolicy(Qt::ClickFocus);
			return s;
		};

		auto *ax = makeAngle();
		auto *ay = makeAngle();
		auto *az = makeAngle();

		auto addAngleRow = [&](const QString &name, QDoubleSpinBox *s) {
			auto *r = new QWidget(w);
			auto *rl = new QHBoxLayout(r);
			rl->setContentsMargins(0, 0, 0, 0);
			rl->setSpacing(8);
			auto *la = new QLabel(name, r);
			la->setMinimumWidth(130);
			rl->addWidget(la, 0);
			rl->addWidget(s, 1);
			l->addWidget(r);
		};

		addAngleRow("Поворот X:", ax);
		addAngleRow("Поворот Y:", ay);
		addAngleRow("Поворот Z:", az);

		if (viewport_) {
			const auto r = viewport_->GetRotationDegrees();
			ax->setValue(r.x);
			ay->setValue(r.y);
			az->setValue(r.z);
		}

		if (viewport_) {
			connect(viewport_, &s21::GlWidget::RotationChanged, this,
					[ax, ay, az](double x, double y, double z) {
						const QSignalBlocker bx(ax);
						const QSignalBlocker by(ay);
						const QSignalBlocker bz(az);

						ax->setValue(x);
						ay->setValue(y);
						az->setValue(z);
					});
		}

		auto apply = [this, ax, ay, az]() {
			if (!viewport_) return;
			viewport_->SetRotationDegrees(ax->value(), ay->value(), az->value());
		};

		connect(ax, qOverload<double>(&QDoubleSpinBox::valueChanged), this, [apply](double) { apply(); });
		connect(ay, qOverload<double>(&QDoubleSpinBox::valueChanged), this, [apply](double) { apply(); });
		connect(az, qOverload<double>(&QDoubleSpinBox::valueChanged), this, [apply](double) { apply(); });

		connect(ax, &QDoubleSpinBox::editingFinished, this, [this, apply]() {
			apply();
			if (viewport_) viewport_->setFocus();
		});
		connect(ay, &QDoubleSpinBox::editingFinished, this, [this, apply]() {
			apply();
			if (viewport_) viewport_->setFocus();
		});
		connect(az, &QDoubleSpinBox::editingFinished, this, [this, apply]() {
			apply();
			if (viewport_) viewport_->setFocus();
		});

		auto *reset = new QPushButton("Сброс поворота", w);
		reset->setFocusPolicy(Qt::NoFocus);
		connect(reset, &QPushButton::clicked, this, [this, ax, ay, az] {
			if (!viewport_) return;
			viewport_->ResetRotation();
			ax->setValue(0.0);
			ay->setValue(0.0);
			az->setValue(0.0);
			viewport_->setFocus();
		});
		l->addWidget(reset);


		l->addStretch(1);
		return w;
	}

	QWidget *Interface::MakePageEdges_() {
		auto *w = new QWidget();
		auto *l = new QVBoxLayout(w);
		l->setContentsMargins(0, 0, 0, 0);
		l->setSpacing(10);

		auto *dashed = new QCheckBox("Пунктир");
		connect(dashed, &QCheckBox::toggled, this, [this](bool on) {
			if (viewport_) viewport_->SetEdgesDashed(on);
			if (viewport_) viewport_->setFocus();
		});

		auto *width = new QDoubleSpinBox();
		width->setRange(1.0, 10.0);
		width->setSingleStep(1.0);
		width->setValue(1.0);
		connect(width, qOverload<double>(&QDoubleSpinBox::valueChanged), this, [this](double v) {
			if (viewport_) viewport_->SetEdgeWidth(static_cast<float>(v));
		});

		auto *dashPeriod = new QDoubleSpinBox();
		dashPeriod->setRange(0.005, 1.0);
		dashPeriod->setSingleStep(0.005);
		dashPeriod->setValue(0.05);
		connect(dashPeriod, qOverload<double>(&QDoubleSpinBox::valueChanged), this, [this](double v) {
			if (viewport_) viewport_->SetDashPeriod(static_cast<float>(v));
		});

		auto *dashFill = new QDoubleSpinBox();
		dashFill->setRange(0.05, 0.95);
		dashFill->setSingleStep(0.05);
		dashFill->setValue(0.50);
		connect(dashFill, qOverload<double>(&QDoubleSpinBox::valueChanged), this, [this](double v) {
			if (viewport_) viewport_->SetDashFill(static_cast<float>(v));
		});

		auto *pickEdge = new QPushButton("Цвет рёбер...");
		pickEdge->setFocusPolicy(Qt::NoFocus);
		connect(pickEdge, &QPushButton::clicked, this, [this] {
			QColorDialog dlg;
			dlg.setOption(QColorDialog::DontUseNativeDialog, true);
			dlg.setOption(QColorDialog::ShowAlphaChannel, false);
			if (dlg.exec() == QDialog::Accepted) {
				if (viewport_) viewport_->SetEdgeColor(dlg.selectedColor());
			}
			if (viewport_) viewport_->setFocus();
		});

		auto addRow = [&](const QString &name, QWidget *ctrl) {
			auto *r = new QWidget(w);
			auto *rl = new QHBoxLayout(r);
			rl->setContentsMargins(0, 0, 0, 0);
			rl->setSpacing(8);
			auto *la = new QLabel(name, r);
			la->setMinimumWidth(130);
			rl->addWidget(la, 0);
			rl->addWidget(ctrl, 1);
			l->addWidget(r);
		};

		l->addWidget(dashed);
		addRow("Толщина:", width);
		addRow("Шаг пунктира:", dashPeriod);
		addRow("Заполнение:", dashFill);
		l->addWidget(pickEdge);

		l->addStretch(1);
		return w;
	}

	QWidget *Interface::MakePageVertices_() {
		auto *w = new QWidget();
		auto *l = new QVBoxLayout(w);
		l->setContentsMargins(0, 0, 0, 0);
		l->setSpacing(10);

		l->addWidget(new QLabel("Отображение вершин:", w));

		auto *mode = new QComboBox(w);
		mode->addItem("Отсутствует", 0);
		mode->addItem("Круг", 1);
		mode->addItem("Квадрат", 2);

		connect(mode, qOverload<int>(&QComboBox::currentIndexChanged), this, [this, mode](int) {
			if (!viewport_) return;
			viewport_->SetVertexMode(mode->currentData().toInt());
			viewport_->setFocus();
		});

		auto *size = new QDoubleSpinBox(w);
		size->setRange(1.0, 30.0);
		size->setSingleStep(1.0);
		size->setValue(6.0);
		connect(size, qOverload<double>(&QDoubleSpinBox::valueChanged), this, [this](double v) {
			if (viewport_) viewport_->SetVertexSize(static_cast<float>(v));
		});

		auto *pick = new QPushButton("Цвет вершин...", w);
		pick->setFocusPolicy(Qt::NoFocus);
		connect(pick, &QPushButton::clicked, this, [this] {
			QColorDialog dlg;
			dlg.setOption(QColorDialog::DontUseNativeDialog, true);
			dlg.setOption(QColorDialog::ShowAlphaChannel, false);
			if (dlg.exec() == QDialog::Accepted) {
				if (viewport_) viewport_->SetVertexColor(dlg.selectedColor());
			}
			if (viewport_) viewport_->setFocus();
		});

		l->addWidget(mode);

		auto *row = new QWidget(w);
		auto *rowL = new QHBoxLayout(row);
		rowL->setContentsMargins(0, 0, 0, 0);
		rowL->setSpacing(8);
		auto *lab = new QLabel("Размер:", row);
		lab->setMinimumWidth(130);
		rowL->addWidget(lab, 0);
		rowL->addWidget(size, 1);
		l->addWidget(row);

		l->addWidget(pick);

		l->addWidget(new QLabel("Масштаб модели:", w));

		auto* scaleRow = new QWidget(w);
		auto* scaleRowL = new QHBoxLayout(scaleRow);
		scaleRowL->setContentsMargins(0, 0, 0, 0);
		scaleRowL->setSpacing(8);

		auto* scaleLab = new QLabel("Масштаб:", scaleRow);
		scaleLab->setMinimumWidth(130);

		auto* bMinus = new QPushButton("−", scaleRow);
		auto* bPlus  = new QPushButton("+", scaleRow);
		bMinus->setFixedWidth(40);
		bPlus->setFixedWidth(40);
		bMinus->setFocusPolicy(Qt::NoFocus);
		bPlus->setFocusPolicy(Qt::NoFocus);

		auto* scaleSpin = new QDoubleSpinBox(scaleRow);
		scaleSpin->setRange(0.001, 15.0);
		scaleSpin->setDecimals(3);
		scaleSpin->setSingleStep(0.1);
		scaleSpin->setKeyboardTracking(false);
		scaleSpin->setFocusPolicy(Qt::ClickFocus);

		if (viewport_) scaleSpin->setValue(viewport_->GetScale());

		connect(scaleSpin, qOverload<double>(&QDoubleSpinBox::valueChanged), this,
		        [this](double v) {
		          if (viewport_) viewport_->SetScale(v);
		        });

		connect(bPlus, &QPushButton::clicked, this, [this]() {
		  if (viewport_) viewport_->ScaleBy(1.10);   // +10%
		  if (viewport_) viewport_->setFocus();
		});

		connect(bMinus, &QPushButton::clicked, this, [this]() {
		  if (viewport_) viewport_->ScaleBy(1.0 / 1.10); // -~9.09%
		  if (viewport_) viewport_->setFocus();
		});

		if (viewport_) {
		  connect(viewport_, &s21::GlWidget::ScaleChanged, this,
		          [scaleSpin](double s) {
		            const QSignalBlocker b(scaleSpin);
		            scaleSpin->setValue(s);
		          });
		}

		scaleRowL->addWidget(scaleLab, 0);
		scaleRowL->addWidget(bMinus, 0);
		scaleRowL->addWidget(scaleSpin, 1);
		scaleRowL->addWidget(bPlus, 0);

		l->addWidget(scaleRow);

		auto* resetScale = new QPushButton("Сброс масштаба", w);
		resetScale->setFocusPolicy(Qt::NoFocus);
		connect(resetScale, &QPushButton::clicked, this, [this]() {
		  if (viewport_) viewport_->ResetScale();
		  if (viewport_) viewport_->setFocus();
		});
		l->addWidget(resetScale);

		l->addStretch(1);
		return w;
	}

	QWidget *Interface::MakePageBackground_() {
		auto *w = new QWidget();
		auto *l = new QVBoxLayout(w);
		l->setContentsMargins(0, 0, 0, 0);
		l->setSpacing(10);

		l->addWidget(new QLabel("Цвет фона:", w));

		auto *pickBg = new QPushButton("Выбрать цвет...", w);
		pickBg->setFocusPolicy(Qt::NoFocus);
		connect(pickBg, &QPushButton::clicked, this, [this] {
			QColorDialog dlg;
			dlg.setOption(QColorDialog::DontUseNativeDialog, true);
			dlg.setOption(QColorDialog::ShowAlphaChannel, false);
			if (dlg.exec() == QDialog::Accepted) {
				if (viewport_) viewport_->SetBackgroundColor(dlg.selectedColor());
			}
			if (viewport_) viewport_->setFocus();
		});

		l->addWidget(pickBg);

		l->addSpacing(6);
		l->addWidget(new QLabel("Заливка граней:", w));

		auto *fillEnabled = new QCheckBox("Показывать заливку", w);
		fillEnabled->setChecked(false);
		fillEnabled->setChecked(viewport_ ? viewport_->GetFillEnabled() : false);

		if (viewport_) {
			connect(viewport_, &s21::GlWidget::FillEnabledChanged, this,
					[fillEnabled](bool on) {
						const QSignalBlocker b(fillEnabled);
						fillEnabled->setChecked(on);
					});
		}

		auto *fillTransparent = new QCheckBox("Прозрачная заливка", w);
		fillTransparent->setChecked(false);

		connect(fillEnabled, &QCheckBox::toggled, this,
				[this, fillTransparent](bool on) {
					if (viewport_) viewport_->SetFillEnabled(on);
					fillTransparent->setEnabled(on);
					if (viewport_) viewport_->setFocus();
				});

		connect(fillTransparent, &QCheckBox::toggled, this,
				[this](bool on) {
					if (!viewport_) return;
					viewport_->SetFillOpaque(!on);
					viewport_->setFocus();
				});

		l->addWidget(fillEnabled);
		l->addWidget(fillTransparent);

		auto *pickFill = new QPushButton("Цвет заливки...", w);
		pickFill->setFocusPolicy(Qt::NoFocus);

		connect(fillEnabled, &QCheckBox::toggled, this,
				[pickFill](bool on) { pickFill->setEnabled(on); });

		connect(pickFill, &QPushButton::clicked, this,
				[this, fillTransparent] {
					QColorDialog dlg;
					dlg.setOption(QColorDialog::DontUseNativeDialog, true);
					dlg.setOption(QColorDialog::ShowAlphaChannel, false);
					if (dlg.exec() == QDialog::Accepted) {
						const float a = fillTransparent->isChecked() ? 0.0f : 1.0f;
						if (viewport_) viewport_->SetFillColor(dlg.selectedColor(), a);
					}
					if (viewport_) viewport_->setFocus();
				});

		l->addWidget(pickFill);

		l->addStretch(1);
		return w;
	}

	QByteArray Interface::ImageToBytes_(const QImage &img,
										const char *format,
										int quality) {
		QByteArray data;
		QBuffer buffer(&data);
		if (!buffer.open(QIODevice::WriteOnly)) return {};

		if (!img.save(&buffer, format, quality)) return {};
		return data;
	}

	void Interface::SaveRenderedImage_() {
		if (!viewport_) return;

		const QImage frame = viewport_->grabFramebuffer(); // читает из FBO
		if (frame.isNull()) return;

		const QString filter = "BMP (*.bmp);;JPEG (*.jpg *.jpeg)"; // выбор пути пользователя (спецификация Qt)
		const QString path = QFileDialog::getSaveFileName(this,
														tr("Сохранить изображение"),
														QString(),
														filter);
		if (path.isEmpty()) return;

		const QString suf = QFileInfo(path).suffix().toLower();
		const char *fmt = nullptr;
		int quality = -1;

		if (suf == "bmp") {
			fmt = "BMP";
			quality = -1;
		} else if (suf == "jpg" || suf == "jpeg") {
			fmt = "JPG";
			quality = 95;
		} else {
			fmt = "BMP";
			quality = -1;
		}
		const bool ok = frame.save(path, fmt, quality); // берем из FBO
		if (!ok) {
			QMessageBox::warning(this, tr("Ошибка"),
								tr("Не удалось сохранить изображение в файл."));
		}

		if (viewport_) viewport_->setFocus();
	}

	QWidget *Interface::MakePageFile_() {
		auto *w = new QWidget();
		auto *l = new QVBoxLayout(w);
		l->setContentsMargins(0, 0, 0, 0);
		l->setSpacing(10);

		l->addWidget(new QLabel("Операции с файлами:", w));

		auto *openModel = new QPushButton("Загрузить модель...", w);
		auto *saveImg = new QPushButton("Сохранить изображение...", w);
		openModel->setFocusPolicy(Qt::NoFocus);
		saveImg->setFocusPolicy(Qt::NoFocus);

		connect(saveImg, &QPushButton::clicked, this, [this] {
			SaveRenderedImage_();
		});

		connect(openModel, &QPushButton::clicked, this, [this] {
			if (!viewport_) return;

			const QString path = QFileDialog::getOpenFileName(
				this,
				tr("Открыть модель"),
				QString(),
				tr("Wavefront OBJ (*.obj);;Все файлы (*.*)")
			);

			if (path.isEmpty()) {
				viewport_->setFocus();
				return;
			}

			const bool ok = viewport_->LoadModelFromObjFile(path);
			if (!ok) {
				QMessageBox::warning(this, tr("Ошибка"),
									tr("Не удалось загрузить модель. Проверьте формат OBJ и доступность файла."));
			}

			viewport_->setFocus();
		});

		l->addWidget(openModel);
		l->addWidget(saveImg);
		l->addStretch(1);
		return w;
	}
}

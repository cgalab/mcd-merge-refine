#ifdef ENABLE_VIEW

#include <cmath>
#include <QDesktopWidget>
#include "view.h"
#include "decomp.h"

MainWindow::MainWindow(QWidget *parent)
: QMainWindow(parent)
, m_scene{this} {
    ui.setupUi(this);

    ui.graphicsView->setScene(&m_scene);
    ui.graphicsView->installEventFilter(&m_navigation);
    ui.graphicsView->viewport()->installEventFilter(&m_navigation);
    ui.graphicsView->setRenderHint(QPainter::Antialiasing);

    QObject::connect(ui.actionToggleEdges, &QAction::triggered, this,
                     [this]() {
                         m_bShowEdges = !m_bShowEdges;
                         emit toggleShowEdges(m_bShowEdges);
                     });
    QObject::connect(ui.actionPrev, &QAction::triggered, this,
                     [this]() {
                         emit prev();
                     });
    QObject::connect(ui.actionNext, &QAction::triggered, this,
                     [this]() {
                         emit next();
                     });

    resize(1280, 1024);
    centerWidget();
}

void MainWindow::addItem(QGraphicsItem *item) {
    if (typeid(*item) == typeid(ConvexDecompGraphicsItem)) {
        const auto cdgi = static_cast<ConvexDecompGraphicsItem *> (item);
        QObject::connect(this, &MainWindow::next, cdgi, 
                         &ConvexDecompGraphicsItem::onNextEv);
        QObject::connect(this, &MainWindow::prev, cdgi, 
                         &ConvexDecompGraphicsItem::onPrevEv);
    }
    
    m_scene.addItem(item);
}

void MainWindow::centerWidget() {
    auto rect = geometry();
    rect.moveCenter(QApplication::desktop()->availableGeometry().center());
    setGeometry(rect);
}

#include "moc_view.cpp"

#endif

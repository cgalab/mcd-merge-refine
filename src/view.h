#ifndef VIEW_H
#define VIEW_H

#ifdef ENABLE_VIEW

#include <QGraphicsScene>
#include <QGraphicsItem>
#include <CGAL/Qt/GraphicsViewNavigation.h>
#include "ui_mcd-ref-merge.h"

class MainWindow : public QMainWindow {
    Q_OBJECT

public:
    MainWindow(QWidget *parent = nullptr);

    void addItem(QGraphicsItem *item);
    
signals:
    void toggleShowEdges(bool bIsVisible);
    void prev();
    void next();

private:
    void centerWidget();
    
    Ui::MainWindow ui;
    CGAL::Qt::GraphicsViewNavigation m_navigation;
    QGraphicsScene m_scene;
    bool m_bShowEdges{true};
};

#endif

#endif /* VIEW_H */


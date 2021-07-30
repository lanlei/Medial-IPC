#ifndef CUSTOMGRAPHICS_H
#define CUSTOMGRAPHICS_H

/**
 * 自定义多边形控件 作者:赵彦博(QQ:408815041 zyb920@hotmail.com) 2019-3-28
 * 1:自定义随意绘制多边形
 * 2:产生闭合形状后可单击选中移动整个多边形
 * 3:可拉动某个点
 * 4:支持多个多边形
 * 5:鼠标右键退出绘制
 * 6:可设置各种颜色
 */

#include <QWidget>

#ifdef quc
#if (QT_VERSION < QT_VERSION_CHECK(5,7,0))
#include <QtDesigner/QDesignerExportWidget>
#else
#include <QtUiPlugin/QDesignerExportWidget>
#endif

class QDESIGNER_WIDGET_EXPORT CustomGraphics : public QWidget
#else
class CustomGraphics : public QWidget
#endif

{
    Q_OBJECT
    Q_PROPERTY(bool selectDotVisible READ getSelectDotVisible WRITE setSelectDotVisible)
    Q_PROPERTY(int dotRadius READ getDotRadius WRITE setDotRadius)
    Q_PROPERTY(int lineWidth READ getLineWidth WRITE setLineWidth)

    Q_PROPERTY(QColor dotColor READ getDotColor WRITE setDotColor)
    Q_PROPERTY(QColor lineColor READ getLineColor WRITE setLineColor)
    Q_PROPERTY(QColor polygonColor READ getPolygonColor WRITE setPolygonColor)
    Q_PROPERTY(QColor selectColor READ getSelectColor WRITE setSelectColor)

public:
    typedef struct {
        QVector<QPoint> pos;
        bool selected;
    } Polygon;

    explicit CustomGraphics(QWidget *parent = 0);

protected:
    void mousePressEvent(QMouseEvent *e);
    void mouseMoveEvent(QMouseEvent *e);
    void mouseReleaseEvent(QMouseEvent *e);
    void paintEvent(QPaintEvent *);
    void drawPolygon(QPainter *p, const Polygon &v);
    void drawLines(QPainter *p, const QList<QPoint> &list, bool isFirst = true);

private:
    bool selectDotVisible;      //选中点可见
    int dotRadius;              //点的半径
    int lineWidth;              //线条宽度

    QColor dotColor;            //点的颜色
    QColor lineColor;           //线条颜色
    QColor polygonColor;        //多边形颜色
    QColor selectColor;         //选中颜色

    QPoint tempPoint;           //临时点
    QList<QPoint> tempPoints;   //点集合
    QList<Polygon> tempPolygons;//多边形集合

    bool pressed;               //鼠标是否按下
    QPoint lastPoint;           //鼠标按下处的坐标
    QPoint ellipsePos;          //保存按下点的坐标
    int selectedEllipseIndex;   //选中点的index
    Polygon pressedPolygon;     //保存按下时多边形的原始坐标
    int selectedIndex;          //选中多边形的index

private:
    //计算两点间的距离
    double length(const QPoint &p1, const QPoint &p2);
    //检测是否选中多边形
    bool checkPoint(const QVector<QPoint> &points, int x, int y);

public:
    bool getSelectDotVisible()  const;
    int getDotRadius()          const;
    int getLineWidth()          const;

    QColor getDotColor()        const;
    QColor getLineColor()       const;
    QColor getPolygonColor()    const;
    QColor getSelectColor()     const;

    QSize sizeHint()            const;
    QSize minimumSizeHint()     const;

public Q_SLOTS:
    void setSelectDotVisible(bool selectDotVisible);
    void setDotRadius(int dotRadius);
    void setLineWidth(int lineWidth);

    void setDotColor(const QColor &dotColor);
    void setLineColor(const QColor &lineColor);
    void setPolygonColor(const QColor &polygonColor);
    void setSelectColor(const QColor &selectColor);

    //清除临时绘制的
    void clearTemp();
    //清除所有
    void clearAll();
};

#endif // CUSTOMGRAPHICS_H

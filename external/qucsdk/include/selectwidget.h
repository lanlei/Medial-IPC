#ifndef SELECTWIDGET_H
#define SELECTWIDGET_H

/**
 * 描点跟随窗体控件 作者:feiyangqingyun(QQ:517216493) 2018-9-9
 * 1:可设置是否绘制描点
 * 2:可设置边距
 * 3:可设置描点颜色
 * 4:可设置描点尺寸
 * 5:可设置描点样式 正方形+圆形
 * 6:可设置选中边框宽度
 * 7:支持上下左右按键移动窗体
 * 8:支持delete键删除窗体
 * 9:支持八个描点改变窗体大小尺寸
 * 10:支持拖动窗体
 */

#include <QWidget>

#ifdef quc
#if (QT_VERSION < QT_VERSION_CHECK(5,7,0))
#include <QtDesigner/QDesignerExportWidget>
#else
#include <QtUiPlugin/QDesignerExportWidget>
#endif

class QDESIGNER_WIDGET_EXPORT SelectWidget : public QWidget
#else
class SelectWidget : public QWidget
#endif

{
    Q_OBJECT
    Q_ENUMS(PointStyle)

    Q_PROPERTY(bool drawPoint READ getDrawPoint WRITE setDrawPoint)
    Q_PROPERTY(int padding READ getPadding WRITE setPadding)
    Q_PROPERTY(int borderWidth READ getBorderWidth WRITE setBorderWidth)
    Q_PROPERTY(int pointSize READ getPointSize WRITE setPointSize)
    Q_PROPERTY(QColor pointColor READ getPointColor WRITE setPointColor)
    Q_PROPERTY(PointStyle pointStyle READ getPointStyle WRITE setPointStyle)

public:
    enum PointStyle {
        PointStyle_Rect = 0,        //正方形
        PointStyle_Circle = 1       //圆形
    };

    explicit SelectWidget(QWidget *parent = 0);

protected:
    bool eventFilter(QObject *watched, QEvent *event);
    void resizeEvent(QResizeEvent *);
    void mouseMoveEvent(QMouseEvent *);
    void paintEvent(QPaintEvent *);
    void drawRect(QPainter *painter);
    void drawCircle(QPainter *painter);
    void drawBorder(QPainter *painter);

private:
    bool drawPoint;                 //绘制描点
    int padding;                    //窗体的边距
    int borderWidth;                //边框宽度
    int pointSize;                  //描点的尺寸
    QColor pointColor;              //描点的颜色
    PointStyle pointStyle;          //描点的形状

    bool pressed;                   //鼠标按下
    bool pressedLeft;               //鼠标按下左侧
    bool pressedRight;              //鼠标按下右侧
    bool pressedTop;                //鼠标按下上侧
    bool pressedBottom;             //鼠标按下下侧
    bool pressedLeftTop;            //鼠标按下左上侧
    bool pressedRightTop;           //鼠标按下右上侧
    bool pressedLeftBottom;         //鼠标按下左下侧
    bool pressedRightBottom;        //鼠标按下右下侧

    int rectX, rectY, rectW, rectH; //窗体坐标+宽高
    QPoint lastPos;                 //鼠标按下处坐标

    QRectF rectLeft;                //左侧区域
    QRectF rectRight;               //右侧区域
    QRectF rectTop;                 //上侧区域
    QRectF rectBottom;              //下侧区域
    QRectF rectLeftTop;             //左上侧区域
    QRectF rectRightTop;            //右上侧区域
    QRectF rectLeftBottom;          //左下侧区域
    QRectF rectRightBottom;         //右下侧区域
    QWidget *widget;                //跟随的窗体

public:
    bool getDrawPoint()             const;
    int getPadding()                const;
    int getBorderWidth()            const;
    int getPointSize()              const;
    QColor getPointColor()          const;
    PointStyle getPointStyle()      const;

    QWidget *getWidget()            const;
    QSize sizeHint()                const;
    QSize minimumSizeHint()         const;

public Q_SLOTS:
    //设置是否绘制描点
    void setDrawPoint(bool drawPoint);
    //设置跟随窗体的边距
    void setPadding(int padding);
    //设置线条宽度
    void setBorderWidth(int borderWidth);
    //设置描点的尺寸
    void setPointSize(int pointSize);
    //设置描点的颜色
    void setPointColor(const QColor &pointColor);
    //设置描点的样式
    void setPointStyle(const PointStyle &pointStyle);

    //设置跟随窗体
    void setWidget(QWidget *widget);

Q_SIGNALS:
    void widgetPressed(QWidget *widget);
    void widgetRelease(QWidget *widget);
    void widgetDelete(QWidget *widget);
};

#endif // SELECTWIDGET_H

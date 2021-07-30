#ifndef SHADOWCLOCK_H
#define SHADOWCLOCK_H

/**
 * 光晕时钟控件 作者:雨田哥(QQ:3246214072) 整理:feiyangqingyun(QQ:517216493) 2019-10-07
 * 1:可设置圆弧半径宽度
 * 2:可设置光晕宽度
 * 3:可设置光晕颜色
 * 4:可设置文本颜色
 * 5:采用动画机制平滑进度展示时间
 */

#include <QWidget>

#ifdef quc
#if (QT_VERSION < QT_VERSION_CHECK(5,7,0))
#include <QtDesigner/QDesignerExportWidget>
#else
#include <QtUiPlugin/QDesignerExportWidget>
#endif

class QDESIGNER_WIDGET_EXPORT ShadowClock : public QWidget
#else
class ShadowClock : public QWidget
#endif

{
    Q_OBJECT
    Q_PROPERTY(int radiusWidth READ getRadiusWidth WRITE setRadiusWidth)
    Q_PROPERTY(int shadowWidth READ getShadowWidth WRITE setShadowWidth)

    Q_PROPERTY(QColor textColor READ getTextColor WRITE setTextColor)
    Q_PROPERTY(QColor shadowColor READ getShadowColor WRITE setShadowColor)

public:
    explicit ShadowClock(QWidget *parent = 0);
    ~ShadowClock();

protected:
    void paintEvent(QPaintEvent *);
    void drawArc(QPainter *painter, int radius, qreal angle);
    void drawText(QPainter *painter);

private:
    int radiusWidth;            //半径宽度
    int shadowWidth;            //光晕宽度

    QColor textColor;           //文本颜色
    QColor shadowColor;         //光晕颜色

public:
    int getRadiusWidth()        const;
    int getShadowWidth()        const;

    QColor getTextColor()       const;
    QColor getShadowColor()     const;

    QSize sizeHint()            const;
    QSize minimumSizeHint()     const;

public Q_SLOTS:
    //设置半径宽度+光晕宽度
    void setRadiusWidth(int radiusWidth);
    void setShadowWidth(int shadowWidth);

    //设置文本颜色+光晕颜色
    void setTextColor(const QColor &textColor);
    void setShadowColor(const QColor &shadowColor);
};

#endif // SHADOWCLOCK_H

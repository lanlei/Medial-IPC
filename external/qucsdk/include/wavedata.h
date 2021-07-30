#ifndef WAVEDATA_H
#define WAVEDATA_H

/**
 * 音量采样值波形控件 作者:feiyangqingyun(QQ:517216493) 2017-9-10
 * 1:可设置采样深度
 * 2:可设置当前位置线条宽度/线条颜色
 * 3:可设置前景色/背景色
 * 4:可设置数据展示样式,线条样式/柱状样式/平滑样式
 */

#include <QWidget>

#ifdef quc
#if (QT_VERSION < QT_VERSION_CHECK(5,7,0))
#include <QtDesigner/QDesignerExportWidget>
#else
#include <QtUiPlugin/QDesignerExportWidget>
#endif

class QDESIGNER_WIDGET_EXPORT WaveData : public QWidget
#else
class WaveData : public QWidget
#endif

{
    Q_OBJECT
    Q_ENUMS(WaveStyle)
    Q_PROPERTY(double deep READ getDeep WRITE setDeep)
    Q_PROPERTY(bool showLine READ getShowLine WRITE setShowLine)
    Q_PROPERTY(int lineWidth READ getLineWidth WRITE setLineWidth)
    Q_PROPERTY(QColor lineColor READ getLineColor WRITE setLineColor)
    Q_PROPERTY(QColor foreground READ getForeground WRITE setForeground)
    Q_PROPERTY(QColor background READ getBackground WRITE setBackground)
    Q_PROPERTY(WaveStyle waveStyle READ getWaveStyle WRITE setWaveStyle)

public:
    enum WaveStyle {
        WaveStyle_Line = 0,     //线条样式
        WaveStyle_Smooth = 1,   //平滑样式
        WaveStyle_Bar = 2       //柱状样式
    };

    explicit WaveData(QWidget *parent = 0);

protected:
    void mousePressEvent(QMouseEvent *);
    void paintEvent(QPaintEvent *);
    void drawBg(QPainter *painter);
    void drawData(QPainter *painter);
    void drawLine(QPainter *painter);

private:    
    double deep;                //采集深度
    bool showLine;              //显示线条
    int lineWidth;              //线条宽度
    QColor lineColor;           //线条颜色
    QColor foreground;          //前景色
    QColor background;          //背景色
    WaveStyle waveStyle;        //数据样式

    int length;                 //采样点长度
    int position;               //当前位置
    QVector<float> data;        //采样点数据

public:    
    double getDeep()            const;
    bool getShowLine()          const;
    int getLineWidth()          const;
    QColor getLineColor()       const;
    QColor getForeground()      const;
    QColor getBackground()      const;
    WaveStyle getWaveStyle()    const;

    QSize sizeHint()            const;
    QSize minimumSizeHint()     const;

public slots:
    //设置深度
    void setDeep(double deep);

    //设置是否显示线条
    void setShowLine(bool showLine);
    //设置线条宽度
    void setLineWidth(int lineWidth);
    //设置线条颜色
    void setLineColor(const QColor &lineColor);

    //设置前景色
    void setForeground(const QColor &foreground);
    //设置背景色
    void setBackground(const QColor &background);

    //设置数据样式
    void setWaveStyle(const WaveStyle &waveStyle);

    //设置总长度
    void setLength(int length);
    //设置当前位置
    void setPosition(int position);

    //设置当前数据
    void setData(const QVector<float> &data);
    //清空数据
    void clearData();

signals:
    void positionChanged(int position);
};

#endif // WAVEDATA_H

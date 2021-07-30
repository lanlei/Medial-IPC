#ifndef SMOOTHCURVECREATOR_H
#define SMOOTHCURVECREATOR_H

#include <QList>
#include <QPointF>
#include <QPainterPath>

class SmoothCurveCreator
{
public:
    static QPainterPath createSmoothCurve(const QVector<QPointF> &points);

private:
	static void calculateFirstControlPoints(double *&result, const double *rhs, int n);
    static void calculateControlPoints(const QVector<QPointF> &knots, QVector<QPointF> *firstControlPoints, QVector<QPointF> *secondControlPoints);
};

#endif // SMOOTHCURVECREATOR_H

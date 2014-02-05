#ifndef GLOSSATRON_H
#define GLOSSATRON_H

#include <QObject>
#include <QList>
#include <QRect>
#include <QVariant>
#include <QStringList>

QT_BEGIN_NAMESPACE
class Image;
class QXmlStreamWriter;
class QString;
class QLine;
QT_END_NAMESPACE

#include "../../PalatoglossatronQt/interfaces.h"
#include "fftw3.h"

class GlossatronPlugin : public QObject, public TrackingInterface
{
    Q_OBJECT
    Q_PLUGIN_METADATA(IID "palatoglossatron.qt.trackinginterface/1.0" FILE "glossatron.json")
    Q_INTERFACES(TrackingInterface)

public:
    GlossatronPlugin();
    ~GlossatronPlugin();

    QString name() const;
    bool initialized();
    bool needsGrid();
    void setGrid(QList<QLine>* grid);
    QList<QPoint*> trace(QImage *img);
    void receiveXMLSettings(QString xml);
    void createFilterFFT();
    QString settingsToXML();

    QStringList settingsTraceNames() const;
    bool setSettingTrace(int i, QList<QPoint*>* trace);

    void settings();

    QList<QList<QPoint*>*> gridPoints;

private:
    double sobel(fftw_complex *img, quint32 i, quint32 j);

    quint32 nx, ny; // FFT parameters
    QRect roi;
    bool sobelStep;
    bool success;

    double *gaussianLaplacian;
    fftw_plan theplan, nalpeht, filterPlan;
    void setFFTImageSizes(quint32 x, quint32 y);
    fftw_complex *in, *filterResult;
    fftw_complex *out;
    fftw_complex *filter_512_1024;
    fftw_complex *filter_other;

    QList<int> upperBound, lowerBound;

    QStringList settingsLabels;
    QList<QVariant> settingsValues;

    QStringList settingsTraces;
};

#endif

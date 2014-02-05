#include <QtGui>

#include <math.h>
#include <stdlib.h>
#include <QMessageBox>
#include <QImage>
#include <QXmlStreamWriter>
#include <QString>
#include <QtDebug>
#include <QLine>
#include <QTime>
#include <QRect>

#include "glossatron.h"
#include "dataentrywidget.h"

#define IMG(X,Y) ( qRed( *(((QRgb*)img->scanLine(Y)) + X) ) )
#define IMAGE(X,Y) (*(image+width*(Y)+(X)))
#define GAUSSIAN(X,Y) (  *(((QRgb*)gaussian.scanLine(Y)) + X)  )
#define Z(i) (*(buf+((i-1)*3)+2))
#define GRD(i,j) (gridPoints.at(i)->at(j))
#define GRD_X(i,j) (gridPoints.at(i)->at(j)->x())
#define GRD_Y(i,j) (gridPoints.at(i)->at(j)->y())

GlossatronPlugin::GlossatronPlugin()
{
    settingsTraces << "Upper bound";
    settingsTraces << "Lower bound";

    settingsLabels << "Perform Sobel Step (yes/no)";
    settingsValues << "yes";

    success = true;

    filter_other = 0;

    roi.setTopLeft(QPoint(150,62));
    roi.setBottomRight(QPoint(709,393));

    setFFTImageSizes(480,720);

    // initialize FFTW
    in = (fftw_complex*)fftw_malloc(sizeof(fftw_complex)*nx*ny);
    out = ( fftw_complex* )fftw_malloc(sizeof(fftw_complex)*nx*ny);
    filterResult = (fftw_complex*)fftw_malloc(sizeof(fftw_complex)*nx*ny);
    gaussianLaplacian = (double*)fftw_malloc(sizeof(double)*57*57);
    if(in==NULL || out==NULL || filterResult==NULL || gaussianLaplacian==NULL)
    {
	success = false;
	return;
    }

    QFile gaussianLaplacianFile(":/gaulap_double.bin");
    if(!gaussianLaplacianFile.open(QIODevice::ReadOnly ))
    {
	qDebug() << "Gaussian/Laplacian filter could not be opened.";
	success = false;
	return;
    }
    else
    {
	QDataStream gaussianStream(&gaussianLaplacianFile);
	gaussianStream.setByteOrder(QDataStream::LittleEndian);
	gaussianStream.setFloatingPointPrecision(QDataStream::DoublePrecision);

	for(quint32 i=0; i < 57*57; i++)
	{
	    gaussianStream >> *(gaussianLaplacian + i);
	}
    }

    theplan = fftw_plan_dft_2d(nx,ny,in,out, 1, FFTW_ESTIMATE );
    nalpeht = fftw_plan_dft_2d(nx,ny,out,filterResult, -1 , FFTW_ESTIMATE );

    // store the filter size we're most likely to encounter
    filter_512_1024 = (fftw_complex*)fftw_malloc(sizeof(fftw_complex)*512*1024);
    if(filter_512_1024==NULL)
    {
	success = false;
	return;
    }
    QFile filter_512_1024_file(":/filter_512_1024.bin");
    if(!filter_512_1024_file.open(QIODevice::ReadOnly ))
    {
	qDebug() << "filter_512_1024_file could not be opened. Performance will likely be somewhat slower at first.";
	fftw_free( filter_512_1024 );
	filter_512_1024 = NULL;
    }
    else
    {
	QDataStream filterStream(&filter_512_1024_file);
	double tmp;
	filterStream.setByteOrder(QDataStream::LittleEndian);
	filterStream.setFloatingPointPrecision(QDataStream::DoublePrecision);

	for(quint32 i=0; i < 512*1024; i++)
	{
	    filterStream >> tmp;
	    (filter_512_1024+i)[0][0] = tmp;
	    filterStream >> tmp;
	    (filter_512_1024+i)[0][1] = tmp;
	}
    }

}

GlossatronPlugin::~GlossatronPlugin()
{
    fftw_destroy_plan(theplan);
    fftw_destroy_plan(filterPlan);
    fftw_destroy_plan(nalpeht);

    fftw_free( in );
    fftw_free( out );
    fftw_free( filterResult );
    free(gaussianLaplacian);
    if(filter_other != 0)
    {
	fftw_free( filter_other );
    }
}

QString GlossatronPlugin::name() const
{
    return "Glossatron";
}

bool GlossatronPlugin::initialized()
{
    return success;
}

QList<QPoint*> GlossatronPlugin::trace(QImage *img)
{
    QList<QPoint*> list;

    // so, nx is the number of rows, ny the number of columns
    if( nx < (quint32)img->height() || ny < (quint32)img->width() ) { setFFTImageSizes(img->height(), img->width()); }

    quint32 i,j;
    for(i=0; i< (quint32)(img->height()); i++)
    {
	for(j=0; j< (quint32)(img->width()); j++)
	{
	    (in+ny*i+j)[0][0] = qRed(*(((QRgb*)img->scanLine(i)) + j));
	    (in+ny*i+j)[0][1] = 0;
	}
	for(j=img->width(); j< ny; j++)
	{
	    (in+ny*i+j)[0][0] = 0;
	    (in+ny*i+j)[0][1] = 0;
	}
    }
    for(i=img->height(); i< nx; i++)
    {
	for(j=0; j< ny; j++)
	{
	    (in+ny*i+j)[0][0] = 0;
	    (in+ny*i+j)[0][1] = 0;
	}
    }

    theplan = fftw_plan_dft_2d(nx,ny,in,out, 1, FFTW_ESTIMATE );
    fftw_execute(theplan);

    fftw_complex *filter;

    qint32 offset;

    if( nx == 512 && ny == 1024 && filter_512_1024!=NULL )
    {
	filter = filter_512_1024;
	offset = -28; // bit of a hack?
    }
    else
    {
	createFilterFFT();
	filter = filter_other;
	offset = 28; // bit of a hack?
    }
    //    filter = filter_512_1024;

    double tmpReal;
    // here do the multiplication
    for(i=0; i<nx; i++)
    {
	for(j=0; j<ny; j++)
	{
	    tmpReal = (out+ny*i+j)[0][0] * (filter+ny*i+j)[0][0] - (out+ny*i+j)[0][1] * (filter+ny*i+j)[0][1];
	    (out+ny*i+j)[0][1] = (out+ny*i+j)[0][0] * (filter+ny*i+j)[0][1] + (out+ny*i+j)[0][1] * (filter+ny*i+j)[0][0];
	    (out+ny*i+j)[0][0] = tmpReal;
	}
    }

    nalpeht = fftw_plan_dft_2d(nx,ny,out,filterResult, -1 , FFTW_ESTIMATE );
    fftw_execute(nalpeht);

    double max;
    qint32 maxind;
    for(i=0; i<(unsigned)gridPoints.count(); i++) // for every radius
    {
	max = 0;
	maxind = -1;

	if( lowerBound.count() > 0 && lowerBound.count() == upperBound.count() )
	{
	    for(j= lowerBound.at(i); j< (unsigned)upperBound.at(i); j++)
	    {
		if( (filterResult + ny*(GRD_Y(i,j) + offset) + GRD_X(i,j) + offset )[0][0] > max )
		{
		    maxind = j;
		    max = (filterResult + ny*(GRD_Y(i,j) + offset) + GRD_X(i,j) + offset )[0][0];
		}
	    }
	}
	else
	{
	    for(j= 0; j< (unsigned)gridPoints.at(i)->count(); j++)
	    {
		if( (filterResult + ny*(GRD_Y(i,j) + offset) + GRD_X(i,j) + offset )[0][0] > max )
		{
		    maxind = j;
		    max = (filterResult + ny*(GRD_Y(i,j) + offset) + GRD_X(i,j) + offset )[0][0];
		}
	    }
	}


	if(maxind != -1)
	{
	    if(sobelStep)
	    {
		while( maxind >= 0 &&  sobel(filterResult, GRD_Y(i,maxind) , GRD_X(i,maxind) ) < sobel(filterResult, GRD_Y(i,maxind-1) , GRD_X(i,maxind-1) ) )
		{
		    //		    qDebug() << i << "sobel adjustment";
		    maxind--;
		}
	    }
	    list << new QPoint( gridPoints.at(i)->at(maxind)->x(), gridPoints.at(i)->at(maxind)->y() );
	}
	else
	{
	    list << new QPoint(-1,-1);
	}
    }

    /*
    fid = fopen("in.bin","wb");
    fwrite(in,sizeof(double),2*nx*ny,fid);
    fclose(fid);

    fid = fopen("filterResult.bin","wb");
    fwrite(filterResult,sizeof(double),2*nx*ny,fid);
    fclose(fid);

    FILE *fid = fopen("filterResult.bin","wb");
    fwrite(filterResult,sizeof(double),2*nx*ny,fid);
    fclose(fid);
*/

    //    qDebug() << "done";
    return list;
}

void GlossatronPlugin::createFilterFFT()
{
    filter_other = (fftw_complex*)fftw_malloc(sizeof(fftw_complex)*nx*ny);
    fftw_complex *padded_filter = (fftw_complex*)fftw_malloc(sizeof(fftw_complex)*nx*ny);
    if(padded_filter==NULL || filter_other==NULL) { return; }

    quint32 i,j;
    for(i=0; i< 57; i++)
    {
	for(j=0; j< 57; j++)
	{
	    (padded_filter+ny*i+j)[0][0] = *(gaussianLaplacian + 57*i + j );
	    (padded_filter+ny*i+j)[0][1] = 0;
	}
	for(j=57; j< ny; j++)
	{
	    (padded_filter+ny*i+j)[0][0] = 0;
	    (padded_filter+ny*i+j)[0][1] = 0;
	}
    }
    for(i=57; i< nx; i++)
    {
	for(j=0; j< ny; j++)
	{
	    (padded_filter+ny*i+j)[0][0] = 0;
	    (padded_filter+ny*i+j)[0][1] = 0;
	}
    }
    /*
    FILE *fid = fopen("paddedFilter.bin","wb");
    fwrite(padded_filter,sizeof(double),2*nx*ny,fid);
    fclose(fid);
*/
    filterPlan = fftw_plan_dft_2d(nx,ny,padded_filter,filter_other, 1, FFTW_ESTIMATE );
    fftw_execute(filterPlan);
    /*
    FILE *fid = fopen("filter_other.bin","wb");
    fwrite(filter_other,sizeof(double),2*nx*ny,fid);
    fclose(fid);

    QDataStream filterStream(new QFile("filter_512_1024.bin"));
    filterStream.setByteOrder(QDataStream::LittleEndian);
    filterStream.setFloatingPointPrecision(QDataStream::DoublePrecision);

    for(quint32 i=0; i < 512*1024; i++)
    {
	filterStream << (filter_other+i)[0][0];
	filterStream << (filter_other+i)[0][1];
    }
*/

    fftw_free(padded_filter);
}

double GlossatronPlugin::sobel(fftw_complex *img, quint32 i, quint32 j)
{
    double one, two;

    // http://en.wikipedia.org/wiki/Sobel_filter
    one = -1*(img+(i-1)*ny+(j-1))[0][0] + -2*(img+(i-1)*ny+(j))[0][0] + -1*(img+(i-1)*ny+(j+1))[0][0] +
	  1*(img+(i+1)*ny+(j-1))[0][0] + 2*(img+(i+1)*ny+(j))[0][0] + 1*(img+(i+1)*ny+(j+1))[0][0];
    two = -1*(img+(i-1)*ny+(j-1))[0][0] + -1*(img+(i-1)*ny+(j+1))[0][0] +
	  -2*(img+(i)*ny+(j-1))[0][0] + 2*(img+(i)*ny+(j+1))[0][0] +
	  -1*(img+(i+1)*ny+(j-1))[0][0] + 1*(img+(i+1)*ny+(j+1))[0][0];

    return sqrt( pow(one,2) + pow(two,2) );
}

bool GlossatronPlugin::needsGrid()
{
    return true;
}

void GlossatronPlugin::setGrid(QList<QLine>* grid)
{
    int i, j;
    int dx, dy;
    float m, b;
    QPoint p0, p1, tmp;

    //    QList<QList<QPoint*>*> gridPoints;
    for(i=0; i<gridPoints.count(); i++)
    {
	qDeleteAll(gridPoints.at(i)->begin(), gridPoints.at(i)->end());
	gridPoints.at(i)->clear();
    }
    qDeleteAll(gridPoints.begin(), gridPoints.end());
    gridPoints.clear();

    for(i=0; i<grid->count(); i++)
    {
	gridPoints << new QList<QPoint*>;

	// http://www.cs.unc.edu/~mcmillan/comp136/Lecture6/Lines.html
	if( grid->at(i).x1() > grid->at(i).x2() ) { p1 = grid->at(i).p1(); p0 = grid->at(i).p2(); }
	else { p0 = grid->at(i).p1(); p1 = grid->at(i).p2(); }

	// vertical lines
	if( p0.x() == p1.x() )
	{
	    if( p0.y() < p1.y() ) { tmp=p1; p1=p0; p0=tmp; } // switch them so that p2 is lower
	    for(j=p1.y() ; j<= p0.y(); j++)
	    {
		*(gridPoints.last()) << new QPoint(p0.x(),j);
	    }
	}
	// horizontal lines
	else if( p0.y() == p1.y() )
	{
	    if( p0.x() < p1.x() ) { tmp=p1; p1=p0; p0=tmp; } // switch them so that p2 is lower
	    for(j=p1.x() ; j<= p0.x(); j++)
	    {
		*(gridPoints.last()) << new QPoint(j,p0.y());
	    }
	}
	else // sloped lines
	{
	    dx = p1.x() - p0.x(); dy = p1.y() - p0.y();
	    if(abs(dx) > abs(dy))
	    {
		m = (float)dy / (float)dx;
		b = p0.y() - m*p0.x();
		dx = (dx < 0) ? -1 : 1;
		while (p0.x() != p1.x())
		{
		    p0.rx() += dx;
		    *(gridPoints.last()) << new QPoint( p0.x(), round( m*p0.x() + b ) );
		}
	    }
	    else
	    {
		m = (float)dx / (float)dy;
		b = p0.x() - m*p0.y();
		dy = (dy < 0) ? -1 : 1;
		while (p0.y() != p1.y())
		{
		    p0.ry() += dy;
		    *(gridPoints.last()) << new QPoint( round(m*p0.y() + b) , p0.y() );
		}

	    }
	}
    }
    /*
    QImage test(QSize(720,480),QImage::Format_Mono);
    test.fill(0);
    for(i=0; i<gridPoints.count(); i++)
    {
	for(j=0; j<gridPoints.at(i)->count(); j++)
	{
	    test.setPixel(*(gridPoints.at(i)->at(j)),1);
	}
    }
    test.save("radii.png","png");
*/
}

void GlossatronPlugin::setFFTImageSizes(quint32 x, quint32 y)
{
    quint32 i;

    i=1;
    do
    {
	nx = (quint32)floor(pow(2,i));
	i++;
    }
    while( nx < x );

    i=1;
    do
    {
	ny = (quint32)floor(pow(2,i));
	i++;
    }
    while( ny < y );
}

QString GlossatronPlugin::settingsToXML()
{
    //    qDebug() << "GlossatronPlugin::settingsToXML";

    QString theXml;
    QXmlStreamWriter writer(&theXml);

    writer.writeTextElement("sobelStep",QString::number(sobelStep));

    writer.writeStartElement("roi");

    writer.writeEmptyElement("topLeft");
    writer.writeAttribute("x",QString::number(roi.topLeft().x()));
    writer.writeAttribute("y",QString::number(roi.topLeft().y()));

    writer.writeEmptyElement("bottomRight");
    writer.writeAttribute("x",QString::number(roi.bottomRight().x()));
    writer.writeAttribute("y",QString::number(roi.bottomRight().y()));

    writer.writeEndElement(); // roi

    for(int i=0; i< upperBound.count(); i++)
	writer.writeTextElement("upper-bound-index",QString::number(upperBound.at(i)));

    for(int i=0; i< lowerBound.count(); i++)
	writer.writeTextElement("lower-bound-index",QString::number(lowerBound.at(i)));

    return theXml;
}

void GlossatronPlugin::receiveXMLSettings(QString xml)
{
    QString name;

    xml = "<wrapper>" + xml + "</wrapper>";

    QXmlStreamReader reader(xml);
    QXmlStreamAttributes attr;

    while (!reader.atEnd()) {
	if(reader.readNext() == QXmlStreamReader::StartElement)
	{
	    name = reader.name().toString();

	    if(name == "sobelStep") {
		if(reader.readElementText() == "1")
		{
		    sobelStep = true;
		}
		else
		{
		    sobelStep = false;
		}
	    }
	    if(name == "topLeft") {
		attr = reader.attributes();
		if(attr.hasAttribute("x") && attr.hasAttribute("y"))
		{
		    roi.setTopLeft(QPoint( attr.value("x").toString().toInt() , attr.value("y").toString().toInt() ));
		}
		else { qDebug() << "Line " << reader.lineNumber() << ", Column " << reader.columnNumber() << ": " << "Error in topLeft tag. No x or y attribute."; return; }
	    }
	    if(name == "bottomRight") {
		attr = reader.attributes();
		if(attr.hasAttribute("x") && attr.hasAttribute("y"))
		{
		    roi.setBottomRight(QPoint( attr.value("x").toString().toInt() , attr.value("y").toString().toInt() ));
		}
		else { qDebug() << "Line " << reader.lineNumber() << ", Column " << reader.columnNumber() << ": " << "Error in bottomRight tag. No x or y attribute."; return; }
	    }
	    if(name == "upper-bound-index") {
		upperBound << reader.readElementText().toInt();
	    }
	    if(name == "lower-bound-index") {
		lowerBound << reader.readElementText().toInt();
	    }
	}
    }

    return;
}

void GlossatronPlugin::settings()
{
    DataEntryWidget dew(&settingsLabels, &settingsValues, "", 0);
    if( dew.exec() == QDialog::Accepted)
    {	
	for(int i=0; i<settingsValues.count(); i++)
	{
	    settingsValues.replace(i, dew.values.at(i));
	}
    }
}

QStringList GlossatronPlugin::settingsTraceNames() const
{
    return settingsTraces;
}

bool GlossatronPlugin::setSettingTrace(int i, QList<QPoint*>* trace)
{
    //    qDebug() << i << trace;
    //    QList<QList<QPoint*>*> gridPoints;

    if(trace->count() != gridPoints.count() ) { qDebug() << "Unequal number of trace points and grid lines."; return false; }

    if( i == 0 ) // upper
    {
	upperBound.clear();
	for(int j = 0; j < trace->count(); j++) // cycle through the grid points
	{
	    int index = -1;
	    for(int k=0; k<gridPoints.at(j)->count(); k++)
	    {
//		if( *(gridPoints.at(j)->at(k)) == *(trace->at(j)) )
		// check for approximate equality since these things are subject to rounding errors

		QPoint tmp = *(gridPoints.at(j)->at(k)) - *(trace->at(j));
		if(tmp.manhattanLength() <= 2)
		{
		    index = k;
		    break;
		}
	    }
	    if( index == -1 )
		return false;
	    else
		upperBound << index;
	}

	if( lowerBound.count() == upperBound.count() )
	{
	    for( int j=0; j < lowerBound.count(); j++ )
	    {
		if( lowerBound.at(j) > upperBound.at(j) )
		{
		    int tmp = lowerBound.at(j);
		    lowerBound[j] = upperBound.at(j);
		    upperBound[j] = tmp;
		}
	    }
	}

	return true;
    }
    else if (i == 1) // lower
    {
	lowerBound.clear();
	for(int j = 0; j < trace->count(); j++) // cycle through the grid points
	{
	    int index = -1;
	    for(int k=0; k<gridPoints.at(j)->count(); k++)
	    {
		QPoint tmp = *(gridPoints.at(j)->at(k)) - *(trace->at(j));
		if(tmp.manhattanLength() <= 2)
		{
		    index = k;
		    break;
		}
	    }

	    if( index == -1 )
		return false;
	    else
		lowerBound << index;
	}

	if( lowerBound.count() == upperBound.count() )
	{
	    for( int j=0; j < lowerBound.count(); j++ )
	    {
		if( lowerBound.at(j) > upperBound.at(j) )
		{
		    int tmp = lowerBound.at(j);
		    lowerBound[j] = upperBound.at(j);
		    upperBound[j] = tmp;
		}
	    }
	}
/*
	qDebug() << "lowerBound";
	for(int k=0; k < lowerBound.count(); k++)
	    qDebug() << lowerBound.at(k);
*/
	return true;
    }

    return false;
}

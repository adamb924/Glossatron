#ifndef DATAENTRYWIDGET_H
#define DATAENTRYWIDGET_H

#include <QDialog>
#include <QWidget>
#include <QList>
#include <QLineEdit>
#include <QList>
#include <QVariant>

class QStringList;

class DataEntryWidget : public QDialog
{
Q_OBJECT
public:
    explicit DataEntryWidget(QStringList* f, QList<QVariant>* v, QString label, QWidget *parent);

    QList<QLineEdit*> edits;
    QList<QVariant> values;

signals:

public slots:
    void accept();

};

#endif // DATAENTRYWIDGET_H

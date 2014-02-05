TEMPLATE = lib
QT += widgets gui
CONFIG += plugin
# INCLUDEPATH += /Users/Adam/Documents/QtWork/PalatoglossatronQt
HEADERS = glossatron.h \
    dataentrywidget.h
SOURCES = glossatron.cpp \
    dataentrywidget.cpp
TARGET = $$qtLibraryTarget(pgqt_glossatron)
# DESTDIR = /Users/Adam/Documents/QtWork/PalatoglossatronQt/plugins
RESOURCES += glossatron.qrc
LIBS += -L./ \
    -lfftw3-3 \
    -lm

#-------------------------------------------------
#
# Project created by QtCreator 2014-12-08T18:07:53
#
#-------------------------------------------------

QT       += core gui opengl widgets declarative

greaterThan(QT_MAJOR_VERSION, 4): QT += widgets

TARGET = flow_sim
TEMPLATE = app

INCLUDEPATH += C:\tools\boost_1_56_0

LIBS += "-LC:\tools\boost_1_56_0\stage\lib"


SOURCES += main.cpp\
        mainwindow.cpp \
    ui_vals.cpp \
    mesh.cpp \
    glwidget.cpp \
    scene.cpp \
    sail.cpp \
    mast.cpp

HEADERS  += mainwindow.h \
    ui_vals.h \
    mesh.h \
    glwidget.h \
    scene.h \
    sail.h \
    mast.h

FORMS    += mainwindow.ui

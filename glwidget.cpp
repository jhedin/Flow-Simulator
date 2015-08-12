#include "glwidget.h"
#include "ui_vals.h"
#include "scene.h"

GLWidget::GLWidget( QWidget *parent ) : QGLWidget( parent ) {
    // ADD YOUR CODE
    sc = new scene();
}

GLWidget::~GLWidget() {

}

void GLWidget::initializeGL() {

    glClearColor(1.0f, 1.0f, 1.0f, 1.0f);
    glShadeModel( GL_FLAT );
    glMatrixMode( GL_MODELVIEW );
    glLoadIdentity();
}

void GLWidget::paintGL() {
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    //glLoadIdentity();
    glLineWidth(1);
    glColor3f(0.0, 0.0, 0.0);
    glEnable(GL_LINE_SMOOTH);
    sc->render();

}

void GLWidget::resizeGL( int width, int height ) {
    glMatrixMode( GL_PROJECTION );
    glLoadIdentity();
    glOrtho(-1.0,user->x_len + 1,-1.0,user->y_len + 1,-10.0,10.0);
    glFlush();
    glMatrixMode(GL_MODELVIEW);
    glViewport( 0, 0, (GLint)width, (GLint)height );
    h = height;
    w = width;
}

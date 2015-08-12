#ifndef GLWIDGET_H
#define GLWIDGET_H

#include <QWidget>

#include <QtOpenGL/QGLWidget>
#include <GL/gl.h>
#include <GL/glu.h>

class scene;
class ui_vals;

class GLWidget : public QGLWidget {
    Q_OBJECT

public:
    /// CONSTRUCTOR
    GLWidget(QWidget *parent = 0);

    /// DESTRUCTOR
    ~GLWidget();

    QTimer* timer;
    scene * sc;
    ui_vals * user;
    // ADD YOUR CODE

private:
    int h,w;


protected:
    /// OPENGL
    void initializeGL();
    void paintGL();
    void resizeGL( int width, int height );

signals:
    void sizeChanged(double width, double height);

public slots:
    // ADD YOUR CODE

private slots:


protected slots:


};

#endif // GLWIDGET_H

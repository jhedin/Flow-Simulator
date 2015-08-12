#include "mainwindow.h"
#include "ui_mainwindow.h"
#include "scene.h"
#include "mesh.h"
#include "mast.h"
#include "sail.h"
#include "glwidget.h"
#include "ui_vals.h"

MainWindow::MainWindow(QWidget *parent) :
    QMainWindow(parent),
    ui(new Ui::MainWindow)
{
     glWidget = new GLWidget();
     user = new ui_vals();

     glWidget->user = user;
     glWidget->sc->user = user;
     glWidget->sc->ms->user = user;
     glWidget->sc->mt->user = user;
     glWidget->sc->s->user = user;

     glWidget->sc->init();

    ui->setupUi(this);
    ui->VIEW->setWidget(glWidget);

    time = 0;
    rendered_time = 0;

    startTimer(33);//(int)(user->dt * 10000.0));
}

MainWindow::~MainWindow()
{
    delete ui;
}

void MainWindow::timerEvent(QTimerEvent *)
{
    rendered_time += 1;

    while(time < rendered_time / 1000.0)
    {
        glWidget->sc->step();
        time += user->dt;
    }

    glWidget->updateGL();
}

void MainWindow::on_reset_clicked()
{
    glWidget->sc->init();
}

void MainWindow::on_dt_valueChanged(double arg1)
{
    user->dt = arg1;
    glWidget->sc->init();
}

void MainWindow::on_ul_valueChanged(double arg1)
{
    user->ul = arg1;
    glWidget->sc->init();
}

void MainWindow::on_anu_valueChanged(double arg1)
{
    user->anu = arg1;
    glWidget->sc->init();
}

void MainWindow::on_i_bar_valueChanged(double arg1)
{
    user->dt = arg1;
    glWidget->sc->init();
}

void MainWindow::on_j_bar_valueChanged(double arg1)
{
    user->j_bar = arg1;
    glWidget->sc->init();
}

void MainWindow::on_x_len_valueChanged(double arg1)
{
    user->x_len = arg1;
    glWidget->sc->init();
}

void MainWindow::on_y_len_valueChanged(double arg1)
{
    user->y_len = arg1;
    glWidget->sc->init();
}

void MainWindow::on_P0_valueChanged(double arg1)
{
    user->P0 = arg1;
    glWidget->sc->init();
}

void MainWindow::on_u0_valueChanged(double arg1)
{
    user->u0 = arg1;
    glWidget->sc->init();
}

void MainWindow::on_v0_valueChanged(double arg1)
{
    user->v0 = arg1;
    glWidget->sc->init();
}

void MainWindow::on_n_lines_valueChanged(double arg1)
{
     user->streamlines = arg1;
     glWidget->sc->init();
}

void MainWindow::on_mast_x_valueChanged(double arg1)
{
    user->mast_x = arg1;
    glWidget->sc->init();
}

void MainWindow::on_mast_y_valueChanged(double arg1)
{
    user->mast_y = arg1;
    glWidget->sc->init();
}

void MainWindow::on_mast_r_valueChanged(double arg1)
{
    user->mast_r= arg1;
    glWidget->sc->init();
}

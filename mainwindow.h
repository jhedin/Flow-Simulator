#ifndef MAINWINDOW_H
#define MAINWINDOW_H

#include <QMainWindow>

class GLWidget;
class ui_vals;

namespace Ui {
class MainWindow;
}

class MainWindow : public QMainWindow
{
    Q_OBJECT

public:
    explicit MainWindow(QWidget *parent = 0);
    ~MainWindow();

    GLWidget * glWidget;
    ui_vals * user;

protected:
    void timerEvent(QTimerEvent * );

private slots:
    void on_reset_clicked();

    void on_dt_valueChanged(double arg1);

    void on_ul_valueChanged(double arg1);

    void on_anu_valueChanged(double arg1);

    void on_i_bar_valueChanged(double arg1);

    void on_j_bar_valueChanged(double arg1);

    void on_x_len_valueChanged(double arg1);

    void on_y_len_valueChanged(double arg1);

    void on_P0_valueChanged(double arg1);

    void on_u0_valueChanged(double arg1);

    void on_v0_valueChanged(double arg1);

    void on_n_lines_valueChanged(double arg1);

    void on_mast_x_valueChanged(double arg1);

    void on_mast_y_valueChanged(double arg1);

    void on_mast_r_valueChanged(double arg1);

private:
    Ui::MainWindow *ui;
    double time;
    double rendered_time;
};

#endif // MAINWINDOW_H

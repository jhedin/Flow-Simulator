#include "ui_vals.h"

ui_vals::ui_vals(QObject *parent) :
    QObject(parent)
{
    dt = 0.0002;
    frame_rate = 30;

    ur = 80; //*(2.0/3.0);
    ul = 80;
    vb = 0;
    vt = 0;

    anu = 1.0;
    i_bar = 50;
    j_bar = 50;
    x_len = 15;
    y_len = 15;
    P0 = 0;
    u0 = 50;
    v0 = 0;

    streamlines = 50;

    mast_x = 7.5;
    mast_y = 7.5;
    mast_r = 1.5;

    P_max = 90;
    P_min = -40;

}

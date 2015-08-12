#ifndef MESH_H
#define MESH_H

#include <QObject>
#include <vector>
#include <math.h>
#include "boost/numeric/ublas/matrix.hpp"
#include <boost/numeric/ublas/io.hpp>
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/vector_proxy.hpp>
#include <boost/numeric/ublas/lu.hpp>
#include <boost/numeric/ublas/triangular.hpp>
#include <GL/gl.h>
#include <GL/glu.h>

class ui_vals;

typedef struct {
    int ip;
    int jp;
} pixel_str;

class mesh : public QObject
{
    Q_OBJECT
public:
    explicit mesh(QObject *parent = 0);
    ui_vals * user;

    // grabs initialization values from the user and creates/fills P, u, and v from them
    void init();

    void step();

    void draw_contour();

    // fills phi using one of tow assumptions, then draws a contour map of it
    void draw_streamlines_object();
    void draw_streamlines_fluid();
    void draw_mast();
    // place colour dots at the right locations in the glwidget
    void draw_pressure_field();

signals:

public slots:

private:

    boost::numeric::ublas::matrix<double> P; // Pressure Field
    boost::numeric::ublas::matrix<double> u, v; // Velocity Fields in x and y
    double dx, dy; // size of a cell
    boost::numeric::ublas::matrix<double> psi;

    // values for calculating the new state of the flow
    boost::numeric::ublas::matrix<double> u_b, v_b;
    boost::numeric::ublas::matrix<double> D; // Diffusivity Field

    // leaving out perturbations for now, todo: make it so mouse clicks cause perturbations
    double perturb; // velocity perturbation to speed up the development of turbulent effects

    // optimization vars
    double D_max;
    double D_test;
    double beta; // constant based on dt, dx, dy

    double sum_entry, sum_exit;

    double psi_max;
    boost::numeric::ublas::matrix<double> streamline_vals;

    std::vector<pixel_str> mast_pixels;
    std::vector<pixel_str> mast_fill;

    // calculates u_b and v_b
    void explicit_step();

    // updates D, P, u, and v
    // returns D_max; this step optimizes the values
    // dont run it more than 100 times per timestep
    double implicit_step();
    void mesh_mast();
    void mesh_pixel_bounds(pixel_str px);
    // x, y, and r are real coordinates, need to map them into is and js
    void circle_px(double r, double xm, double ym, std::vector<pixel_str>*pixels);
    void row(int x, int y, int width, std::vector<pixel_str>*pixels);
    void circle_fill(double r, double xm, double ym, std::vector<pixel_str>*pixels);





};

#endif // MESH_H

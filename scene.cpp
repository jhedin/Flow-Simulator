#include "scene.h"
#include "ui_vals.h"
#include "mesh.h"
#include "mast.h"
#include "sail.h"

scene::scene(QObject *parent) :
    QObject(parent)
{
    mt = new mast();
    ms = new mesh();
    s = new sail();
}

void scene::step()
{
    ms->step();
}

void scene::init()
{
    ms->init();
}

void scene::render()
{
    ms->draw_pressure_field();
    ms->draw_streamlines_object();
    ms->draw_contour();
    ms->draw_mast();
}

#ifndef SCENE_H
#define SCENE_H

#include <QObject>

class mesh;
class sail;
class mast;
class ui_vals;

class scene : public QObject
{
    Q_OBJECT
public:
    explicit scene(QObject *parent = 0);

    ui_vals * user;

    mesh * ms;
    sail * s;
    mast * mt;

    void step();
    void init();
    void render();

signals: 

};

#endif // SCENE_H

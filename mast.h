#ifndef MAST_H
#define MAST_H

#include <QObject>
class ui_vals;

class mast : public QObject
{
    Q_OBJECT
public:
    explicit mast(QObject *parent = 0);
    ui_vals * user;

    void render();


signals:

public slots:

};

#endif // MAST_H

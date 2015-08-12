#ifndef SAIL_H
#define SAIL_H

#include <QObject>

class ui_vals;

class sail : public QObject
{
    Q_OBJECT
public:
    explicit sail(QObject *parent = 0);
    ui_vals * user;

signals:

public slots:

};

#endif // SAIL_H

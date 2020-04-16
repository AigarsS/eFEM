#ifndef POINTLOAD_H
#define POINTLOAD_H

#include <string>

class PointLoad
{
public:
    std::string name;
    int coordX;
    double value;

    PointLoad(std::string name, int coordX, double value);
    ~PointLoad();
};

#endif
#ifndef POINTLOAD_H
#define POINTLOAD_H

#include <iostream>

class PointLoad
{
    private:
        std::string name;
        int coordX;
        double value;

    public:
        PointLoad(std::string name, int coordX, double value);
        ~PointLoad();
        int getCoordX();
        double getLoadValue();
};

#endif
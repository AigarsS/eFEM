#ifndef LOAD_H
#define LOAD_H

#include <iostream>

class Load
{
    private:
        int loadType;
        int coordX_0, coordX_1;
        double value_0, value_1;

    public:
        Load(int loadType, int coordX_0, double value_0, int coordX_1, double value_1);
        ~Load();
        int getCoordX0();
        int getCoordX1();
        double getLoadValue();
        int getLoadType();
};

#endif
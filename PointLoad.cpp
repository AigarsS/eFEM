#include "PointLoad.h"


PointLoad::PointLoad(std::string name, int CoordX, double Value){
    this->name = name;
    this->coordX = coordX;
    this->value = value;
};

PointLoad::~PointLoad(){};
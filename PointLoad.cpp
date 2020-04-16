#include "PointLoad.h"


PointLoad::PointLoad(std::string name, int coordX, double value){
    this->name = name;
    this->coordX = coordX;
    this->value = value;
};

PointLoad::~PointLoad(){};

int PointLoad::getCoordX(){
    return this->coordX;
};

double PointLoad::getLoadValue(){
    return this->value;
};
#include "Load.h"


Load::Load(int loadType, int coordX_0, double value_0, int coordX_1, double value_1){
    this->loadType = loadType;
    this->coordX_0 = coordX_0;
    this->coordX_1 = coordX_1;
    this->value_0 = value_0;
    this->value_1 = value_1;
};

Load::~Load(){};

int Load::getCoordX0(){
    return this->coordX_0;
};

int Load::getCoordX1(){
    return this->coordX_1;
};

double Load::getLoadValue(){
    return this->value_0;
};

int Load::getLoadType(){
    return this->loadType;
};
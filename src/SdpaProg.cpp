//
//  sdpa.cpp
//  Gravity
//  interface for SDPA
//  Created by Guanglei Wang on 11/7/17.
//
//
#include "SdpaProg.h"


SdpaProgram::SdpaProgram(){
    _sdpa_model = new SDPA();
}


SdpaProgram::~SdpaProgram() {
    delete _sdpa_model;
}


void SdpaProgram::read_model(char* filename_input, char* filename_param){
    cout << "SDPA solver" << endl;
    _sdpa_model->readInput(filename_input);
    _sdpa_model->readInput(filename_param);
    _sdpa_model->setDisplay(); // by default, stdout
    _sdpa_model->initializeSolve();
    _sdpa_model->solve();
    _sdpa_model->getDisplay();
    cout << "SDPA solver" << endl;
    fprintf(stdout, "xVec = \n");
    _sdpa_model->printResultXMat();
}


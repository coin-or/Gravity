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
    cout << filename_input << endl;
    FILE* file_output;
    file_output= fopen("../../myresult", "w");
    _sdpa_model->readParameter(filename_param);
    _sdpa_model->readInput(filename_input, file_output);

    _sdpa_model->setDisplay(); // by default, stdout
    _sdpa_model->initializeSolve();
    _sdpa_model->solve();
    _sdpa_model->getDisplay();
    cout << "SDPA outputs solution with" << endl;
    fprintf(stdout, "XMat = \n");
    _sdpa_model->printResultXMat();
    fprintf(stdout, "xVect = \n");
    _sdpa_model->printResultXVec();
    fprintf(stdout, "YMat = \n");
    _sdpa_model->printResultYMat();

    fclose(file_output);
}


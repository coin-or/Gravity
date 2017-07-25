//
//  Sdpa.h
//  Gravity
//
//  Created by Guanglei Wang on 11/7/17.
//
//

#ifndef SdpaProgram_h
#define SdpaProgram_h

#include <stdio.h>
#include <assert.h>
#ifdef USE_SDPA
#include <sdpa_call.h>     /* SDPA callable library interface */
#endif
#include <gravity/model.h>

class SdpaProgram {
public:
   // Model* _model;
    SdpaProgram();
    ~SdpaProgram();
    void read_model(char* filename_input, char* filename_param);
    
private:
    SDPA* _sdpa_model;
};

#endif /* sdpa_h */

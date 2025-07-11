#include <adf.h>
#include "aie_api/aie.hpp"
#include "aie_api/aie_adf.hpp"
#include "aie_api/utils.hpp"

#ifndef FUNCTION_KERNELS_H
#define FUNCTION_KERNELS_H

#define __FSIM__

#define V_SIZE 32
#define NUM_BUNCHES 7


using namespace adf;

// common hyperparameters
static const int16 N_PASS           = 16;
static const int16 MIN_MASS         = 60;
static const int16 MAX_MASS         = 100;
static const float PI               = 3.1415926535;
static const float MINDELTAR2       = 0.5 * 0.5;

// physics constants                
static const float PI_MASS          = 0.13957039;
static const float GAMMA_MASS       = 0.;
                                    
// hadron hyperparameters           
static const int16 PI_HAD           = 720;

static const float PT_CONV_HAD      =  0.25;
static const float ANG_CONV_HAD     =  PI / 720;
static const float ANG_CONV2_HAD    =  (PI / 720) * (PI / 720);

static const int16 MIN_PT_HAD       =  100; // 25 GeV
static const float MAXISO_HAD       =  0.3;
static const int32 MINDR2_HAD       =  0; // 0
static const int32 MAXDR2_HAD       =  13131; // 0.5 * 0.5                                   

// egamma hyperparameters
static const float PT_CONV_EGAMMA   = 0.03125;
static const float ANG_CONV_EGAMMA  = PI / 4096;
static const float ANG_CONV2_EGAMMA = (PI / 4096) * (PI / 4096);

static const int16 MIN_PT_EGAMMA    = 640; // 20 GeV
static const float MAXISO_EGAMMA    = 0.3;

/*
Attention! The following two hyperparameters are
calculateed using ANG_CONV2_HAD
*/
static const int32 MINDR2_EGAMMA    = 21; // 0.02 * 0.02
static const int32 MAXDR2_EGAMMA    = 13131; // 0.5 * 0.5


void wPiGamma(input_stream<int16> * __restrict in0, input_stream<int16> * __restrict in1, output_stream<int32> * __restrict out);

#endif
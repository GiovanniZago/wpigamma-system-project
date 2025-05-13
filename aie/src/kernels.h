#include <adf.h>
#include "aie_api/aie.hpp"
#include "aie_api/aie_adf.hpp"
#include "aie_api/utils.hpp"

#ifndef FUNCTION_KERNELS_H
#define FUNCTION_KERNELS_H

#define V_SIZE 32
#define NUM_BUNCHES 7


using namespace adf;

// hyperparameters
static const int16 N_PASS = 16;
static const int16 MIN_PT_HAD = 800; // 25 GeV
static const int16 MIN_PT_EGAMMA = 640; //  20 GeV
static const int16 PI = 4096;
static const int16 MPI = -4096;
static const int16 TWOPI = 8192;
static const int16 MTWOPI = -8192;
static const int32 MINDR2 = 0;
static const int32 MAXDR2 = 424972;
static const int32 MINDELTAR2 = 424972;
static const int32 MAXISO_PI = 10;
static const int32 MAXISO_EGAMMA = 10;
static const int16 MIN_MASS = 60;
static const int16 MAX_MASS = 100;


// conversion factors
static const float PT_CONV = 0.03125;
static const float ANG_CONV = 3.1415926535 / PI;
static const float ANG_CONV2 = (3.1415926535 / PI) * (3.1415926535 / PI);

// physics constants
static const float PI_MASS = 0.13957039;
static const float GAMMA_MASS = 0.;

void wPiGamma(input_stream<int16> * __restrict in0, input_stream<int16> * __restrict in1, output_stream<int16> * __restrict out);

#endif
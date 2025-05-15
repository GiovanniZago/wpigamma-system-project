#include <adf.h>
#include "kernels.h"

using namespace adf;

class WPiGammaGraph : public graph 
{
    private:
        kernel wPiGamma_k;

    public:
        input_plio in0;
        input_plio in1;
        output_plio out;

        WPiGammaGraph() 
        {
            wPiGamma_k = kernel::create(wPiGamma);

            in0 = input_plio::create("in0", plio_32_bits, "foo0.txt", 360);
            in1 = input_plio::create("in1", plio_32_bits, "foo1.txt", 360);
            out = output_plio::create("out", plio_32_bits, "out.txt", 360);

            // PL inputs
            connect<stream>(in0.out[0], wPiGamma_k.in[0]);
            connect<stream>(in1.out[0], wPiGamma_k.in[1]);

            // PL outputs
            connect<stream>(wPiGamma_k.out[0], out.in[0]);

            // sources and runtime ratios
            source(wPiGamma_k) = "kernels.cpp";
            runtime<ratio>(wPiGamma_k) = 1;
        }
};
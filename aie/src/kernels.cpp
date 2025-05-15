#include "kernels.h"

void wPiGamma(input_stream<int16> * __restrict in0, input_stream<int16> * __restrict in1, output_stream<int32> * __restrict out)
{
    // assume events of 32 candidates
    aie::vector<int16, V_SIZE> pts_had  = { aie::broadcast<int16, V_SIZE>(0) };
    aie::vector<int16, V_SIZE> etas_had = { aie::broadcast<int16, V_SIZE>(0) };
    aie::vector<int16, V_SIZE> phis_had = { aie::broadcast<int16, V_SIZE>(0) }; 
    aie::vector<int16, V_SIZE> pids_had = { aie::broadcast<int16, V_SIZE>(0) };

    aie::vector<int16, V_SIZE> pts_eg  = { aie::broadcast<int16, V_SIZE>(0) };
    aie::vector<int16, V_SIZE> etas_eg = { aie::broadcast<int16, V_SIZE>(0) };
    aie::vector<int16, V_SIZE> phis_eg = { aie::broadcast<int16, V_SIZE>(0) }; 
    aie::vector<int16, V_SIZE> pids_eg = { aie::broadcast<int16, V_SIZE>(0) };

    // control variables
    bool skip_event = false;

    // auxiliari variables
    const aie::vector<int16, V_SIZE> zeros_vector  = aie::zeros<int16, V_SIZE>();
    const aie::vector<int16, V_SIZE> pi_vector     = aie::broadcast<int16, V_SIZE>(PI_HAD);
    const aie::vector<int16, V_SIZE> mpi_vector    = aie::broadcast<int16, V_SIZE>(-PI_HAD);
    const aie::vector<int16, V_SIZE> twopi_vector  = aie::broadcast<int16, V_SIZE>(2 * PI_HAD);
    const aie::vector<int16, V_SIZE> mtwopi_vector = aie::broadcast<int16, V_SIZE>(-2 * PI_HAD);

    // read input variables
    pts_had.insert(0, readincr_v<V_SIZE>(in0));
    etas_had.insert(0, readincr_v<V_SIZE>(in1));
    phis_had.insert(0, readincr_v<V_SIZE>(in0));
    pids_had.insert(0, readincr_v<V_SIZE>(in1));
    pts_eg.insert(0, readincr_v<V_SIZE>(in0));
    etas_eg.insert(0, readincr_v<V_SIZE>(in1));
    phis_eg.insert(0, readincr_v<V_SIZE>(in0));    

    // pt cut and the pid filter for hadrons to identify pions
    aie::mask<V_SIZE> select_had = aie::ge(pts_had, MIN_PT_HAD) & aie::ge(pids_had, (int16) 2) & aie::le(pids_had, (int16) 5);
    unsigned int      count      = select_had.count();

    if (count < 1) skip_event = true;

    // pt cut for egammas
    aie::mask<V_SIZE> select_eg = aie::ge(pts_eg, MIN_PT_EGAMMA);
                      count     = select_eg.count();
                      
    if (count < 1) skip_event = true;

    #if defined (__FSIM__)  
    aie::print(select_had, true, "select_had: ");
    aie::print(select_eg, true, "select_eg: ");
    #endif

    // calculate isolation of pions and egamma
    aie::mask<V_SIZE> iso_had;
    aie::mask<V_SIZE> iso_eg;
    
    for (unsigned int ii = 0; ii < V_SIZE; ii++)
    {
        if (skip_event) continue;

        /*
        Pion isolation
        */
        if (!select_had.test(ii)) continue;

        int16 pt_cur_had = pts_had[ii];
        int16 eta_cur_had = etas_had[ii];
        int16 phi_cur_had = phis_had[ii];

        aie::vector<int16, V_SIZE> d_eta        = aie::sub(eta_cur_had, etas_had);
        
        aie::vector<int16, V_SIZE> d_phi        = aie::sub(phi_cur_had, phis_had);
        aie::vector<int16, V_SIZE> d_phi_ptwopi = aie::add(d_phi, twopi_vector); // d_phi + 2 * pi
        aie::vector<int16, V_SIZE> d_phi_mtwopi = aie::add(d_phi, mtwopi_vector); // d_phi - 2 * pi
        aie::mask<V_SIZE>          is_gt_pi     = aie::gt(d_phi, pi_vector);
        aie::mask<V_SIZE>          is_lt_mpi    = aie::lt(d_phi, mpi_vector);
                                   d_phi        = aie::select(d_phi, d_phi_ptwopi, is_lt_mpi); // select element from d_phi if element is geq of -pi, otherwise from d_phi_ptwopi
                                   d_phi        = aie::select(d_phi, d_phi_mtwopi, is_gt_pi); // select element from d_phi if element is leq of pi, otherwise from d_phi_mtwopi
        
        aie::accum<acc48, V_SIZE>  acc          = aie::mul_square(d_eta); // acc = d_eta ^ 2
                                   acc          = aie::mac_square(acc, d_phi); // acc = acc + d_phi ^ 2
        aie::vector<int32, V_SIZE> dr2          = acc.to_vector<int32>(0); // convert accumulator into vector
        
        aie::mask<V_SIZE>          cone_mask    = aie::ge(dr2, MINDR2_HAD) & aie::le(dr2, MAXDR2_HAD);
        aie::vector<int16, V_SIZE> pt_to_sum    = aie::select(zeros_vector, pts_had, cone_mask);
        int32                      pt_sum       = aie::reduce_add(pt_to_sum);

        if (pt_sum <= MAXISO_HAD * pt_cur_had) iso_had.set(ii);

        #if defined (__FSIM__)  
        printf("\tpt_sum_had: %d\n", pt_sum);
        printf("\tpt_cur_had: %d\n", pt_cur_had);
        aie::print(cone_mask, true, "\t\tcone_mask_had: ");
        #endif

        /*
        Egamma isolation 
        */
        if (!select_eg.test(ii)) continue;

        // convert pt, eta, phi of the egamma candidate into the hadron format in order to allow a comparison
        int16 pt_cur_eg  = static_cast<int16>((pts_eg[ii] * PT_CONV_EGAMMA) / PT_CONV_HAD);
        int16 eta_cur_eg = static_cast<int16>((etas_eg[ii] * ANG_CONV_EGAMMA) / ANG_CONV_HAD);
        int16 phi_cur_eg = static_cast<int16>((phis_eg[ii] * ANG_CONV_EGAMMA) / ANG_CONV_HAD);

                                   d_eta        = aie::sub(eta_cur_eg, etas_had);
                                   
                                   d_phi        = aie::sub(phi_cur_eg, phis_had);
                                   d_phi_ptwopi = aie::add(d_phi, twopi_vector); // d_phi + 2 * pi
                                   d_phi_mtwopi = aie::add(d_phi, mtwopi_vector); // d_phi - 2 * pi
                                   is_gt_pi     = aie::gt(d_phi, pi_vector);
                                   is_lt_mpi    = aie::lt(d_phi, mtwopi_vector);
                                   d_phi        = aie::select(d_phi, d_phi_ptwopi, is_lt_mpi); // select element from d_phi if element is geq of -pi, otherwise from d_phi_ptwopi
                                   d_phi        = aie::select(d_phi, d_phi_mtwopi, is_gt_pi); // select element from d_phi if element is leq of pi, otherwise from d_phi_mtwopi
                                   
                                   acc          = aie::mul_square(d_eta); // acc = d_eta ^ 2
                                   acc          = aie::mac_square(acc, d_phi); // acc = acc + d_phi ^ 2
                                   dr2          = acc.to_vector<int32>(0); // convert accumulator into vector
                                   
                                   cone_mask    = aie::ge(dr2, MINDR2_EGAMMA) & aie::le(dr2, MAXDR2_EGAMMA);
                                   pt_to_sum    = aie::select(zeros_vector, pts_eg, cone_mask);
                                   pt_sum       = aie::reduce_add(pt_to_sum);

        if (pt_sum <= MAXISO_EGAMMA * pt_cur_eg) iso_eg.set(ii);

        #if defined (__FSIM__)  
        printf("\tpt_sum_eg: %d\n", pt_sum);
        printf("\tpt_cur_eg: %d\n", pt_cur_eg);
        aie::print(cone_mask, true, "\t\tcone_mask_eg: ");
        #endif
    }

    #if defined (__FSIM__)  
    aie::print(iso_had, true, "iso_had: ");
    aie::print(iso_eg, true, "iso_eg: ");
    #endif
    
    // loop over hadrons
    for (unsigned int ii = 0; ii < V_SIZE; ii++)
    {
        if (!select_had.test(ii)) continue;
        if (!iso_had.test(ii)) continue;

        // loop over egammas
        for (unsigned int jj = 0; jj < V_SIZE; jj++)
        {
            if (!select_eg.test(jj)) continue;
            if (!iso_eg.test(jj)) continue;

            // check deltaR
            int16 pt_cur_had   = pts_had[ii];
            int16 eta_cur_had  = etas_had[ii];
            int16 phi_cur_had  = phis_had[ii];

            int16 pt_cur_eg    = pts_eg[jj];
            int16 eta_cur_eg   = etas_eg[jj];
            int16 phi_cur_eg   = phis_eg[jj];

            float d_eta_scalar = eta_cur_had * ANG_CONV_HAD - eta_cur_eg * ANG_CONV_EGAMMA;
            float d_phi_scalar = phi_cur_had * ANG_CONV_HAD - phi_cur_eg * ANG_CONV_EGAMMA;
                  d_phi_scalar = (d_phi_scalar <= PI) ? ((d_phi_scalar >= -PI) ? d_phi_scalar : d_phi_scalar + 2 * PI) : d_phi_scalar - 2 * PI;

            float dr2_scalar   = d_eta_scalar * d_eta_scalar + d_phi_scalar * d_phi_scalar;
            if (dr2_scalar < MINDELTAR2) continue;

            // calculate invariant mass
            float px_had       = pt_cur_had * PT_CONV_HAD * aie::cos(phi_cur_had * ANG_CONV_HAD);
            float py_had       = pt_cur_had * PT_CONV_HAD * aie::sin(phi_cur_had * ANG_CONV_HAD);
            float x            = eta_cur_had * ANG_CONV_HAD;
            float sinh         = x + ((x * x * x) / 6);
            float pz_had       = pt_cur_had * PT_CONV_HAD * sinh;
            float e_had        = aie::sqrt(px_had * px_had + py_had * py_had + pz_had * pz_had + PI_MASS * PI_MASS);

            float px_eg        = pt_cur_eg * PT_CONV_EGAMMA * aie::cos(phi_cur_eg * ANG_CONV_EGAMMA);
            float py_eg        = pt_cur_eg * PT_CONV_EGAMMA * aie::sin(phi_cur_eg * ANG_CONV_EGAMMA);
                  x            = eta_cur_eg * ANG_CONV_EGAMMA;
                  sinh         = x + ((x * x * x) / 6);
            float pz_eg        = pt_cur_eg * PT_CONV_EGAMMA * sinh;
            float e_eg         = aie::sqrt(px_eg * px_eg + py_eg * py_eg + pz_eg * pz_eg + GAMMA_MASS * GAMMA_MASS);

            float px_tot       = px_had + px_eg;
            float py_tot       = py_had + py_eg;
            float pz_tot       = pz_had + pz_eg;
            float e_tot        = e_had + e_eg;

            float invariant_mass = aie::sqrt(e_tot * e_tot - px_tot * px_tot - py_tot * py_tot - pz_tot * pz_tot);
            if ((invariant_mass < MIN_MASS) | (invariant_mass > MAX_MASS)) continue;

            writeincr(out, ii);
            writeincr(out, jj);
            writeincr(out, invariant_mass);

            break;
        }

        break;
    }
}
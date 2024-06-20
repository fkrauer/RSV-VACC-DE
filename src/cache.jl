# Prep caches

function make_cache_fit(n_ages::Int64, n_levels::Int64, n_strata::Int64)

    n_states = 8 #9 # N of states with demographic changes
    return (
        dualcache(zeros(n_ages)), #n
        #dualcache(zeros(n_ages, n_ages)), #c
        dualcache(zeros(n_ages, n_levels, n_states, n_strata)), #ageing_out
        dualcache(zeros(n_ages, n_levels, n_states, n_strata)), # ageing_in
        dualcache(zeros(n_ages, n_levels, n_states, n_strata)), # births
        dualcache(zeros(n_ages, n_levels, n_states, n_strata)), # deaths
        dualcache(zeros(n_ages)), # I_tot
        dualcache(zeros(n_ages)), # A_tot
        dualcache(zeros(n_ages, n_ages)), # inf_mat
        dualcache(zeros(n_ages)), # λ
        dualcache(zeros(n_ages, n_levels, n_strata)), # S_to_E
        dualcache(zeros(n_ages, n_levels, n_strata)), # E_to_I
        dualcache(zeros(n_ages, n_levels, n_strata)), # E_to_A
        dualcache(zeros(n_ages, n_levels, n_strata)), # I_to_R1
        dualcache(zeros(n_ages, n_levels, n_strata)), # A_to_R1
        dualcache(zeros(n_ages, n_levels, n_strata)), # M1_to_M2
        dualcache(zeros(n_ages, n_levels, n_strata)), # M2_to_S
        
        dualcache(zeros(n_ages, n_levels)), # R1_to_R2
        dualcache(zeros(n_ages, n_levels)), # waning_out_R2
        dualcache(zeros(n_ages, n_levels)), # R2_to_S

        dualcache(zeros(n_ages, 3)), # mAB 
        dualcache(zeros(n_ages)) # mAB waning


    )
end

function make_cache_vacc(n_ages::Int64, n_levels::Int64, n_strata::Int64)

    n_states = 9 # N of states with demographic changes
    return (
        dualcache(zeros(n_ages)), #n
        #dualcache(zeros(n_ages, n_ages)), #c
        dualcache(zeros(n_ages, n_levels, n_states, n_strata)), #ageing_out
        dualcache(zeros(n_ages, n_levels, n_states, n_strata)), # ageing_in
        dualcache(zeros(n_ages, n_levels, n_states, n_strata)), # births
        dualcache(zeros(n_ages, n_levels, n_states, n_strata)), # deaths
        dualcache(zeros(n_ages)), # I_tot
        dualcache(zeros(n_ages)), # A_tot
        dualcache(zeros(n_ages, n_ages)), # inf_mat
        dualcache(zeros(n_ages)), # λ
        dualcache(zeros(n_ages, n_levels, n_strata)), # S_to_E
        dualcache(zeros(n_ages, n_levels, n_strata)), # E_to_I
        dualcache(zeros(n_ages, n_levels, n_strata)), # E_to_A
        dualcache(zeros(n_ages, n_levels, n_strata)), # I_to_R1
        dualcache(zeros(n_ages, n_levels, n_strata)), # A_to_R1
        dualcache(zeros(n_ages, n_levels, n_strata)), # M1_to_M2
        dualcache(zeros(n_ages, n_levels, n_strata)), # M2_to_S

        dualcache(zeros(n_ages, n_levels, n_strata)), # R1_to_R2 *
        dualcache(zeros(n_ages, n_levels, n_strata)), # waning_out_R2 *
        dualcache(zeros(n_ages, n_levels, n_strata)), # R2_to_S *

        dualcache(zeros(n_ages, n_levels, n_states)), # prop_preg
        dualcache(zeros(n_ages, n_levels, n_states)), # T3_in
        dualcache(zeros(n_ages, n_levels, n_states, n_strata)), # T3_out
        dualcache(zeros(n_ages, n_levels, n_strata)), # hospitalization

        dualcache(zeros(n_ages, n_levels, n_states)), # mAB *
        dualcache(zeros(n_ages, n_levels, n_strata)), # W_to_E

        dualcache(zeros(n_ages, n_levels, n_strata)), # TX_waning1
        dualcache(zeros(n_ages, n_levels, n_strata)), # TX_waning2
        dualcache(zeros(n_ages, n_levels, n_strata)), # TX_waning2_

        dualcache(zeros(n_ages, n_levels)), # vacc_S_OA
        dualcache(zeros(n_ages, n_levels)), #vacc_S_PED
        dualcache(zeros(n_ages, n_levels)), # vacc_S_MAT

        dualcache(zeros(n_ages, n_levels)), # cumul_P
        dualcache(zeros(n_ages, n_levels, n_strata)) # cumul_V

    )
end



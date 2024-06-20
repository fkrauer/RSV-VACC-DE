function load_cb_intervention(θ, 
                            iₚₐ, 
                            iₚₗ,
                            iᵥₗ,
                            iᵥₐ;
                            vₙ::String="none",
                            vₛ::Bool=true, 
                            simyears::Int64,
                            modelfunc::String,
                            switch::Int64
                            )
    
    modelfunc in ("MSEIARS5_fit", "MSEIARS5_vacc") || throw(ArgumentError("model function name incorrectly specified"))

    # start Interventions ----------------------------------------

    coverage_p = θ.pₚ
    coverage_v = θ.v

    if !iszero(coverage_p)

        # start and stop of season ("switch")
        tstarts_p_switch = collect(range(θ.tₚ, step=365.0, length=simyears))
        tstops_p_switch = collect(range(θ.tₚ + θ.dₚ, step=365.0, length=simyears))
        # time points of mAB treatments ("move")
        tstarts_p_move = sort(reduce(vcat, collect.(range.(collect(range(θ.tₚ, step=θ.interval_p, length=Int(θ.nₚ))), step=365.0, length=simyears))))

        # Flip the switch on
        function affect_start_p_switch!(integrator)
            integrator.u[1, 1, switch, 2] = 1.0
        end
    
        # Flip the switch off
        function affect_stop_p_switch!(integrator)
            integrator.u[1, 1, switch, 2] = 0.0
        end

        # Move individuals to the mAB treatment arm/stratum
        function affect_start_p_move!(integrator)
    
            # Keep track of number of individuals to immunise
            mAB_M1 = integrator.u[iₚₐ, iₚₗ, 1, 1] * coverage_p
            mAB_M2 = integrator.u[iₚₐ, iₚₗ, 2, 1] * coverage_p
            mAB_S = integrator.u[iₚₐ, iₚₗ, 3, 1] * coverage_p    

            # keep track of number of doses (D) (only in the vacc model)
            if modelfunc == "MSEIARS5_vacc"
    
                ## beneficial doses
                integrator.u[iₚₐ, iₚₗ, 14, 2] += integrator.u[iₚₐ, iₚₗ, 1, 1] * coverage_p # M1 untreated arm
                integrator.u[iₚₐ, iₚₗ, 14, 2] += integrator.u[iₚₐ, iₚₗ, 2, 1] * coverage_p # M2 untreated arm
                integrator.u[iₚₐ, iₚₗ, 14, 2] += integrator.u[iₚₐ, iₚₗ, 3, 1] * coverage_p # S untreated arm
                ## wasted doses in the treated arm
                integrator.u[iₚₐ, iₚₗ, 14, 2] += integrator.u[iₚₐ, iₚₗ, 1, 2] * coverage_p # M1 treated arm
                integrator.u[iₚₐ, iₚₗ, 14, 2] += integrator.u[iₚₐ, iₚₗ, 2, 2] * coverage_p # M2 treated arm
                integrator.u[iₚₐ, iₚₗ, 14, 2] += integrator.u[iₚₐ, iₚₗ, 3, 2] * coverage_p # S treated arm
                integrator.u[iₚₐ, iₚₗ, 14, 2] += integrator.u[iₚₐ, iₚₗ, 4, 2] * coverage_p # E treated arm
                #integrator.u[iₚₐ, iₚₗ, 14, 2] += integrator.u[iₚₐ, iₚₗ, 5, 2] * coverage_p # I treated arm --> visibly sick people are never immunised
                integrator.u[iₚₐ, iₚₗ, 14, 2] += integrator.u[iₚₐ, iₚₗ, 6, 2] * coverage_p # A treated arm
                integrator.u[iₚₐ, iₚₗ, 14, 2] += integrator.u[iₚₐ, iₚₗ, 7, 2] * coverage_p # R1 treated arm
                integrator.u[iₚₐ, iₚₗ, 14, 2] += integrator.u[iₚₐ, iₚₗ, 8, 2] * coverage_p # R2 treated arm
                integrator.u[iₚₐ, iₚₗ, 14, 2] += integrator.u[iₚₐ, iₚₗ, 9, 2] * coverage_p # W treated arm
                ## wasted doses in the untreated arm
                integrator.u[iₚₐ, iₚₗ, 14, 2] += integrator.u[iₚₐ, iₚₗ, 4, 1] * coverage_p # E untreated arm
                #integrator.u[iₚₐ, iₚₗ, 14, 2] += integrator.u[iₚₐ, iₚₗ, 5, 1] * coverage_p # I untreated arm --> visibly sick people are never immunised
                integrator.u[iₚₐ, iₚₗ, 14, 2] += integrator.u[iₚₐ, iₚₗ, 6, 1] * coverage_p # A untreated arm
                integrator.u[iₚₐ, iₚₗ, 14, 2] += integrator.u[iₚₐ, iₚₗ, 7, 1] * coverage_p # R1 untreated arm
                integrator.u[iₚₐ, iₚₗ, 14, 2] += integrator.u[iₚₐ, iₚₗ, 8, 1] * coverage_p # R2 untreated arm


            end

            # Immunise with mAB: untreated M1, M2 and S to treated M1, M2, S
            integrator.u[iₚₐ, iₚₗ, 1, 1] -= mAB_M1 # Out of M1 in untreated arm
            integrator.u[iₚₐ, iₚₗ, 1, 2] += mAB_M1 # Into M1 in mAB treatment arm

            integrator.u[iₚₐ, iₚₗ, 2, 1] -= mAB_M2 # Out of M2 in untreated arm
            integrator.u[iₚₐ, iₚₗ, 2, 2] += mAB_M2 # Into M2 in mAB treatment arm

            integrator.u[iₚₐ, iₚₗ, 3, 1] -= mAB_S # Out of S in untreated arm
            integrator.u[iₚₐ, iₚₗ, 3, 2] += mAB_S # Into S in mAB treatment arm

        end

        # keep track of immunised individuals (not the same as n of Doses, as Palivizumab can be up to 5 doses per individual)
        function affect_p_countTX!(integrator)

            # beneficial doses
            integrator.u[iₚₐ, iₚₗ, 15, 2] += integrator.u[iₚₐ, iₚₗ, 1, 1] * coverage_p # M1 untreated arm
            integrator.u[iₚₐ, iₚₗ, 15, 2] += integrator.u[iₚₐ, iₚₗ, 2, 1] * coverage_p # M2 untreated arm
            integrator.u[iₚₐ, iₚₗ, 15, 2] += integrator.u[iₚₐ, iₚₗ, 3, 1] * coverage_p # S untreated arm
            # wasted doses 
            integrator.u[iₚₐ, iₚₗ, 15, 2] += integrator.u[iₚₐ, iₚₗ, 4, 1] * coverage_p # E untreated arm
            #integrator.u[iₚₐ, iₚₗ, 15, 2] += integrator.u[iₚₐ, iₚₗ, 5, 1] * coverage_p # I untreated arm --> visibly sick people are never immunised
            integrator.u[iₚₐ, iₚₗ, 15, 2] += integrator.u[iₚₐ, iₚₗ, 6, 1] * coverage_p # A untreated arm
            integrator.u[iₚₐ, iₚₗ, 15, 2] += integrator.u[iₚₐ, iₚₗ, 7, 1] * coverage_p # R1 untreated arm
            integrator.u[iₚₐ, iₚₗ, 15, 2] += integrator.u[iₚₐ, iₚₗ, 8, 1] * coverage_p # R2 untreated arm

            integrator.u[iₚₐ, iₚₗ, 15, 2] += integrator.u[iₚₐ, iₚₗ, 1, 2] * coverage_p # M1 treated arm
            integrator.u[iₚₐ, iₚₗ, 15, 2] += integrator.u[iₚₐ, iₚₗ, 2, 2] * coverage_p # M2 treated arm
            integrator.u[iₚₐ, iₚₗ, 15, 2] += integrator.u[iₚₐ, iₚₗ, 3, 2] * coverage_p # S treated arm
            integrator.u[iₚₐ, iₚₗ, 15, 2] += integrator.u[iₚₐ, iₚₗ, 4, 2] * coverage_p # E treated arm
            #integrator.u[iₚₐ, iₚₗ, 15, 2] += integrator.u[iₚₐ, iₚₗ, 5, 2] * coverage_p # I treated arm --> visibly sick people are never immunised
            integrator.u[iₚₐ, iₚₗ, 15, 2] += integrator.u[iₚₐ, iₚₗ, 6, 2] * coverage_p # A treated arm
            integrator.u[iₚₐ, iₚₗ, 15, 2] += integrator.u[iₚₐ, iₚₗ, 7, 2] * coverage_p # R1 treated arm
            integrator.u[iₚₐ, iₚₗ, 15, 2] += integrator.u[iₚₐ, iₚₗ, 8, 2] * coverage_p # R2 treated arm
            integrator.u[iₚₐ, iₚₗ, 15, 2] += integrator.u[iₚₐ, iₚₗ, 9, 2] * coverage_p # W treated arm

        end

    
        cbstart_p_switch = PresetTimeCallback(tstarts_p_switch, affect_start_p_switch!, save_positions=falses(length(tstarts_p_switch)))
        cbstop_p_switch = PresetTimeCallback(tstops_p_switch, affect_stop_p_switch!, save_positions=falses(length(tstops_p_switch)))        
        cbstart_p_move = PresetTimeCallback(tstarts_p_move, affect_start_p_move!, save_positions=falses(length(tstarts_p_move)))
        if modelfunc == "MSEIARS5_vacc"
            cbstart_p_countTX = PresetTimeCallback(tstarts_p_switch, affect_p_countTX!, save_positions=falses(length(tstarts_p_switch)))
        else
            cbstart_p_countTX = nothing
        end

    end

    if !iszero(coverage_v)
    
        if vₛ # seasonal

            interval = θ.tiᵥ # interval between seasonal vaccinations
            start_v = θ.tᵥ # start of season
            duration_v = θ.dᵥ # duration of season

            tstarts_v_switch = collect(range(start_v, step=365.0, length=simyears)) # start times of seasonal window ("switch")
            tstops_v_switch = collect(range(start_v + duration_v, step=365.0, length=simyears)) # stop times of seasonal window ("switch")
            tstarts_v_move = collect(range(start_v, step = interval * 365.0, stop = simyears * 365.0)) # start time(s) of vaccinating individuals ("move")

            function affect_start_v_switch!(integrator)
    
                # Flip the switch on
                integrator.u[1, 1, switch, 3] = 1.0
    
            end

            function affect_start_v_move!(integrator)
    
                if vₙ == "maternal" 
                    from_stratum = 5
                    to_stratum = 6
                else 
                    from_stratum = 1
                    to_stratum = 3
                end
        
                # Keep track of number of individuals to immunise
                vacc_S_ = integrator.u[iᵥₐ, iᵥₗ, 3, from_stratum] * coverage_v

                # keep track of number of doses 
    
                ## Beneficial doses
                integrator.u[iᵥₐ, iᵥₗ, 14, to_stratum] += vacc_S_
  
                ## Wasted doses in the already treated arm
                integrator.u[iᵥₐ, iᵥₗ, 14, to_stratum] += integrator.u[iᵥₐ, iᵥₗ, 1, to_stratum] * coverage_v #M1 treated arm
                integrator.u[iᵥₐ, iᵥₗ, 14, to_stratum] += integrator.u[iᵥₐ, iᵥₗ, 2, to_stratum] * coverage_v #M2 treated arm
                integrator.u[iᵥₐ, iᵥₗ, 14, to_stratum] += integrator.u[iᵥₐ, iᵥₗ, 3, to_stratum] * coverage_v #S treated arm
                integrator.u[iᵥₐ, iᵥₗ, 14, to_stratum] += integrator.u[iᵥₐ, iᵥₗ, 4, to_stratum] * coverage_v #E treated arm
                #integrator.u[iᵥₐ ,iᵥₗ, 14, to_stratum] += integrator.u[iᵥₐ, iᵥₗ, 5, to_stratum] * coverage_v #I treated arm
                integrator.u[iᵥₐ, iᵥₗ, 14, to_stratum] += integrator.u[iᵥₐ, iᵥₗ, 6, to_stratum] * coverage_v #A treated arm
                integrator.u[iᵥₐ, iᵥₗ, 14, to_stratum] += integrator.u[iᵥₐ, iᵥₗ, 7, to_stratum] * coverage_v #R1 treated arm
                integrator.u[iᵥₐ, iᵥₗ, 14, to_stratum] += integrator.u[iᵥₐ, iᵥₗ, 8, to_stratum] * coverage_v  #R2 treated arm
                integrator.u[iᵥₐ, iᵥₗ, 14, to_stratum] += integrator.u[iᵥₐ, iᵥₗ, 9, to_stratum] * coverage_v  #W treated arm
                # wasted doses in the untreated arm
                integrator.u[iᵥₐ ,iᵥₗ, 14, to_stratum] += integrator.u[iᵥₐ, iᵥₗ, 1, from_stratum] * coverage_v #M1 untreated arm
                integrator.u[iᵥₐ, iᵥₗ, 14, to_stratum] += integrator.u[iᵥₐ, iᵥₗ, 2, from_stratum] * coverage_v #M2 untreated arm
                integrator.u[iᵥₐ, iᵥₗ, 14, to_stratum] += integrator.u[iᵥₐ, iᵥₗ, 4, from_stratum] * coverage_v #E untreated arm
                #integrator.u[iᵥₐ ,iᵥₗ, 14, to_stratum] += integrator.u[iᵥₐ, iᵥₗ, 5, from_stratum] * coverage_v #I untreated arm
                integrator.u[iᵥₐ, iᵥₗ, 14, to_stratum] += integrator.u[iᵥₐ, iᵥₗ, 6, from_stratum] * coverage_v #A untreated arm
                integrator.u[iᵥₐ, iᵥₗ, 14, to_stratum] += integrator.u[iᵥₐ, iᵥₗ, 7, from_stratum] * coverage_v #R1 untreated arm
                integrator.u[iᵥₐ, iᵥₗ, 14, to_stratum] += integrator.u[iᵥₐ, iᵥₗ, 8, from_stratum] * coverage_v  #R2 untreated arm
                
                # keep track of number of immunised individuals (technically it's the same as state D, but for analogy with mAB, we repeat it here)
                
                ## Beneficial doses
                integrator.u[iᵥₐ, iᵥₗ, 15, to_stratum] += vacc_S_

                ## Wasted doses 
                integrator.u[iᵥₐ, iᵥₗ, 15, to_stratum] += integrator.u[iᵥₐ, iᵥₗ, 1, from_stratum] * coverage_v #M1 untreated arm
                integrator.u[iᵥₐ, iᵥₗ, 15, to_stratum] += integrator.u[iᵥₐ, iᵥₗ, 2, from_stratum] * coverage_v #M2 untreated arm
                integrator.u[iᵥₐ, iᵥₗ, 15, to_stratum] += integrator.u[iᵥₐ, iᵥₗ, 4, from_stratum] * coverage_v #E untreated arm
                #integrator.u[iᵥₐ ,iᵥₗ, 15, to_stratum] += integrator.u[iᵥₐ, iᵥₗ, 5, from_stratum] * coverage_v #I untreated arm
                integrator.u[iᵥₐ, iᵥₗ, 15, to_stratum] += integrator.u[iᵥₐ, iᵥₗ, 6, from_stratum] * coverage_v #A untreated arm
                integrator.u[iᵥₐ, iᵥₗ, 15, to_stratum] += integrator.u[iᵥₐ, iᵥₗ, 7, from_stratum] * coverage_v #R1 untreated arm
                integrator.u[iᵥₐ, iᵥₗ, 15, to_stratum] += integrator.u[iᵥₐ, iᵥₗ, 8, from_stratum] * coverage_v  #R2 untreated arm

                integrator.u[iᵥₐ, iᵥₗ, 15, to_stratum] += integrator.u[iᵥₐ, iᵥₗ, 1, to_stratum] * coverage_v #M1 treated arm
                integrator.u[iᵥₐ, iᵥₗ, 15, to_stratum] += integrator.u[iᵥₐ, iᵥₗ, 2, to_stratum] * coverage_v #M2 treated arm
                integrator.u[iᵥₐ, iᵥₗ, 15, to_stratum] += integrator.u[iᵥₐ, iᵥₗ, 3, to_stratum] * coverage_v #S treated arm
                integrator.u[iᵥₐ, iᵥₗ, 15, to_stratum] += integrator.u[iᵥₐ, iᵥₗ, 4, to_stratum] * coverage_v #E treated arm
                #integrator.u[iᵥₐ ,iᵥₗ, 15, to_stratum] += integrator.u[iᵥₐ, iᵥₗ, 5, to_stratum] * coverage_v #I treated arm
                integrator.u[iᵥₐ, iᵥₗ, 15, to_stratum] += integrator.u[iᵥₐ, iᵥₗ, 6, to_stratum] * coverage_v #A treated arm
                integrator.u[iᵥₐ, iᵥₗ, 15, to_stratum] += integrator.u[iᵥₐ, iᵥₗ, 7, to_stratum] * coverage_v #R1 treated arm
                integrator.u[iᵥₐ, iᵥₗ, 15, to_stratum] += integrator.u[iᵥₐ, iᵥₗ, 8, to_stratum] * coverage_v  #R2 treated arm
                integrator.u[iᵥₐ, iᵥₗ, 15, to_stratum] += integrator.u[iᵥₐ, iᵥₗ, 9, to_stratum] * coverage_v  #W treated arm

                # Immunise: untreated S to treated S
                integrator.u[iᵥₐ, iᵥₗ, 3, from_stratum] -= vacc_S_
                integrator.u[iᵥₐ, iᵥₗ, 3, to_stratum] += vacc_S_

            end
        
            function affect_stop_v_switch!(integrator)
                # Flip the switch off
                integrator.u[1, 1, switch, 3] = 0.0
            end
    
            cbstart_v_switch = PresetTimeCallback(tstarts_v_switch, affect_start_v_switch!, save_positions=falses(length(tstarts_v_switch)))
            cbstop_v_switch = PresetTimeCallback(tstops_v_switch, affect_stop_v_switch!, save_positions=falses(length(tstops_v_switch)))
            cbstart_v_move = PresetTimeCallback(tstarts_v_move, affect_start_v_move!, save_positions=falses(length(tstarts_v_move)))

        else # all-year
            tstarts_v_switch = [0.0]

            function affect_start_v_cont!(integrator)
    
                # Flip the switch on
                integrator.u[1, 1, switch, 3] = 1.0
            end

            cbstart_v_switch = PresetTimeCallback(tstarts_v_switch, affect_start_v_cont!, save_positions=falses(length(tstarts_v_switch)))
            cbstart_v_move = nothing
            cbstop_v_switch = nothing
        end

    end

    if !iszero(coverage_p) && !iszero(coverage_v) # mAB and vacc
        return CallbackSet(cbstart_p_switch, cbstop_p_switch, cbstart_p_move, cbstart_p_countTX, cbstart_v_switch, cbstart_v_move, cbstop_v_switch)
    elseif !iszero(coverage_p) && iszero(coverage_v) # only mAB
        return CallbackSet(cbstart_p_switch, cbstop_p_switch, cbstart_p_move, cbstart_p_countTX)
    elseif iszero(coverage_p) && !iszero(coverage_v) # vacc only
        return CallbackSet(cbstart_v_switch, cbstart_v_move, cbstop_v_switch) 
    else
        return nothing
    end

end



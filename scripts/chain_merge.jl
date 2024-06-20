using ArgParse
using DrWatson
include(scriptsdir("utils.jl"))


function chain_merge(runs; load_model=false, ignorecommit=false)

    # Load the arguments early so we can checkout the correct version.
    runs_args = map(runs) do rundir
        DrWatson.wload(joinpath(rundir, "args.jld2"))
    end;

    # Checks
    @info "(♻) Comparing model runs for consistency..."

    if length(unique((x["betafunc"] for x in runs_args))) > 1
        println("(⚠) The betafunc is not identical among all runs. Exit")
        exit(1)
    end

    unique_gitcommit = unique((x["gitcommit"] for x in runs_args))
    @info "Checking git commit(s) $(unique_gitcommit)..."
    if length(unique_gitcommit) > 1
        println("(⚠) Runs come from different commits; do you want to continue? [y/N]: ")
        answer = readline()
        if lowercase(answer) != "y"
            println("Exiting...")
            exit(1)
        end

        gitcommit = runs_args[1]["gitcommit"]
        println("We'll check out $(gitcommit) then.")
    end

    if !ignorecommit
        unique_gitpatch = unique((getkey(x, "gitpatch", "") for x in runs_args))
        if length(unique_gitpatch) > 1
            @warn "(⚠) Oh no!!! Runs come come from different patches. Applying the first one..."
        end

        @info "(♻) Checking out $(first(unique_gitcommit))..."
        checkout!(runs_args[1])
        @info "(✓) Checkout completed!"
    end

    # Load the chains.
    chains = map(runs) do rundir
        deserialize(joinpath(rundir, "chain.jls"))
    end;

    # Load the models.
    if load_model
        models = map(runs) do rundir
            deserialize(joinpath(rundir, "model.jls"))
        end;
        model = models[1]
    else
        model = nothing
    end

    # Combine the chains if so specified.
    @warn "Assuming first args/model is representative of all args/models."
    run_args = runs_args[1]
    
    @info "Combining $(length(chains)) chains into a single chain."
    chain = MCMCChains.chainscat(chains...)

    # Drop the warmup.
    nadapts = run_args["nadapts"] 
    chain = chain[nadapts + 1:end] 
    @info "Dropping $(nadapts) adaptation steps; $(length(chain)) samples remaining."

    return chain, run_args, model

end


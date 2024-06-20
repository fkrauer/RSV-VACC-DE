using DrWatson, DataFrames, CSV, DynamicPPL, Serialization, Dates

# Get all available runs
available_runs = filter(
    isdir ∘ Base.Fix1(projectdir, "output"),
    readdir(projectdir("output"))
)

# Get args
run_args = collect_results(projectdir("output"), subfolders=true, valid_filetypes=["args.jld2"])
@. run_args.path = replace(run_args.path, "args.jld2" => "")
run_args = select(run_args, sort(names(run_args)))

# Get Git commit ID
gitcommit = map(available_runs) do rundir
    try
        DrWatson.wload(joinpath(projectdir("output"), rundir, "args.jld2"))["gitcommit"]
    catch
        @warn "file not available, skipping"
    end
end
gitcommit = DataFrame(gitcommit = gitcommit)
gitcommit.name = available_runs

# Calculate runtime
starttime = map(available_runs) do rundir
        mtime(joinpath(projectdir("output"), rundir, "args.jld2"))
        #Dates.unix2datetime(mtime(joinpath(projectdir("output"), rundir, "args.jld2")))
end

endtime = map(available_runs) do rundir
        mtime(joinpath(projectdir("output"), rundir, "trace.jls"))
        #Dates.unix2datetime(mtime(joinpath(projectdir("output"), rundir, "trace.jls")))
end

runtime = round.((endtime .- starttime) ./ 3600, digits=1)
runtime[runtime .< 0] .= 0.0
runtime = hcat(runtime, available_runs)
runtime = DataFrame(runtime, [:runtime, :name])

# Get ess from results
ess_min = collect_results(projectdir("output"), subfolders=true, valid_filetypes=["result.jld2"])
@. ess_min.path = replace(ess_min.path, "result.jld2" => "")

# Get initial parameter states sampled from prior
runs = filter(
    isdir ∘ Base.Fix1(projectdir, "output"),
    readdir(projectdir("output"))
)

inits = map(runs) do rundir
    path = joinpath(projectdir("output"), rundir, "theta_init.jls")
    if isfile(path)
        var_info = deserialize(path)
        DynamicPPL.invlink!(var_info, SampleFromPrior())
        foo = DynamicPPL.getval(var_info, keys(var_info))
        DataFrame(theta_init = NamedTuple{Symbol.(Tuple(keys(var_info)))}(foo), name = rundir)
    else
        DataFrame(theta_init = "", name = rundir)
    end
end
inits = reduce(vcat, inits, cols = :union)

# Combine
out = outerjoin(run_args, ess_min, on=:path)
out = outerjoin(runtime, out, on=:name)
out = outerjoin(out, gitcommit, on=:name)
out = outerjoin(out, inits, on=:name)
out = select!(out, Not(:path))

CSV.write("output/runs.txt", out, transform=(col, val) -> something(val, missing)) 
CSV.write("output/runs.csv", out, transform=(col, val) -> something(val, missing)) 
# Unicode symbols don't render if written to .csv or .xlsx, unfortunately

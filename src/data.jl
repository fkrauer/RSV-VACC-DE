function load_data(; n_agegroups=25)

    # Epi data
    epidata = DataFrame(CSV.File(joinpath(dirname(pathof(RSVVACCDE)), "..", "data", string("epidata_", n_agegroups, ".csv"))))
    sort!(epidata, [:groupno])

    # 1. Time series lab data (AGI)
    agi = DataFrame(CSV.File(joinpath(dirname(pathof(RSVVACCDE)), "..", "data/data_AGI_ts_age.csv")))
    agi_mat = agi[:, 2:end] |> Tables.matrix
    ts_obsindex = findall(!ismissing, vec(agi_mat))

    # age-stratified weekly data
    data_ts_age = zeros(Int64, length(ts_obsindex), 1)
    data_ts_age .= vec(agi_mat)[ts_obsindex]
    times_obsindex = agi.time 

    # total weekly data for fitting
    data_ts = map(sum, skipmissing.(eachrow(agi_mat)))
    data_ts = reshape(data_ts, length(data_ts),1)

    # proportional case load by age group for fitting
    data_ts_prop = DataFrame(CSV.File(joinpath(dirname(pathof(RSVVACCDE)), "..", "data/data_AGI_prop.csv")))
    data_ts_prop = data_ts_prop.num
    N_ts = sum(data_ts_prop)

    # 2. Hospitalization data (TK)
    tk = DataFrame(CSV.File(joinpath(dirname(pathof(RSVVACCDE)), "..", "data", string("data_TK", n_agegroups, "_hosp_quarter.csv"))))
    tk = tk[:, 2:end] |> Tables.matrix
    hosp_obsindex = findall(!ismissing, vec(tk))
    # age-stratified quarterly data
    data_hosp_age = zeros(Int64, length(hosp_obsindex), 1)
    data_hosp_age .= vec(tk)[hosp_obsindex]
    # total quarterly data for fitting
    data_hosp = map(sum, skipmissing.(eachrow(tk)))
    data_hosp = data_hosp[findall(data_hosp .!= 0.0),:]  
    # proportional case load by age group for fitting
    data_hosp_prop = DataFrame(CSV.File(joinpath(dirname(pathof(RSVVACCDE)), "..", "data", string("data_TK", n_agegroups, "_prop.csv"))))
    data_hosp_prop = data_hosp_prop.num   
    N_hosp = sum(data_hosp_prop)

    # Number of days per quarter, to match daily incidence to quarter hospital data
    dim_q = repeat([91, 91, 92, 91],4) 

    # 3. seroconversion data from Netherlands (Andeweg et al)
    seroconversion = DataFrame(CSV.File(joinpath(dirname(pathof(RSVVACCDE)), "..", "data", "seroconversion.csv")))
    data_noconv = Array{Int}(undef, (11, 1))
    data_noconv .= seroconversion.noconv
    N_noconv = Array{Int}(undef, (11, 1))
    N_noconv .= seroconversion.N

    # 4. contact matrix
    c = CSV.File(joinpath(dirname(pathof(RSVVACCDE)), "..", "data", string("contacts_", n_agegroups, "_raw.csv")), header=false) |> Tables.matrix

    return epidata, data_ts, data_ts_age, data_ts_prop, N_ts, times_obsindex, ts_obsindex, data_hosp, data_hosp_age, data_hosp_prop, N_hosp, hosp_obsindex, dim_q, data_noconv, N_noconv, c
end

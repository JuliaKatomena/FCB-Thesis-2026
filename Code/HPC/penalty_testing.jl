#import Pkg
#using Pkg
#Pkg.activate(@__DIR__)


#cd(@__DIR__)

import Pkg
Pkg.activate(raw"D:/DTU/OR/VS_codes/FCB")
Pkg.status()

using JuMP, Plots, HiGHS, Statistics, LinearAlgebra, LaTeXStrings, DataFrames, CSV, XLSX, StatsPlots, NPZ, UnPack, JLD2, PyCall, Random
include("constants.jl")
include("opt_models.jl")
include("outputs_processing.jl")

outdir = joinpath(@__DIR__, "export_HPC")
mkpath(outdir)  # creates folder if missing (and parents too)

file_list = ["seasonal_scenarios_0125.csv"]

Pw, centers, Prob, nS_all = preprocess_scenario_files(file_list; season="Winter")         # Considering that they all have the same uncertainty

hour = 1
prob = Vector{Float64}(undef, nS_all)
prob = round.(Prob[1][:, hour], digits = 5)     # in this case we prob is only constituted of 1 element

k = 10
fold_size = div(nS_all, k)
Pw_hour = Pw[hour]
n_hours = size(Pw)[1]

fold = 1

# Assign the entire complex NamedTuple type definition to a variable.
BMResultType = @NamedTuple{
    profits::@NamedTuple{pol_in::Matrix{Vector{Float64}}, flex_in::Matrix{Vector{Float64}}, reg_in::Matrix{Vector{Float64}}, pol_out::Matrix{Vector{Float64}}, flex_out::Matrix{Vector{Float64}}, reg_out::Matrix{Vector{Float64}}},
    bids::@NamedTuple{B_pol_in::Matrix{Float64}, B_pol_up_in::Matrix{Float64}, B_pol_down_in::Matrix{Float64}, B_flex_up_in::Matrix{Float64}, B_flex_down_in::Matrix{Float64}, B_pol_out::Matrix{Float64}, B_pol_up_out::Matrix{Float64}, B_pol_down_out::Matrix{Float64}, B_flex_up_out::Matrix{Float64}, B_flex_down_out::Matrix{Float64}, DA_bids_insample::Vector{Any}, DA_bids_outsample::Vector{Any}, B_bids_all_insample::Vector{Any}, B_bids_all_outsample::Vector{Any}},
    aggregates::@NamedTuple{pol_out_avg::Matrix{Float64}, flex_out_avg::Matrix{Float64}, reg_out_avg::Matrix{Float64}, pol_in_sums::Matrix{Float64}, pol_out_sums::Matrix{Float64}, flex_in_sums::Matrix{Float64}, flex_out_sums::Matrix{Float64}, reg_in_sums::Matrix{Float64}, reg_out_sums::Matrix{Float64}, pol_in_avg::Matrix{Float64}, flex_in_avg::Matrix{Float64}, reg_in_avg::Matrix{Float64}, mean_pol_bids_in::Matrix{Float64}, mean_pol_bids_out::Matrix{Float64}, mean_pol_bids_up_in::Matrix{Float64}, mean_pol_bids_up_out::Matrix{Float64}, mean_pol_bids_down_in::Matrix{Float64}, mean_pol_bids_down_out::Matrix{Float64}, mean_flex_bids_up_in::Matrix{Float64}, mean_flex_bids_up_out::Matrix{Float64}, mean_flex_bids_down_in::Matrix{Float64}, mean_flex_bids_down_out::Matrix{Float64}, mean_pol_DA_bids_out::Matrix{Float64}, mean_flex_DA_bids_out::Matrix{Float64}, mean_reg_DA_bids_out::Matrix{Float64}},
    history::@NamedTuple{lambdas::Vector{Any}, history_insample::Vector{Any}, system_imbalance::Vector{Any}},
    split_info::@NamedTuple{in_sample_prob::Vector{Vector{Float64}}, out_sample_prob::Vector{Vector{Float64}}, in_sample_Pw::Vector{Matrix{Float64}}, out_sample_Pw::Vector{Matrix{Float64}}}
}

function extract_vals_from_history(run)
    h = run[:history][:history_insample]
    
    # x[1] w Julii to odpowiednik x[0] w Pythonie (indeksowanie od 1)
    iter_bin    = [x[1].iter for x in h]
    iter_bin_LP = [x[2].iter for x in h]

    obj_bin    = [round(x[1].obj; digits = 2) for x in h]
    obj_bin_LP = [round(x[2].obj; digits = 2) for x in h]

    return iter_bin, iter_bin_LP, obj_bin, obj_bin_LP
end


penalty = 100:100:2000

all_Cf_R_BM   = []
all_Cf_R_GDCA = []


iter_bin_all = []
iter_bin_LP_all = []
obj_bin_all = []
obj_bin_LP_all = []

time_taken = []


results_RI_NI = reserve_RI_NI(NI_up_og, NI_down_og, Cf_Rup_og,  Cf_Rdown_og)
results_RI = reserve_RI(Cf_Rup_og, Cf_Rdown_og)
rf_up, rf_down = results_RI_NI[end].rf_up, results_RI_NI[end].rf_down
q_up, q_down = get_your_qs(results_RI_NI[end], results_RI[end])
Lambda_r_up, Lambda_r_down, fig = merit_order_price_reserve(rf_up, rf_down, Cf_Rup_og, Cf_Rdown_og; dir = "down regulation")
Lambda_R = (Lambda_r_up, Lambda_r_down)
R_bids = [rf_up, rf_down]

res_BM = run_kfold_simulation("Benchmark", Pw, nS_all, prob, 0, 0, fold, fold_size,
                            rf_up, rf_down, NI_up_og, NI_down_og, Cf_Rup_og, Cf_Rdown_og, R_bids, Lambda_R, hour)
push!(all_Cf_R_BM, res_BM)

for p in penalty

    time_1 = time()
    R_GDCA= run_kfold_simulation("TA_bin_GDCA", Pw, nS_all, prob, q_up, q_down, fold, fold_size,
                            rf_up, rf_down, NI_up_og, NI_down_og, Cf_Rup_og, Cf_Rdown_og, R_bids, Lambda_R, hour; history = res_BM.history.history_insample, max_iter = 400, penalty_t = p)
    time_2 = time()
    
    iter_bin, iter_bin_LP, obj_bin, obj_bin_LP = extract_vals_from_history(R_GDCA)

    push!(time_taken, time_2 - time_1)
    push!(iter_bin_all, iter_bin)
    push!(iter_bin_LP_all, iter_bin_LP)
    push!(obj_bin_all, obj_bin)
    push!(obj_bin_LP_all, obj_bin_LP)
    push!(all_Cf_R_GDCA, R_GDCA)

    jobid = get(ENV, "LSB_JOBID", "manual")
    checkpoint = joinpath(outdir, "checkpoint_" * jobid * ".jld2")

    data_partial = Dict(
        "BM" => all_Cf_R_BM,
        "GDCA" => all_Cf_R_GDCA,
        "q_up" => q_up,
        "q_down" => q_down,
        "iter_bin_all" => iter_bin_all,
        "iter_bin_LP_all" => iter_bin_LP_all,
        "obj_bin_all" => obj_bin_all,
        "obj_bin_LP_all" => obj_bin_LP_all,
        "time_taken" => time_taken,
        "penalty" => p,
    )

    tmp = checkpoint * ".tmp"
    @save tmp data_partial
    mv(tmp, checkpoint; force=true)     
        
end


# after the loop

data = Dict(
        "BM" => all_Cf_R_BM,
        "GDCA" => all_Cf_R_GDCA,
        "q_up"=> q_up,
        "q_down" => q_down,
        "iter_bin" => iter_bin_all,
        "iter_bin_LP" => iter_bin_LP_all,
        "obj_bin" => obj_bin_all,
        "obj_bin_LP" => obj_bin_LP_all,
        "time_taken" => time_taken,
        "penalty" => penalty,
)


jobid = get(ENV, "LSB_JOBID", "manual")
@save joinpath(outdir, "penalty_test_" * jobid * ".jld2") data



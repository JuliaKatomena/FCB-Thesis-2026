import Pkg
using Pkg
Pkg.activate(@__DIR__)


cd(@__DIR__)


using JuMP, Plots, HiGHS, Statistics, LinearAlgebra, LaTeXStrings, DataFrames, CSV, XLSX, StatsPlots, NPZ, UnPack, JLD2, PyCall, Random
include("constants.jl")
include("opt_models.jl")
include("outputs_processing.jl")

#outdir = joinpath(@__DIR__, "export_HPC")
#mkpath(outdir)  # creates folder if missing (and parents too)

file_list = ["seasonal_scenarios_0125.csv"]

Pw, centers, Prob, nS_all = preprocess_scenario_files(file_list; season="Spring")         # Considering that they all have the same uncertainty

hour = 1
prob = Vector{Float64}(undef, nS_all)
prob = round.(Prob[1][:, hour], digits = 5)     # in this case we prob is only constituted of 1 element

k = 10
fold_size = div(nS_all, k)
Pw_hour = Pw[hour]
n_hours = size(Pw)[1]

# Assign the entire complex NamedTuple type definition to a variable.
BMResultType = @NamedTuple{
    profits::@NamedTuple{pol_in::Matrix{Vector{Float64}}, flex_in::Matrix{Vector{Float64}}, reg_in::Matrix{Vector{Float64}}, pol_out::Matrix{Vector{Float64}}, flex_out::Matrix{Vector{Float64}}, reg_out::Matrix{Vector{Float64}}},
    bids::@NamedTuple{B_pol_in::Matrix{Float64}, B_pol_up_in::Matrix{Float64}, B_pol_down_in::Matrix{Float64}, B_flex_up_in::Matrix{Float64}, B_flex_down_in::Matrix{Float64}, B_pol_out::Matrix{Float64}, B_pol_up_out::Matrix{Float64}, B_pol_down_out::Matrix{Float64}, B_flex_up_out::Matrix{Float64}, B_flex_down_out::Matrix{Float64}, DA_bids_insample::Vector{Any}, DA_bids_outsample::Vector{Any}, B_bids_all_insample::Vector{Any}, B_bids_all_outsample::Vector{Any}},
    aggregates::@NamedTuple{pol_out_avg::Matrix{Float64}, flex_out_avg::Matrix{Float64}, reg_out_avg::Matrix{Float64}, pol_in_sums::Matrix{Float64}, pol_out_sums::Matrix{Float64}, flex_in_sums::Matrix{Float64}, flex_out_sums::Matrix{Float64}, reg_in_sums::Matrix{Float64}, reg_out_sums::Matrix{Float64}, pol_in_avg::Matrix{Float64}, flex_in_avg::Matrix{Float64}, reg_in_avg::Matrix{Float64}, mean_pol_bids_in::Matrix{Float64}, mean_pol_bids_out::Matrix{Float64}, mean_pol_bids_up_in::Matrix{Float64}, mean_pol_bids_up_out::Matrix{Float64}, mean_pol_bids_down_in::Matrix{Float64}, mean_pol_bids_down_out::Matrix{Float64}, mean_flex_bids_up_in::Matrix{Float64}, mean_flex_bids_up_out::Matrix{Float64}, mean_flex_bids_down_in::Matrix{Float64}, mean_flex_bids_down_out::Matrix{Float64}, mean_pol_DA_bids_out::Matrix{Float64}, mean_flex_DA_bids_out::Matrix{Float64}, mean_reg_DA_bids_out::Matrix{Float64}},
    history::@NamedTuple{lambdas::Vector{Any}, history_insample::Vector{Any}, system_imbalance::Vector{Any}},
    split_info::@NamedTuple{in_sample_prob::Vector{Vector{Float64}}, out_sample_prob::Vector{Vector{Float64}}, in_sample_Pw::Vector{Matrix{Float64}}, out_sample_Pw::Vector{Matrix{Float64}}}
}

function extract_iters_from_run(run)
    h = run[:history][:history_insample]
    
    # x[1] w Julii to odpowiednik x[0] w Pythonie (indeksowanie od 1)
    iter_bin    = [x[1].iter for x in h]
    iter_bin_LP = [x[2].iter for x in h]
    
    return iter_bin, iter_bin_LP
end


cost_Rup_SCGT = 5.0:12.0 #5.6.7.   8.9.  10.11.  12
cost_Rup_CCGT = 4.0:11.0 #4.5.6.  7.8.  9.10.  11
cost_Rdown_SCGT = 2.0:9.0 #2.3.4. 5.6. 7.8. 9
cost_Rdown_CCGT = 1.0:8.0 #1.2.3. 4.5. 6.7. 8

penalty = 1000
size_Cf_R = length(cost_Rup_SCGT)


all_Cf_R_BM   = []   
all_Cf_R_GDCA = []
iter_bin_all = []
iter_bin_LP_all = []

Cf_Rup = Vector{Vector{Float64}}(undef, size_Cf_R)
Cf_Rdown = Vector{Vector{Float64}}(undef, size_Cf_R)
q_up = Vector{Float64}(undef, size_Cf_R)
q_down = Vector{Float64}(undef, size_Cf_R)


for i in 1:size_Cf_R

    cost_Rup_SCGT_i = cost_Rup_SCGT[i]
    cost_Rup_CCGT_i = cost_Rup_CCGT[i]
    cost_Rdown_SCGT_i = cost_Rdown_SCGT[i]
    cost_Rdown_CCGT_i = cost_Rdown_CCGT[i]
    
    Cf_Rup[i] = vcat(fill(cost_Rup_SCGT_i, nSCGT), fill(cost_Rup_CCGT_i, nCCGT))
    Cf_Rdown[i] = vcat(fill(cost_Rdown_SCGT_i, nSCGT), fill(cost_Rdown_CCGT_i, nCCGT))
    
    results_RI_NI = reserve_RI_NI(NI_up_og, NI_down_og, Cf_Rup[i],  Cf_Rdown[i])
    results_RI = reserve_RI(Cf_Rup[i], Cf_Rdown[i])
    rf_up, rf_down = results_RI_NI[end].rf_up, results_RI_NI[end].rf_down
    q_up[i], q_down[i] = get_your_qs(results_RI_NI[end], results_RI[end])
    Lambda_r_up, Lambda_r_down, fig = merit_order_price_reserve(rf_up, rf_down, Cf_Rup[i], Cf_Rdown[i]; dir = "down regulation")
    _, _, fig_up = merit_order_price_reserve(rf_up, rf_down, Cf_Rup[i], Cf_Rdown[i]; dir = "up regulation")
    Lambda_R = (Lambda_r_up, Lambda_r_down)
    R_bids = [rf_up, rf_down]

    res_BM = run_kfold_simulation("Benchmark", Pw, nS_all, prob, 0, 0, k, fold_size, 
                            rf_up, rf_down, NI_up_og, NI_down_og, Cf_Rup[i], Cf_Rdown[i], R_bids, Lambda_R, hour)
                            
    R_GDCA= run_kfold_simulation("TA_bin_GDCA", Pw, nS_all, prob, q_up[i], q_down[i], k, fold_size, 
                            rf_up, rf_down, NI_up_og, NI_down_og, Cf_Rup[i], Cf_Rdown[i], R_bids, Lambda_R, hour; history = res_BM.history.history_insample, max_iter = 1000, penalty_t = penalty, ind = true, bin = true )
    
    iter_bin, iter_bin_LP = extract_iters_from_run(R_GDCA)
    push!(iter_bin_all, iter_bin)
    push!(iter_bin_LP_all, iter_bin_LP)
    push!(all_Cf_R_BM, res_BM)
    push!(all_Cf_R_GDCA, R_GDCA)


end


# after the loop
data = Dict(
    "meta" => Dict(
        "k" => k,
        "fold_size" => fold_size,
        "nS_all" => nS_all,
        "prob" => prob,
        "cost_Rup_SCGT" => collect(cost_Rup_SCGT),
        "cost_Rdown_SCGT" => collect(cost_Rdown_SCGT),
        "cost_Rup_CCGT" => collect(cost_Rup_CCGT),
        "cost_Rdown_CCGT" => collect(cost_Rdown_CCGT),
    ),
    "BM" => all_Cf_R_BM,
    "GDCA" => all_Cf_R_GDCA,
    "q_up" => q_up,
    "q_down" => q_down,
    "iter_bin_all" => iter_bin_all,
    "iter_bin_LP_all" => iter_bin_LP_all,
    "penalty" => penalty
)

jobid = get(ENV, "LSB_JOBID", "manual")
@save joinpath(outdir, "Cf_" * jobid * ".jld2") data







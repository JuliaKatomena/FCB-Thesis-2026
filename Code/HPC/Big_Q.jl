import Pkg
using Pkg
Pkg.activate(@__DIR__)
cd(@__DIR__)


using JuMP, Plots, HiGHS, Statistics, LinearAlgebra, LaTeXStrings, DataFrames, CSV, XLSX, StatsPlots, NPZ, UnPack, JLD2, PyCall, Random

include("constants.jl")
include("opt_models.jl")
include("outputs_processing.jl")

# -------------------------------------------------------------------
# Output directory
# -------------------------------------------------------------------
outdir = joinpath(@__DIR__, "export_HPC")
mkpath(outdir)

# -------------------------------------------------------------------
# Scenario preprocessing
# -------------------------------------------------------------------
file_list = ["seasonal_scenarios_0125.csv"]

Pw, centers, Prob, nS_all = preprocess_scenario_files(file_list;  season = "Winter")

hour = 1
prob = round.(Prob[1][:, hour], digits = 5)

k = 10
fold_size = div(nS_all, k)
Pw_hour = Pw[hour]
n_hours = size(Pw, 1)

# -------------------------------------------------------------------
# Helper function
# -------------------------------------------------------------------
function extract_iters_from_run(run)
    h = run[:history][:history_insample]

    # Julia is 1-indexed
    iter_bin    = [x[1].iter for x in h]
    iter_bin_LP = [x[2].iter for x in h]

    return iter_bin, iter_bin_LP
end

# -------------------------------------------------------------------
# Parameters
# -------------------------------------------------------------------
penalty = 1000

# -------------------------------------------------------------------
# Reserve calculations
# -------------------------------------------------------------------
results_RI_NI = reserve_RI_NI( NI_up_og, NI_down_og, Cf_Rup_og, Cf_Rdown_og)

results_RI = reserve_RI(Cf_Rup_og, Cf_Rdown_og)

rf_up   = results_RI_NI[end].rf_up
rf_down = results_RI_NI[end].rf_down

q_up, q_down = get_your_qs(  results_RI_NI[end], results_RI[end])
Lambda_r_up, Lambda_r_down, fig =
    merit_order_price_reserve( rf_up, rf_down,Cf_Rup_og, Cf_Rdown_og; dir = "down regulation" )

Lambda_R = (Lambda_r_up, Lambda_r_down)
R_bids   = [rf_up, rf_down]

q = q_up + q_down

# -------------------------------------------------------------------
# Run benchmark
# -------------------------------------------------------------------
res_BM = run_kfold_simulation("Benchmark", Pw, nS_all, prob, 0, 0, k, fold_size,
    rf_up, rf_down, NI_up_og, NI_down_og, Cf_Rup_og, Cf_Rdown_og, R_bids, Lambda_R, hour)
# -------------------------------------------------------------------
# Run GDCA
# -------------------------------------------------------------------
R_GDCA = run_kfold_simulation("TA_bin_GDCA", Pw, nS_all, prob, q, q, k, fold_size,
    rf_up, rf_down, NI_up_og, NI_down_og, Cf_Rup_og, Cf_Rdown_og, R_bids, Lambda_R, hour;
    history  = res_BM.history.history_insample, max_iter = 800, penalty_t = penalty)

iter_bin, iter_bin_LP = extract_iters_from_run(R_GDCA)
# -------------------------------------------------------------------
# Save results
# -------------------------------------------------------------------
data = Dict(
    "BM" => res_BM,
    "GDCA" => R_GDCA,
    "q_up" => q_up,
    "q_down" => q_down,
    "iter_bin_all" => iter_bin,
    "iter_bin_LP_all" => iter_bin_LP,
    "penalty" => penalty
)

jobid = get(ENV, "LSB_JOBID", "manual")
@save joinpath(outdir, "Big_Q_" * jobid * ".jld2") data
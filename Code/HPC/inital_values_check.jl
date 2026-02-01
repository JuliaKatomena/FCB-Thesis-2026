using Pkg
Pkg.activate(@__DIR__)

using JuMP, Plots, HiGHS,  Statistics, LinearAlgebra, LaTeXStrings, DataFrames, CSV, XLSX, StatsPlots, NPZ, UnPack, JLD2, Random

cd(@__DIR__)

include("constants.jl")
include("opt_models.jl")
include("output_pros.jl")

outdir = joinpath(@__DIR__, "export_HPC")
mkpath(outdir)  # creates folder if missing (and parents too)

jobid = get(ENV, "LSB_JOBID", "manual")
checkpoint_file = joinpath(outdir, "OAT_checkpoint_" * jobid * ".jld2")
save_every = 5   # save every 5 model runs (adjust)



results_RI_NI = reserve_RI_NI(NI_up_og, NI_down_og, Cf_Rup_og, Cf_Rdown_og)
results_RI = reserve_RI(Cf_Rup_og, Cf_Rdown_og)
rf_up, rf_down = results_RI_NI[end].rf_up, results_RI_NI[end].rf_down
q_up, q_down = get_your_qs(results_RI_NI[end], results_RI[end])
Lambda_r_up, Lambda_r_down, _ = merit_order_price_reserve(rf_up, rf_down, Cf_Rup_og, Cf_Rdown_og; dir = "down regulation")
Lambda_R = (Lambda_r_up, Lambda_r_down)
R_bids = [rf_up, rf_down]

file_list = ["seasonal_scenarios_0125.csv"]

Pw, centers, Prob, nS_all = preprocess_scenario_files(file_list; season="Spring")         # Considering that they all have the same uncertainty

hour = 1
prob = Vector{Float64}(undef, nS_all)
prob = round.(Prob[1][:, hour], digits = 5)     # in this case we prob is only constituted of 1 element

k = 10
fold_size = div(nS_all, k)
Pw_hour = Pw[hour]
n_hours = size(Pw)[1]

# Benchmark 
#fold 1 by setting k = 1
fold = 1
res_BM =  run_kfold_simulation("Benchmark", Pw, nS_all, prob, 0, 0, fold, fold_size, 
                            rf_up, rf_down,NI_up_og, NI_down_og, Cf_Rup_og,  Cf_Rdown_og, R_bids, Lambda_R, hour)

                            # We only look at the first fold
in_sample_Pw, out_sample_Pw, in_sample_prob, out_sample_prob, onS = sample_split(Pw[hour], nS_all, prob, fold_size, fold)

# We extract the intitial values from the benchmark, this will be our basis 
@unpack k0_up, k0_down, ppB0_up, ppB0_down, z0_up, z0_down, Imb_abs_0, s0, s_up0, s_down0 = extract_initial_values_from_benchmark(res_BM.history.history_insample[fold], fold_size)

base_kwargs = Dict{Symbol,Any}(
    :k0_up => k0_up,
    :k0_down => k0_down,
    :ppB0_up => ppB0_up,
    :ppB0_down => ppB0_down,
    :z0_up => z0_up,
    :z0_down => z0_down,
    :s_up0 => s_up0,
    :s_down0 => s_down0,
)

base_kwargs[:q_up] = q_up
base_kwargs[:q_down] = q_down
base_kwargs[:rf_up] = rf_up
base_kwargs[:rf_down] = rf_down
base_kwargs[:NI_up] = NI_up_og
base_kwargs[:NI_down] = NI_down_og

k0_up_range  = collect(0.0:0.5:1.0)
k0_down_range  = collect(0.0:0.5:1.0)
z0_up_range = collect(0.0:NI_up_og/2:NI_up_og)
z0_down_range  = collect(0.0:NI_down_og/2:NI_down_og)
ppB0_up_range = collect(0.0:50.0:100.0)
ppB0_down_range  = collect(0.0:50.0:100.0)
s0_up_range = collect(0.0:1.0:1.0)
s0_down_range  = collect(0.0:1.0:1.0)

function run_model_experiments_OAT(
    Pw, nS, prob;
    model::Symbol,
    sweep_vars::Dict{Symbol,<:AbstractVector},
    base_kwargs = Dict(),
    checkpoint_file::Union{Nothing,String}=nothing,
    Cf_Rup_og, 
    Cf_Rdown_og,
    Lambda_R,
    save_every::Int=10)

    # ------------------------------------------------------------------
    # Two-sided One-At-A-Time (OAT) sensitivity analysis.
    #
    # For each sweep variable:
    #   - run a baseline at lower bounds
    #   - perturb the variable upward
    #   - run a baseline at upper bounds
    #   - perturb the variable downward
    #
    # Scalars are expanded to arrays using base_kwargs shapes.
    # ------------------------------------------------------------------

    results = []
    var_names = collect(keys(sweep_vars))

    # Helper: expand scalar to array if needed
    function expand(vname, vval)
        if haskey(base_kwargs, vname)
            return fill(vval, size(base_kwargs[vname])...)
        else
            return vval
        end
    end

    function run_and_store!(results, Pw, nS, prob, model;
                            kwargs,
                            sweep_vals,
                            baseline_side,
                            experiment,
                            active_var)

        #println(kwargs)
        DA_bids, obj, h_end = run_model(Pw, nS, prob, model; kwargs...)

        profits_all, _, _, _ = participant_profits(
            h_end, nS, prob, Cf_Rup_og, Cf_Rdown_og; 
            Q_up = kwargs[:q_up], 
            Q_down = kwargs[:q_down],
            Lambda_R = Lambda_R 
        )

        result = Dict(
            :ppDA          => DA_bids[1],
            :pfDA          => DA_bids[2],
            :pgDA          => DA_bids[3],
            :objective     => obj,
            :profit_pol    => profits_all[1], # [pol, flex, reg]
            :profit_flex   => profits_all[2],
            :profit_reg    => profits_all[3],
            :baseline_side => baseline_side,
            :experiment    => experiment,
            :active_var    => active_var
        )

        # Reszta kodu (save sweep_vals, checkpointing) pozostaje bez zmian
        for (vname, vval) in sweep_vals
            result[vname] = vval
        end
        push!(results, result)
    
        # --- checkpoint save ---
        if checkpoint_file !== nothing && (length(results) % save_every == 0)
            tmp = checkpoint_file * ".tmp"
            data_partial = Dict(
                "results" => results,
                "meta" => Dict(
                    "hour" => hour,
                    "fold" => fold,
                    "model" => model,
                    "saved_at_len" => length(results),
                )
            )
            @save tmp data_partial
            mv(tmp, checkpoint_file; force=true)
            println("Init_val_checkpoint: ", checkpoint_file, " (", length(results), " runs)")
        end

    end
    
    # --------------------------------------------------------------
    # LOWER-BOUND BASELINE
    # --------------------------------------------------------------
    baseline_lb_kwargs = deepcopy(base_kwargs)
    baseline_lb_vals   = Dict{Symbol,Any}()

    for (vname, vrange) in sweep_vars
        lb = first(vrange)
        baseline_lb_vals[vname] = lb
        baseline_lb_kwargs[vname] = expand(vname, lb)
        #println("vname:", vname)
        #println("vrange:", vrange)
    end

    #println("→ Running baseline at LOWER bounds")
    #println("baseline vals", baseline_lb_vals)
    run_and_store!(
        results, Pw, nS, prob, model;
        kwargs        = baseline_lb_kwargs,
        sweep_vals   = baseline_lb_vals,
        baseline_side = :LB,
        experiment    = :baseline,
        active_var    = :none)

    # OAT upward from lower bound
    for (vname, vrange) in sweep_vars
        lb = first(vrange)

        for vval in vrange
            vval == lb && continue

            #println("→ OAT from LB: $vname = $vval")

            kwargs = deepcopy(baseline_lb_kwargs)
            kwargs[vname] = expand(vname, vval)

            sweep_vals = deepcopy(baseline_lb_vals)
            sweep_vals[vname] = vval

            run_and_store!(
                results, Pw, nS, prob, model;
                kwargs        = kwargs,
                sweep_vals   = sweep_vals,
                baseline_side = :LB,
                experiment    = :OAT,
                active_var    = vname
            )
        end
    end

    #println("Running upper bounds now")

    # --------------------------------------------------------------
    # UPPER-BOUND BASELINE
    # --------------------------------------------------------------
    baseline_ub_kwargs = deepcopy(base_kwargs)
    baseline_ub_vals   = Dict{Symbol,Any}()

    for (vname, vrange) in sweep_vars
        ub = last(vrange)
        baseline_ub_vals[vname] = ub
        baseline_ub_kwargs[vname] = expand(vname, ub)
    end

    #println("→ Running baseline at UPPER bounds")
    run_and_store!(
        results, Pw, nS, prob, model;
        kwargs        = baseline_ub_kwargs,
        sweep_vals   = baseline_ub_vals,
        baseline_side = :UB,
        experiment    = :baseline,
        active_var    = :none)

    # OAT downward from upper bound
    for (vname, vrange) in sweep_vars
        ub = last(vrange)

        for vval in vrange
            vval == ub && continue

            #println("→ OAT from UB: $vname = $vval")

            kwargs = deepcopy(baseline_ub_kwargs)
            kwargs[vname] = expand(vname, vval)

            sweep_vals = deepcopy(baseline_ub_vals)
            sweep_vals[vname] = vval

            run_and_store!(
            results, Pw, nS, prob, model;
            kwargs        = kwargs,
            sweep_vals   = sweep_vals,
            baseline_side = :UB,
            experiment    = :OAT,
            active_var    = vname
            )
        end
    end

    return results
end


# Variables set to lower bounds. 
results_LP = run_model_experiments_OAT(in_sample_Pw, fold_size, in_sample_prob; model =:bin_LP,
    sweep_vars = Dict(:k0_up => k0_up_range, :k0_down => k0_down_range, :z0_up => z0_up_range, :z0_down => z0_down_range,
                         :s_up0 => s0_up_range, :s_down0 => s0_down_range, :ppB0_up => ppB0_up_range, :ppB0_down => ppB0_down_range),
    base_kwargs=base_kwargs,
    checkpoint_file = checkpoint_file,
    Cf_Rup_og = Cf_Rup_og, 
    Cf_Rdown_og = Cf_Rdown_og,
    Lambda_R = Lambda_R,
    save_every = save_every
    )

objectives = getindex.(results_LP, :objective)
std(objectives)


ppDAs = getindex.(results_LP, :ppDA)
ppDAs_std_devs = round.(mapslices(std, ppDAs; dims=1), digits=2)  # std along each row
ppDAs_means = round.(mapslices(mean, ppDAs; dims=1), digits=2)  # mean along each row

pfDAs = getindex.(results_LP, :pfDA)
pfDAs_std_devs = round.(mapslices(std, pfDAs; dims=1), digits=2)  # std along each row
pfDAs_means = round.(mapslices(mean, pfDAs; dims=1), digits=2)  # mean along each row


pgDAs = getindex.(results_LP, :pgDA)
pgDAs_std_devs = round.(mapslices(std, pgDAs; dims=1), digits=2)  # std along each row
pgDAs_means = round.(mapslices(mean, pgDAs; dims=1), digits=2)  # mean along each row

profits_pol = getindex.(results_LP, :profit_pol)
profits_pol_means = round.(mapslices(mean, profits_pol; dims=1), digits=2)
profits_pol_std = round.(mapslices(std, profits_pol; dims=1), digits=2)

profits_flex = getindex.(results_LP, :profit_flex)
profits_flex_means = round.(mapslices(mean, profits_flex; dims=1), digits=2)
profits_flex_std = round.(mapslices(std, profits_flex; dims=1), digits=2)

profits_reg = getindex.(results_LP, :profit_reg)
profits_reg_means = round.(mapslices(mean, profits_reg; dims=1), digits=2)
profits_reg_std = round.(mapslices(std, profits_reg; dims=1), digits=2)


data = Dict(
    "results" => results_LP,
    "meta" => Dict(
        "hour" => hour,
        "fold" => fold,
        "fold_size" => fold_size,
        "model" => :bin,
    ),
    "mean_ppDA" => ppDAs_means,
    "std_ppDA" => ppDAs_std_devs,
    "mean_pfDA" => pfDAs_means,
    "std_pfDA" => pfDAs_std_devs,
    "mean_pgDA" => pgDAs_means,
    "std_pgDA" => pgDAs_std_devs,
    "mean_profit_pol" => profits_pol_means,
    "std_profit_pol" => profits_pol_std,
    "mean_profit_flex" => profits_flex_means,
    "std_profit_flex" => profits_flex_std,
    "mean_profit_reg" => profits_reg_means,
    "std_profit_reg" => profits_reg_std,

)

final_file = joinpath(outdir, "Init_val_final_" * jobid * ".jld2")
@save final_file data


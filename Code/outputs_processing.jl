####################################################################################

                        # IND. OPT. OUTPUT COLLECTION

function drop_polluter(h::NamedTuple, pol::Int)
    per_polluter = (:ppDA,)
    per_matrix   = (:ppB, :k_up, :ppB_up, :diff_up,
                    :k_down, :ppB_down, :diff_down,
                    :δ_up, :δ_down)

    return NamedTuple{keys(h)}((
        if k in per_matrix
            v[:, setdiff(1:size(v,2), [pol])]
        elseif k in per_polluter
            v[setdiff(1:length(v), [pol])]
        else
            v
        end
        for (k,v) in Base.pairs(h)   # <-- force Base.pairs
    ))
end

####################################################################################

                        # GENERAL OUTPUT COLLECTION

################# Fixing Binaries for LP
function get_fixed_bin(df)
    # make sure it is sorted consistently
    sort!(df, [:Scenario, :Gen])

    # y: reshape into nS × nP
    y_vec = Float64.(df[!, :y])
    Y = reshape(y_vec, nP, nS)'   # size (nS, nP)

    # s_up, s_down: one value per scenario (take first row of each scenario block)
    s_up_vec   = Float64.(df[1:nP:end, :s_up])
    s_down_vec = Float64.(df[1:nP:end, :s_down])

    @assert size(Y) == (nS, nP)
    @assert length(s_up_vec) == nS
    @assert length(s_down_vec) == nS
    return Y, s_up_vec, s_down_vec 
end

function _extract_namedtuple(history_all, vars::Vector{Symbol})
    # Case 1: already a NamedTuple
    if history_all isa NamedTuple
        return history_all
    end

    # Case 2: Vector of stuff → recurse on the last element
    if history_all isa AbstractVector
        return _extract_namedtuple(history_all[end], vars)
    end

    # Case 3: Tuple of things (e.g. (nt1, nt2, ...))
    if history_all isa Tuple
        # Prefer a NamedTuple that has at least one requested var
        for x in history_all
            if x isa NamedTuple
                fns = fieldnames(typeof(x))
                if any(v -> v in fns, vars)
                    return x
                end
            end
        end
        # Fallback: first NamedTuple if nothing matches vars
        for x in history_all
            x isa NamedTuple && return x
        end
        error("No NamedTuple found inside tuple history_all of type $(typeof(history_all))")
    end

    error("Unsupported history_all type: $(typeof(history_all))")
end


# Getting Data Frame from the history
function last_iter_results_new(history_all, nS;
    pol::Union{Nothing, Int}=nothing,
    flex::Union{Nothing, Int}=nothing,
    reg::Union{Nothing, Int}=nothing,
    vars = Symbol[], digits = 4)

    # --- Robustly pick the NamedTuple that actually has the variables ---
    h = _extract_namedtuple(history_all, vars)
    fnames = fieldnames(typeof(h))
    roundmaybe(x) = x isa AbstractFloat ? round(x; digits=digits) : x


    # === Generator types ===
    gen_types = [
        ("Polluter", "pp", nP, pol),
        ("Flexible", "pf", nF, flex),
        ("Regular",  "pg", nR, reg)
    ]

    # === Define generator-level fields ===
    POLLUTER_FIELDS = Set(Symbol.([
        :ppDA, :ppDA_cost, :ppDA_single_cost,
        :ppB, :ppB_up, :ppB_down,
        :ppB_cost, :ppB_single_cost,
        :k_up, :k_down, :δ_up, :δ_down,
        :diff_up, :diff_down, :t, 
        :diff_up_pct, :diff_down_pct
    ]))

    FLEXIBLE_FIELDS = Set(Symbol.([
        :pfDA, :pfDA_cost, :pfDA_single_cost, :pfB,
        :pfB_up, :pfB_down, :pfB_cost, :pfB_single_cost
    ]))

    REGULAR_FIELDS = Set(Symbol.([:pg1DA]))

    COMMON_SCENARIO_FIELDS = Set(Symbol.([
        :s, :s_up, :s_down, :Imb, :Imb_abs, :Lambda_DA, :Lambda_B, :obj
    ]))
    POLLUTER_ONLY_FIELDS = Set(Symbol.([
        :z_up, :z_down
    ]))

    # === Output containers ===
    df_polluters = DataFrame()
    df_flexibles = DataFrame()
    df_regulars  = DataFrame()

    # === Helper to pick correct field set ===
    function get_fieldset(prefix)
        prefix == "pp" && return POLLUTER_FIELDS
        prefix == "pf" && return FLEXIBLE_FIELDS
        prefix == "pg" && return REGULAR_FIELDS
        return Set{Symbol}()
    end

    # === Main loop ===
    for (gen_type, prefix, nG, selected_idx) in gen_types
        if nG == 0
            continue
        end

        # If a specific generator index is provided, only process that one
        if isnothing(selected_idx)
            gens = 1:nG
        else
            gens = [selected_idx]
        end


        GEN_FIELDS = get_fieldset(prefix)
        rows = NamedTuple[]

        for g in gens, w in 1:nS
            iter_val = :iter in fnames ? getfield(h, :iter) : 1
            row = (iter = iter_val, Scenario = w, Gen = g)

            # === 1. Generator-specific variables ===
            for name in vars
                if (name in [:z_up, :z_down]) && prefix != "pp"
                    continue
                end
                if !(name in fnames)
                    continue
                end

                val = getfield(h, name)
                entry = nothing

                if val isa Number
                    entry = roundmaybe(val)
                elseif val isa AbstractVector
                    if name in GEN_FIELDS && nG >= 1 && g <= length(val)
                        entry = roundmaybe(val[g])
                    elseif length(val) == nS && w <= nS
                        entry = roundmaybe(val[w])
                    elseif length(val) == 1
                        entry = roundmaybe(val[1])
                    end
                elseif val isa AbstractMatrix
                    if size(val) == (nS, nG) && w <= nS && g <= nG
                        entry = roundmaybe(val[w, g])
                    elseif size(val) == (nS, 1) && nG == 1
                        entry = roundmaybe(val[w, 1])
                    end
                end

                if entry !== nothing
                    row = merge(row, NamedTuple{(name,)}((entry,)))
                end
            end

            # === 2. Scenario-level fields ===
            scenario_fields = prefix == "pp" ?
                union(COMMON_SCENARIO_FIELDS, POLLUTER_ONLY_FIELDS) :
                COMMON_SCENARIO_FIELDS

            for name in intersect(scenario_fields, Set(vars))
                if !(name in fnames)
                    continue
                end

                val = getfield(h, name)
                entry = nothing

                if val isa Number
                    entry = roundmaybe(val)
                elseif val isa AbstractVector && length(val) >= w
                    entry = roundmaybe(val[w])
                elseif val isa AbstractMatrix && size(val, 1) >= w
                    entry = roundmaybe(val[w, 1])
                end

                if entry !== nothing
                    row = merge(row, NamedTuple{(name,)}((entry,)))
                end
            end

            push!(rows, row)
        end

        df_ref = DataFrame(Tables.columntable(rows))

        if prefix == "pp"
            df_polluters = df_ref
        elseif prefix == "pf"
            df_flexibles = df_ref
        elseif prefix == "pg"
            df_regulars = df_ref
        end
    end

    return df_polluters, df_flexibles, df_regulars
end

#Plotting Convergence
function plot_history_vars(history; vars::Vector{Symbol}, scenario=1, generator=1)
    n_iter = history[end].iter

    scen_range = isa(scenario, Integer) ? [scenario] : collect(scenario)
    gen_range  = isa(generator, Integer) ? [generator] : collect(generator)

    plt = plot(title="Convergence across scenarios and generators",
               xlabel="Iteration", ylabel="Value")

    for var in vars
        for s in scen_range
            for g in gen_range
                values = Float64[]
                for i in 1:n_iter
                    val = getproperty(history[i], var)

                    if ndims(val) == 2
                        push!(values, val[s, g])
                    elseif ndims(val) == 1 && length(val) >= s
                        push!(values, val[s])
                    else
                        push!(values, val)
                    end
                end

                plot!(plt, values, label="$(var), s=$s, g=$g")
            end
        end
    end
    return plt
end

################# Fixing Binaries for LP
function extract_binaries(history, nS, nP; digits=4)
    h = history[end]

    # y: reshape into nS × nP
    if :y in fieldnames(typeof(h))
        y_val = getfield(h, :y)
        if y_val isa AbstractVector
            Y = reshape(Float64.(y_val), nP, nS)'
        elseif y_val isa AbstractMatrix
            Y = Float64.(y_val)
        else
            error("Unsupported format for y: $(typeof(y_val))")
        end
    else
        Y = nothing
    end

    # s_up, s_down: scenario-level values
    s_up_vec   = (:s_up in fieldnames(typeof(h)))   ? Float64.(getfield(h, :s_up))   : nothing
    s_down_vec = (:s_down in fieldnames(typeof(h))) ? Float64.(getfield(h, :s_down)) : nothing

    return Y, s_up_vec, s_down_vec
end


function sample_split( Pw_all,nS_all, prob_all, fold_size, k)
  

    in_indices = ((k-1)*fold_size + 1):(k*fold_size)
    out_indices = setdiff(1:nS_all, in_indices)

    out_sample_Pw = Pw_all[out_indices, :]
    in_sample_Pw = Pw_all[in_indices, :]

    out_sample_prob = prob_all[out_indices]/sum(prob_all[out_indices])
    in_sample_prob = prob_all[in_indices]/sum(prob_all[in_indices])

    out_size = length(out_indices)

    return in_sample_Pw, out_sample_Pw, in_sample_prob, out_sample_prob, out_size

end 
function extract_initial_values_from_benchmark(history, inS)
    ppB0      = history.ppB
    ppB0_up   = history.ppB_up
    ppB0_down = history.ppB_down
    z0_up     = sum(ppB0_up, dims = 2)
    z0_down   = sum(ppB0_down, dims = 2)

    # Avoid division by zero
    k0_up   = similar(ppB0_up)
    k0_down = similar(ppB0_down)
    for w in 1:inS, p in 1:nP
        k0_up[w, p]   = z0_up[w] != 0   ? ppB0_up[w, p]   / z0_up[w]   : 0.0
        k0_down[w, p] = z0_down[w] != 0 ? ppB0_down[w, p] / z0_down[w] : 0.0
    end

    imb = sum(ppB0, dims = 2)
    Imb_abs_0  = abs.(imb)

    s_up0 = zeros(inS)
    s_down0 = zeros(inS)
    s0 = zeros(inS)

    for w in 1:inS
        if imb[w] > 0 
            s_up0[w] = 0
            s_down0[w] = 1
            s0[w] = -1
        elseif imb[w] < 0
            s_up0[w] = 1   
            s_down0[w] = 0
            s0[w] = 1
        end
    end

    return Dict(
        :ppB0 => ppB0,
        :ppB0_up => ppB0_up,
        :ppB0_down => ppB0_down,
        :k0_up => k0_up,
        :k0_down => k0_down,
        :z0_up => z0_up,
        :z0_down => z0_down,
        :Imb_abs_0 => Imb_abs_0,
        :s0 => s0,
        :s_up0 => s_up0,
        :s_down0 => s_down0
    )
end


####################################################################################

                        # PROFIT COMPUTATION

function participant_profits(history, Cf_Rup, Cf_Rdown; Q_up = nothing, Q_down = nothing, Lambda_DA = nothing, Lambda_B = nothing, Lambda_R = nothing)
    profits_pol = fill(0.0, nP)
    profits_flex = fill(0.0, nF)
    profits_reg = fill(0.0, nR)

    pol_profit = Vector{Vector{Float64}}(undef, nP)  # array of vectors, initially undefined
    flex_profit = Vector{Vector{Float64}}(undef, nF)  # array of vectors, initially undefined

    # Use provided lambdas if available, otherwise from history
    Lambda_DA = isnothing(Lambda_DA) ? history.Lambda_DA : Lambda_DA
    Lambda_B  = isnothing(Lambda_B)  ? history.Lambda_B  : Lambda_B

    # Check if penalty should be included
    include_penalty = !(isnothing(Q_up) || isnothing(Q_down))

    for p in 1:nP
        profit_DA = history.ppDA[p] * (Lambda_DA - Cp[p])
        profit_B = sum(prob[w] * (Lambda_B[w] - Cp[p]) * history.ppB[w, p] for w in 1:nS)
        penalty = include_penalty ? sum(prob[w] * (Q_down * history.k_down[w, p] + Q_up * history.k_up[w, p]) for w in 1:nS) : 0.0
        profits_pol[p] = profit_DA + profit_B - penalty
        pol_profit[p] = [profit_DA, profit_B, - penalty]
    end

    for f in 1:nF
        profit_DA = history.pfDA[f] * (Lambda_DA - Cf[f])
        profit_B = sum(prob[w] * (Lambda_B[w] - Cf[f]) * (history.pfB_up[w, f] - history.pfB_down[w, f]) for w in 1:nS)
        profit_r = sum(history.rf_up[f] * (Lambda_R[1] - Cf_Rup[f]) + history.rf_down[f] * (Lambda_R[2] - Cf_Rdown[f]))
        profits_flex[f] = profit_DA + profit_B + profit_r
        flex_profit[f] = [profit_DA, profit_B, profit_r]
    end

    for r in 1:nR
        profit_DA = history.pgDA[r] * (Lambda_DA - Cr[r])
        profits_reg[r] = profit_DA
    end

    return [profits_pol, profits_flex, profits_reg], pol_profit, flex_profit, profits_reg 
end

# This is the one we are using at last
function profits(nS, prob, DA_bids, B_bids, R_bids, Lambda_DA, Lambda_B, Lambda_R, Cf_Rup, Cf_Rdown; Q_up = nothing, Q_down = nothing, tol = 1e-6)
    #  B_bids_cen_in = [B_bids_pol_cen_in, B_bids_pol_up_cen_in, B_bids_pol_down_cen_in, act_bids_up_cen_in, act_bids_down_cen_in]
    
    pol_profit = Vector{Vector{Float64}}(undef, nP)  # array of vectors, initially undefined
    flex_profit = Vector{Vector{Float64}}(undef, nF)  # array of vectors, initially undefined
    reg_profit = Vector{Vector{Float64}}(undef, nR)  # array of vectors, initially undefined

    # Check if penalty should be included
    include_penalty = !(isnothing(Q_up) || isnothing(Q_down))

    k_up = fill(0.0, nS, nP)
    k_down = fill(0.0, nS, nP)

    sum_plus = 0.0
    sum_minus = 0.0
    sum_none = 0.0
    for w in 1:nS
        total_imbalance = sum(B_bids[1][w, :]) # This is summing over all polluters
        total_positive = sum((B_bids[1][w, p] for p in 1:nP if B_bids[1][w, p] >= 0); init = 0.0)
        total_negative = sum((B_bids[1][w, p] for p in 1:nP if B_bids[1][w, p] < 0); init = 0.0)
        if total_imbalance > tol
            sum_plus += 1
        elseif total_imbalance < - tol
            sum_minus += 1
        else
            sum_none += 1
        end


        for p in 1:nP
            bid =  B_bids[1][w, p]
            
            if total_imbalance > tol && bid > 0 
                k_down[w, p] = bid / total_positive
            elseif total_imbalance < - tol && bid < 0 
                k_up[w, p] = bid  / total_negative
            end
        end
        ##display("kup", k_up[w, :])
        ##display("kdown", k_down[w, :])
    end

    #println("Number of positive imbalances: ", sum_plus)
    #println("Number of negative imbalances: ", sum_minus)
    #println("Number of zero imbalances: ", sum_none)


    for p in 1:nP
        profit_DA = DA_bids[1][p] * (Lambda_DA - Cp[p])
        profit_B = sum(prob[w] * (Lambda_B[w] - Cp[p]) * B_bids[1][w, p] for w in 1:nS)
        penalty_up = include_penalty ? sum(prob[w] * ((Q_up) * k_up[w, p]) for w in 1:nS) : 0.0
        penalty_down = include_penalty ? sum(prob[w] * ((Q_down) * k_down[w, p]) for w in 1:nS) : 0.0
        #profits_pol[p] = profit_DA + profit_B - penalty
        #println("Polluter $p: penalty_up = $penalty_up, penalty_down = $penalty_down")
        pol_profit[p] = [profit_DA, profit_B, - penalty_up, - penalty_down]
    end
    #println("Summed penalty down: ", sum(pol_profit[p][4] for p in 1:nP))
    #println("Summed penalty up: ", sum(pol_profit[p][3] for p in 1:nP))

    for f in 1:nF
        profit_DA = DA_bids[2][f] * (Lambda_DA - Cf[f])
        profit_B = sum(prob[w] * (Lambda_B[w] - Cf[f]) * (B_bids[4][w, f] - B_bids[5][w, f]) for w in 1:nS)
        profit_R = sum(R_bids[1][f] * (Lambda_R[1] - Cf_Rup[f]) + R_bids[2][f] * (Lambda_R[2] - Cf_Rdown[f]))
        #profits_flex[f] = profit_DA + profit_B + profit_r
        flex_profit[f] = [profit_DA, profit_B, profit_R]
    end

    for r in 1:nR
        profit_DA = DA_bids[3][r] * (Lambda_DA - Cr[r])
        reg_profit[r] = [profit_DA]
    end

    return pol_profit, flex_profit, reg_profit, [sum_minus, sum_plus, sum_none]
end

function cen_vs_ind_comparison(p_cen, p_ind, DA_cen, DA_ind, B_cen, B_ind, h;digits=2)
    # First element corresponds to 4 polluters, expand them
    # Helper for rounding any array or scalar
    round_all(x) = x isa AbstractArray ? round.(x, digits=digits) : round(x, digits=digits)

    polluter_rows = DataFrame(
        Participant = ["Polluter $i" for i in 1:nP],
        Cost = Cp,
        Profit_Cen = round_all(p_cen[1]),
        Profit_Ind = round_all(p_ind[1]),
        DA_bids_Cen = round_all(DA_cen[1]),
        DA_bids_Ind = round_all(DA_ind[1]),
        B_bids_Cen = round_all(B_cen[1]),
        B_bids_Ind = round_all(B_ind[1]),
    )
    
    flexible_rows = DataFrame(
        Participant = ["Flexible $i" for i in 1:nF],
        Cost = Cf,
        Profit_Cen = round_all(p_cen[2]),
        Profit_Ind = round_all(p_ind[2]),
        DA_bids_Cen = round_all(DA_cen[2]),
        DA_bids_Ind = round_all(DA_ind[2]),
        B_bids_Cen = round_all(B_cen[2]),
        B_bids_Ind = round_all(B_ind[2]),
    )

    regular_rows = DataFrame(
        Participant = ["Regular $i" for i in 1:nR],
        Cost = Cr,
        Profit_Cen = round_all(p_cen[end]),
        Profit_Ind = round_all(p_ind[end]),
        DA_bids_Cen = round_all(DA_cen[end]),
        DA_bids_Ind = round_all(DA_ind[end]),
        B_bids_Cen = 0.0,
        B_bids_Ind = 0.0,
    )

    # --- Summary rows ---
    total_polluters = DataFrame(
        Participant = ["Total Polluters"],
        Cost = ["-"],
        Profit_Cen = sum(round_all(p_cen[1])),
        Profit_Ind = sum(round_all(p_ind[1])),
        DA_bids_Cen = sum(round_all(DA_cen[1])),
        DA_bids_Ind = sum(round_all(DA_ind[1])),
        B_bids_Cen = sum(round_all(B_cen[1])),
        B_bids_Ind = sum(round_all(B_ind[1])),
    )

    total_flexibles = DataFrame(
        Participant = ["Total Flexibles"],
        Cost = ["-"],
        Profit_Cen = sum(round_all(p_cen[2])),
        Profit_Ind = sum(round_all(p_ind[2])),
        DA_bids_Cen = sum(round_all(DA_cen[2])),
        DA_bids_Ind = sum(round_all(DA_ind[2])),
        B_bids_Cen = sum(round_all(B_cen[2])),
        B_bids_Ind = sum(round_all(B_ind[2])),
    )
    
    total_regular = DataFrame(
        Participant = ["Total Regular"],
        Cost = ["-"],
        Profit_Cen = sum(round_all(p_cen[end])),
        Profit_Ind = sum(round_all(p_ind[end])),
        DA_bids_Cen = sum(round_all(DA_cen[end])),
        DA_bids_Ind = sum(round_all(DA_ind[end])),
        B_bids_Cen = 0.0,
        B_bids_Ind = 0.0,
    )
    # Price information row
    price_row = DataFrame(
        Participant = ["Prices"],
        Cost = "—",
        Profit_Cen = ["DA Price: " * join(round(h.Lambda_DA, digits=digits), ", ")], 
        Profit_Ind = ["Bal Prices: " * join(round.(h.Lambda_B, digits=digits), ", ")], 
        DA_bids_Cen = "—",
        DA_bids_Ind = "—",
        B_bids_Cen = "—",
        B_bids_Ind = "—",
    )

    return vcat(price_row, polluter_rows, total_polluters, flexible_rows, total_flexibles, regular_rows, total_regular)
end

function run_individual_optimizations(h, Cf_Rup, Cf_Rdown; Lambda_R = Lambda_R, 
     polluter_func, flex_func,  polluter_args = NamedTuple())
    profits_ind = [Float64[], Float64[], Float64[]]
    DAbids_ind  = [Float64[], Float64[], Float64[]]
    Bbids_ind   = [Float64[], Float64[]]


    # --- Polluters ---
    for p in 1:nP
        h_new = drop_polluter(h, p)
        history = polluter_func(p, h_new, polluter_args...)
        append!(profits_ind[1], history[end].obj)
        append!(DAbids_ind[1], history[end].ppDA)
        append!(Bbids_ind[1], sum(history[end].ppB[:, 1].*prob))
    end

    # --- Flexible units ---
    for f in 1:nF
        history = flex_func(f, h, h.rf_up, h.rf_down)
        profit_r = sum(h.rf_up[f] * (Lambda_R[1] - Cf_Rup[f]) + h.rf_down[f] * (Lambda_R[2] - Cf_Rdown[f]))
        println("profit from reserve market", profit_r)
        append!(profits_ind[2], history[end].obj + profit_r)
        append!(DAbids_ind[2], history[end].pfDA)
        append!(Bbids_ind[2], sum((history[end].pfB_up[:] - history[end].pfB_down[:]).*prob))
    end

    # --- Regular generator ---
    for r in 1:nR
        history = reg_opt(r, h)
        append!(profits_ind[3], history[end].obj)
        append!(DAbids_ind[3], history[end].pgDA)
    end
    return profits_ind, DAbids_ind, Bbids_ind
end

function plot_balancing_bids(B_cen, B_ind, title; digits=2)
    # Combine into single vector and labels
    participants = vcat(["P. $i" for i in 1:nP], ["F. $i" for i in 1:nF])
  
    B_cen_all = vcat(B_cen[1], B_cen[2])
    B_ind_all = vcat(B_ind[1], B_ind[2])

    # ---- Add total values ----
    total_labels = ["Total Pol.", "Total Flex."]
    total_cen = [sum(B_cen[1]), sum(B_cen[2])]
    total_ind = [sum(B_ind[1]), sum(B_ind[2])]

    participants = vcat(participants, total_labels)
    B_cen_all = vcat(B_cen_all, total_cen)
    B_ind_all = vcat(B_ind_all, total_ind)

    df = DataFrame(
        Participants = participants,
        Centralized = B_cen_all,
        Individual = B_ind_all
    )

    # ---- Numeric x positions ----
    x = 1:nrow(df)
    bar_width = 0.35

    # ---- Plot ----
    StatsPlots.bar(x .- bar_width/2, df.Centralized;
        bar_width,
        label = "Centralized Optimization",
        xlabel = "Participant",
        ylabel = "Bid [MWh]",
        title = title,
        xticks = (x, df.Participants),
        ylims = (minimum(vcat(df.Centralized, df.Individual)) * 1.25,
                maximum(vcat(df.Centralized, df.Individual)) * 1.2),
        color = :steelblue,
        legend = :topleft,
        framestyle = :box,
        rotation = 30,
        bottom_margin = 10Plots.mm) 

    StatsPlots.bar!(x .+ bar_width/2, df.Individual;
        bar_width,
        label = "Individual Optimization",
        color = :orange)

end

# This is also the one being used
function results_cen(DA_bids, R_bids, Pw, nS, prob; Q_up = nothing, Q_down = nothing, Lambda_R = nothing, nb_w = 1)
    Lambda_DA, _ = merit_order_price_DA(DA_bids)

    # We need to run the MOC to get at least the balancing bids for flex
    Lambda_B, B_bids_polluters_cen, act_bids_up, act_bids_down = merit_order_price_B(DA_bids, Pw, nS, R_bids[1], R_bids[2]; nb_w = nb_w) #act_ is for activated_bids
    B_bids_cen = [B_bids_polluters_cen, act_bids_up,  act_bids_down] # Note that there is two balancing bids for flex, 1 for up and 1 for down
    # Conditionally call the version that uses Q_up and Q_down
    if Q_up !== nothing && Q_down !== nothing
        pol_profit, flex_profit, reg_profit = profits(nS, prob, DA_bids, B_bids_cen, R_bids, Lambda_DA, Lambda_B, Lambda_R; Q_up = Q_up, Q_down= Q_down)
    else
        pol_profit, flex_profit, reg_profit = profits(nS, prob, DA_bids, B_bids_cen, R_bids, Lambda_DA, Lambda_B, Lambda_R)
    end

    return pol_profit, flex_profit, reg_profit
end

####################################################################################

                        # COMPUTING Q_UP AND Q_down

function get_your_qs(results_RI_NI, results_RI)
    q_up = results_RI_NI.Lambda_RI_NI_up*results_RI_NI.RI_NI_up - results_RI.Lambda_RI_up*results_RI.RI_up
    q_down = results_RI_NI.Lambda_RI_NI_down*results_RI_NI.RI_NI_down - results_RI.Lambda_RI_down*results_RI.RI_down
    return q_up, q_down
end                          

####################################################################################

                        # ITERATING OVER THE VARIABLES

function run_model_experiments(Pw, nS, prob, models::Vector{Symbol},
                               sweep_vars::Dict{Symbol, <:Any};
                               base_kwargs=Dict())
    results = DataFrame()
    iterations = DataFrame()

    # Create all combinations of sweep variables
    combos = collect(Iterators.product(values(sweep_vars)...))
    var_names = collect(keys(sweep_vars))

    for (combo_idx, combo_vals) in enumerate(combos)
        # Prepare local kwargs for this run
        local_kwargs = deepcopy(base_kwargs)

        # Assign each sweep variable
        for (vname, vval) in zip(var_names, combo_vals)
            if haskey(base_kwargs, vname)
                base_shape = size(base_kwargs[vname])
                # Auto-expand scalars to array shape
                local_kwargs[vname] = isa(vval, Number) ? fill(vval, base_shape...) : vval
            else
                local_kwargs[vname] = vval
            end
        end

        # Save the current sweep values for labeling results later
        sweep_values = Dict(var_names .=> combo_vals)

        for model in models
            println("→ Running $model | " *
                    join(["$v=$(sweep_values[v])" for v in var_names], ", "))

            df_res, df_iter, DA_bids_cen_bin, DA_bids_cen_bin_LP = run_single_model(Pw, nS, prob, model; local_kwargs...)

            # Dynamically add sweep variable columns to results
            for (vname, vval) in sweep_values
                df_res[!, vname] = fill(vval, nrow(df_res))
                df_iter[!, vname] = fill(vval, nrow(df_iter))
            end
            df_res[!, :total] = df_res.polluter .+ df_res.flexible .+ df_res.regular

            append!(results, df_res)
            append!(iterations, df_iter)
        end
    end
    
    plot = profit_vs_initial_values(results, sweep_vars)
    return results, iterations, plot
end

function run_single_model(Pw, nS, prob, model::Symbol; kwargs...)
    if model == :bin
        return run_bin_model(Pw, nS, prob; kwargs...)
    elseif model == :TA_test
        return run_TA_test_model(; kwargs...)
    else
        error("Unknown model: $model")
    end
end

function run_bin_model(Pw, nS, prob; q_up, q_down, rf_up, rf_down, NI_up, NI_down,
                       k0_up, k0_down, z0_up, z0_down,
                       ppB0_up, ppB0_down, s_up0, s_down0,
                       k0_up_singular, k0_down_singular,
                       ppB0_up_singular, ppB0_down_singular)

        # --- run models ---
        history_TA_bin = TA_bin_GDCA(Pw, nS, prob, q_up, q_down, rf_up, rf_down,NI_up, NI_down,  k0_up, k0_down, z0_up, z0_down, ppB0_up, ppB0_down, s_up0, s_down0, penalty_t = 1000, max_iter = 1000, tol = 0.001, alpha = 1.0)
        h = history_TA_bin[end]
        DA_bids_cen_bin = [h.ppDA, h.pfDA, h.pgDA]

        history_bin_LP = TA_bin_LP_GDCA(h, Pw, nS, prob, q_up, q_down, rf_up, rf_down,NI_up, NI_down,  k0_up, k0_down, z0_up, z0_down, ppB0_up, ppB0_down, max_iter = 5000, penalty_t = 1000, tol = 0.01, alpha=1.0)
        h = history_bin_LP[end]
        DA_bids_cen_bin_LP = [h.ppDA, h.pfDA, h.pgDA]

        R_bids = [rf_up, rf_down]
        pol_profit_bin, flex_profit_bin, reg_profit_bin = results_cen(DA_bids_cen_bin, R_bids, out_sample_Pw, onS, out_sample_prob; Lambda_R = Lambda_R, Q_up = q_up, Q_down = q_down)
        pol_profit_bin_LP, flex_profit_bin_LP, reg_profit_bin_LP = results_cen(DA_bids_cen_bin_LP, R_bids, out_sample_Pw, onS, out_sample_prob; Lambda_R = Lambda_R, Q_up = q_up, Q_down = q_down)

        # Individual Opt
    #        profits_bin_LP_ind, _, _ = run_individual_optimizations(h, rf_up, rf_down; polluter_func = polluter_bin_LP, flex_func = flex_bin_LP, polluter_args = (q_up, q_down, k0_up_singular, k0_down_singular, z0_up, z0_down, ppB0_up_singular, ppB0_down_singular))

        # --- summarize profits ---
        df = DataFrame(
            model = ["Cen. Bin", "Cen. Bin-fix LP"],
            polluter = [sum(sum.(pol_profit_bin)), sum(sum.(pol_profit_bin_LP))],
            flexible = [sum(sum.(flex_profit_bin)), sum(sum.(flex_profit_bin_LP))],
            regular  = [sum(sum.(reg_profit_bin)), sum(sum.(reg_profit_bin_LP))]
        )

        df2 = DataFrame(
            model = ["bin", "bin_LP"],
            iterations = [history_TA_bin[end].iter, history_bin_LP[end].iter]
        )

    return df, df2, DA_bids_cen_bin, DA_bids_cen_bin_LP
end

function profit_vs_initial_values(df, sweep_vars::Dict{Symbol, <:Any})
    var_syms = collect(keys(sweep_vars))   # => [:s_up, :s_down]

    x_labels = [join([string(v)*"="*string(row[v]) for v in var_syms], ", ")
            for row in eachrow(df)]

    # --- Numeric x positions ---
    n_combos = length(unique(x_labels))       # number of initial variable combos
    x = 1:n_combos
    bar_width = 0.25                            # width for each bar
    models = unique(df.model)                    # ["Cen. Bin", "Cen. Bin-fix LP", "Ind. Bin-fix LP"]

    # --- Prepare data for each model ---
    # For each model, collect total profits in the order of unique combos
    profits_by_model = Dict(m => [df.total[(x_labels .== combo) .& (df.model .== m)][1] 
                                  for combo in unique(x_labels)] 
                            for m in models)

    # --- Plot ---
    p = bar(
        x .- bar_width, profits_by_model[models[1]];
        bar_width,
        label = models[1],
        color = :steelblue,
        legend = :topright,
        framestyle = :box,
        bottom_margin = 20Plots.mm,
        xrotation = 45,
        left_margin = 10Plots.mm,     # add side margin
        right_margin = 10Plots.mm,    # optional
        top_margin = 10Plots.mm,      # optional
        size = (900, 600)             # increase overall figure size
    )

    bar!(x, profits_by_model[models[2]];
        bar_width,
        label = models[2],
        color = :orange
    )

    # --- Customize axes and ticks ---
    xticks!(x, unique(x_labels))
    ylims!(p, 0, maximum(profits_by_model[models[1]]) * 1.1)
    xlabel!("Initial variable combination")
    ylabel!("Total Profit [€]")
    title!("Total Profit per Model and Initial Values")
    return p
end


function run_kfold_simulation( model::String, Pw, nS_all::Any, prob::Any,
    q_up_in::Any,  q_down_in::Any, k::Any, fold_size::Int,
    rf_up::Any, rf_down::Any, NI_up::Any, NI_down::Any, Cf_Rup::Any, Cf_Rdown::Any, R_bids, Lambda_R, hour;
     history = nothing, max_iter = 100, penalty_up = 100, penalty_down = 100,
      penalty_t = 100, ind::Bool = false, bin::Bool = false)
    # --- 1. Initialize result vectors ---
    profits_outsample = []
    profits_insample = []
    B_bids_insample = []
    B_bids_outsample = []
    DA_bids_insample = []
    DA_bids_outsample = []
    history_insample = [] # Correctly initialized
    lambdas = []
    system_imbalance = []
    profits_outsample_bin = []
    profits_outsample_ind = []
    splits_info = []
    B_bids_all_insample = []
    B_bids_all_outsample = []
 
    # --- 2. K-Fold Loop ---
    for fold in 1:k
        println("Fold: ", fold)
   
        # 2.1. Determine indices and initialize values
        inS = ((fold - 1) * fold_size + 1):(fold * fold_size)
        in_sample_Pw, out_sample_Pw, in_sample_prob, out_sample_prob, onS = sample_split(Pw[hour], nS_all, prob, fold_size, fold)
 
        # 2.2. Model Training (Conditional on 'model' type)
        local h_end, h_end_LP, DA_bids_cen
       
        if model == "TA_bin_Julia"
            @unpack k0_up, k0_down, ppB0_up, ppB0_down, z0_up, z0_down, Imb_abs_0, s0, s_up0, s_down0 = extract_initial_values_from_benchmark(history[fold], fold_size)
            h = TA_bin_Julia(in_sample_Pw, length(inS), in_sample_prob, q_up_in, q_down_in, rf_up, rf_down, NI_up, NI_down, k0_up, k0_down, z0_up, z0_down, ppB0_up, ppB0_down, s_up0, s_down0, max_iter = max_iter, tol = 0.01, tol_m = 0.01, penalty_up = penalty_up, penalty_down = penalty_down)
            h_end_bin = h[end]
            h_LP = TA_bin_LP_Julia(h_end_bin, in_sample_Pw, length(inS), in_sample_prob, q_up_in, q_down_in, rf_up, rf_down, k0_up, k0_down, z0_up, z0_down, ppB0_up, ppB0_down, max_iter = max_iter, tol = 0.01)
            h_end_LP = h_LP[end]
 
            DA_bids_cen = [h_end_LP.ppDA, h_end_LP.pfDA, h_end_LP.pgDA]
            push!(history_insample, (h_end_bin, h_end_LP))
       
        elseif model == "TA_bin_GDCA"
            @unpack k0_up, k0_down, ppB0_up, ppB0_down, z0_up, z0_down, Imb_abs_0, s0, s_up0, s_down0 = extract_initial_values_from_benchmark(history[fold], fold_size)
            time_start = time()
            h = TA_bin_GDCA(in_sample_Pw, length(inS), in_sample_prob, q_up_in, q_down_in, rf_up, rf_down, NI_up, NI_down, k0_up, k0_down, z0_up, z0_down, ppB0_up, ppB0_down, s_up0, s_down0, penalty_t = penalty_t, max_iter = max_iter, tol = 0.01, alpha = 1.0)
            time_end = time()
            println("TA_bin_GDCA time: ", time_end - time_start, " seconds")
            h_end_bin = h[end]
            time_start = time()
            h_LP = TA_bin_LP_GDCA(h_end_bin, in_sample_Pw, length(inS), in_sample_prob, q_up_in, q_down_in, rf_up, rf_down, NI_up, NI_down, k0_up, k0_down, z0_up, z0_down, ppB0_up, ppB0_down, max_iter = max_iter, penalty_t = penalty_t, tol = 0.01, alpha=1.0)
            time_end = time()
            println("TA_bin_LP_GDCA time: ", time_end - time_start, " seconds")
            h_end = h_LP[end]
 
            DA_bids_cen = [h_end.ppDA, h_end.pfDA, h_end.pgDA]
            push!(history_insample, (h_end_bin, h_end))
 
        elseif model == "Benchmark"
            time_start = time()
            h = Benchmark(in_sample_Pw, length(inS), in_sample_prob, rf_up, rf_down, NI_up, NI_down)
            time_end = time()
            println("Benchmark time: ", time_end - time_start, " seconds")
            h_end = h[end]
            DA_bids_cen = [h_end.ppDA, h_end.pfDA, h_end.pgDA]
            push!(history_insample, h_end)
        else
            error("Model '$model' not implemented in the k-fold function.")
        end
       
        # Ensure DA_bids_cen is defined for the next steps
        if !@isdefined(DA_bids_cen)
            error("Model training failed to produce DA bids for fold $fold.")
        end
 
        # 2.3. Market Simulation
       
        # DA Market Price
        Lambda_DA, f = merit_order_price_DA(DA_bids_cen)
 
        # In-sample profits and bids
        Lambda_B_in,  B_bids_pol_cen_in, B_bids_pol_up_cen_in, B_bids_pol_down_cen_in, act_bids_up_cen_in, act_bids_down_cen_in = merit_order_price_B(DA_bids_cen, in_sample_Pw, length(inS), R_bids[1], R_bids[2], NI_up, NI_down)
        B_bids_cen_in = [B_bids_pol_cen_in, B_bids_pol_up_cen_in, B_bids_pol_down_cen_in, act_bids_up_cen_in, act_bids_down_cen_in]
        pol_profit_in, flex_profit_in, reg_profit_in, SI_in = profits(length(inS), in_sample_prob, DA_bids_cen, B_bids_cen_in, R_bids, Lambda_DA, Lambda_B_in, Lambda_R, Cf_Rup, Cf_Rdown; Q_up = q_up_in, Q_down = q_down_in)
 
        # Out-of-sample profits and bids
        Lambda_B_out,  B_bids_pol_cen_out, B_bids_pol_up_cen_out, B_bids_pol_down_cen_out, act_bids_up_cen_out, act_bids_down_cen_out = merit_order_price_B(DA_bids_cen, out_sample_Pw, onS, R_bids[1], R_bids[2], NI_up, NI_down)
        B_bids_cen_out = [B_bids_pol_cen_out, B_bids_pol_up_cen_out, B_bids_pol_down_cen_out, act_bids_up_cen_out, act_bids_down_cen_out]
        pol_profit_out, flex_profit_out, reg_profit_out, SI_out = profits(onS, out_sample_prob, DA_bids_cen, B_bids_cen_out, R_bids, Lambda_DA, Lambda_B_out, Lambda_R, Cf_Rup, Cf_Rdown; Q_up = q_up_in, Q_down = q_down_in)
 
        # Average Balancing Bids
        exp_pol_bids_in, exp_pol_bids_up_in, exp_pol_bids_down_in, exp_flex_bids_up_in, exp_flex_bids_down_in = average_bids_over_scenarios(in_sample_prob, B_bids_cen_in)
        exp_pol_bids_out, exp_pol_bids_up_out, exp_pol_bids_down_out, exp_flex_bids_up_out, exp_flex_bids_down_out = average_bids_over_scenarios(out_sample_prob, B_bids_cen_out)
 
        # 2.4. Collect Results
        push!(profits_outsample, (pol_profit_out, flex_profit_out, reg_profit_out))
        push!(profits_insample, (pol_profit_in, flex_profit_in, reg_profit_in))
        push!(B_bids_insample, (exp_pol_bids_in, exp_pol_bids_up_in, exp_pol_bids_down_in, exp_flex_bids_up_in, exp_flex_bids_down_in))
        push!(B_bids_outsample, (exp_pol_bids_out, exp_pol_bids_up_out, exp_pol_bids_down_out, exp_flex_bids_up_out, exp_flex_bids_down_out))
        push!(DA_bids_insample, DA_bids_cen)
        push!(DA_bids_outsample, DA_bids_cen)
        push!(lambdas, (Lambda_DA, Lambda_B_in, Lambda_B_out))
        push!(system_imbalance, (SI_in, SI_out))
        push!(splits_info, (in_sample_prob, out_sample_prob, in_sample_Pw, out_sample_Pw))
        push!(B_bids_all_insample, (B_bids_pol_up_cen_in, B_bids_pol_down_cen_in))
        push!(B_bids_all_outsample, (B_bids_pol_up_cen_out, B_bids_pol_down_cen_out))

        if fold == k 
            #display(f)  #DA merit order 
        end

        if bin 
            # Binary Optimization
            DA_bids_cen_bin = [h_end_bin.ppDA, h_end_bin.pfDA, h_end_bin.pgDA]
           
            # DA Market Price
            Lambda_DA_bin, _ = merit_order_price_DA(DA_bids_cen_bin)
 
            # Out-of-sample profits and bids
            Lambda_B_out_bin, B_bids_pol_cen_out_bin, B_bids_pol_up_cen_out_bin, B_bids_pol_down_cen_out_bin, act_bids_up_cen_out_bin, act_bids_down_cen_out_bin = merit_order_price_B(DA_bids_cen, out_sample_Pw, onS, R_bids[1], R_bids[2], NI_up, NI_down)
            B_bids_cen_out_bin = [B_bids_pol_cen_out_bin, B_bids_pol_up_cen_out_bin, B_bids_pol_down_cen_out_bin, act_bids_up_cen_out_bin, act_bids_down_cen_out_bin]
            pol_profit_out_bin, flex_profit_out_bin, reg_profit_out_bin, _ = profits(onS, out_sample_prob, DA_bids_cen_bin, B_bids_cen_out_bin, R_bids, Lambda_DA_bin, Lambda_B_out_bin, Lambda_R, Cf_Rup, Cf_Rdown; Q_up = q_up_in, Q_down = q_down_in)
            push!(profits_outsample_bin, (pol_profit_out_bin, flex_profit_out_bin, reg_profit_out_bin))
        end 
 
        if ind
            #----------------------- Binary ---------------------------#
            
            #----------------------- Individual Opt ---------------------------#
            DA_bids_ind = [Float64[] for _ in 1:3]
 
            # --- Polluters ---
            for p in 1:nP
                k0_up_p, k0_down_p, z0_up_p, z0_down_p, ppB0_up_p, ppB0_down_p = k0_up[:, p], k0_down[:, p], z0_up, z0_down, ppB0_up[:, p], ppB0_down[:, p]
                h_new_LP = drop_polluter(h_end, p)
                if model =="TA_bin_GDCA"
                    h_ind = polluter_bin_LP(in_sample_Pw, length(inS), in_sample_prob, p, h_new_LP, q_up_in, q_down_in, NI_up, NI_down, k0_up_p, k0_down_p, z0_up_p, z0_down_p, ppB0_up_p, ppB0_down_p, penalty_t = penalty_t, max_iter = max_iter, tol = 0.01, alpha = 1.0)
                elseif model == "Benchmark"
                    h_ind = polluter_benchmark(p, h_new_LP, in_sample_Pw, length(inS), in_sample_prob)
                end
                append!(DA_bids_ind[1], h_ind[end].ppDA)
            end
            # --- Flexible ---
            for f in 1:nF
                if model =="TA_bin_GDCA"
                    h_ind = flex_bin_LP(in_sample_Pw, length(inS), in_sample_prob, f, h_end, rf_up, rf_down)
                elseif model == "Benchmark"
                    h_ind = flex_benchmark(f, h_end, in_sample_Pw, length(inS), in_sample_prob, rf_up, rf_down)
                end
                append!(DA_bids_ind[2], h_ind[end].pfDA)
            end
            # --- Regular ---
            for r in 1:nR
                h_ind = reg_opt(r, h_end)
                append!(DA_bids_ind[3], h_ind[end].pgDA)
            end
           
            # DA Market Price
            Lambda_DA_ind, _ = merit_order_price_DA(DA_bids_ind)
 
            # Out-of-sample profits and bids
            Lambda_B_out_ind, B_bids_pol_out_ind, B_bids_pol_up_out_ind, B_bids_pol_down_out_ind, act_bids_up_out_ind, act_bids_down_out_ind = merit_order_price_B(DA_bids_ind, out_sample_Pw, onS, R_bids[1], R_bids[2], NI_up, NI_down)
            B_bids_out_ind = [B_bids_pol_out_ind, B_bids_pol_up_out_ind, B_bids_pol_down_out_ind, act_bids_up_out_ind, act_bids_down_out_ind]
            pol_profit_out_ind, flex_profit_out_ind, reg_profit_out_ind, _ = profits(onS, out_sample_prob, DA_bids_ind, B_bids_out_ind, R_bids, Lambda_DA_ind, Lambda_B_out_ind, Lambda_R, Cf_Rup, Cf_Rdown; Q_up = q_up_in, Q_down = q_down_in)
            push!(profits_outsample_ind, (pol_profit_out_ind, flex_profit_out_ind, reg_profit_out_ind))
        end
    end
 
    # --- 3. Aggregate and Post-Process Results ---
    pol_matrix_in  = hcat([p[1] for p in profits_insample]...)
    flex_matrix_in = hcat([p[2] for p in profits_insample]...)
    reg_matrix_in  = hcat([p[3] for p in profits_insample]...)
 
    pol_matrix_out  = hcat([p[1] for p in profits_outsample]...)
    flex_matrix_out = hcat([p[2] for p in profits_outsample]...)
    reg_matrix_out  = hcat([p[3] for p in profits_outsample]...)
 
    B_bids_pol_in = hcat([p[1] for p in B_bids_insample]...)
    B_bids_pol_up_in = hcat([p[2] for p in B_bids_insample]...)
    B_bids_pol_down_in = hcat([p[3] for p in B_bids_insample]...)
    B_bids_flex_up_in = hcat([p[4] for p in B_bids_insample]...)
    B_bids_flex_down_in = hcat([p[5] for p in B_bids_insample]...)
 
    B_bids_pol_out = hcat([p[1] for p in B_bids_outsample]...)
    B_bids_pol_up_out = hcat([p[2] for p in B_bids_outsample]...)
    B_bids_pol_down_out = hcat([p[3] for p in B_bids_outsample]...)
    B_bids_flex_up_out = hcat([p[4] for p in B_bids_outsample]...)
    B_bids_flex_down_out = hcat([p[5] for p in B_bids_outsample]...)
 
    DA_bids_pol_out= hcat([p[1] for p in  DA_bids_outsample]...)   # 10 × k
    DA_bids_flex_out = hcat([p[2] for p in DA_bids_outsample]...)   # 10 × k
    DA_bids_reg_out = hcat([p[3] for p in DA_bids_outsample]...)   # 10 × k
 
    pol_out_avg = expected_revenues(pol_matrix_out)
    flex_out_avg = expected_revenues(flex_matrix_out)
    reg_out_avg = expected_revenues(reg_matrix_out)
 
    pol_in_sums = sum_over_polluters(pol_matrix_in)
    pol_out_sums = sum_over_polluters(pol_matrix_out)
    flex_in_sums = sum_over_polluters(flex_matrix_in)
    flex_out_sums = sum_over_polluters(flex_matrix_out)
    reg_in_sums = sum_over_polluters(reg_matrix_in)
    reg_out_sums = sum_over_polluters(reg_matrix_out)
 
    mean_pol_bids_in = mean(B_bids_pol_in, dims = 2)
    mean_pol_bids_out = mean(B_bids_pol_out, dims = 2)
    mean_pol_bids_up_in = mean(B_bids_pol_up_in, dims = 2)
    mean_pol_bids_up_out = mean(B_bids_pol_up_out, dims = 2)
    mean_pol_bids_down_in = mean(B_bids_pol_down_in, dims = 2)
    mean_pol_bids_down_out = mean(B_bids_pol_down_out, dims = 2)
    mean_flex_bids_up_in = mean(B_bids_flex_up_in, dims = 2)
    mean_flex_bids_up_out = mean(B_bids_flex_up_out, dims = 2)
    mean_flex_bids_down_in = mean(B_bids_flex_down_in, dims = 2)
    mean_flex_bids_down_out = mean(B_bids_flex_down_out, dims = 2)
 
    mean_pol_DA_bids_out = mean(DA_bids_pol_out, dims = 2)
    mean_flex_DA_bids_out = mean(DA_bids_flex_out, dims = 2)
    mean_reg_DA_bids_out = mean(DA_bids_reg_out, dims = 2)

    if bin 
        pol_matrix_out_bin  = hcat([p[1] for p in profits_outsample_bin]...)
        flex_matrix_out_bin = hcat([p[2] for p in profits_outsample_bin]...)
        reg_matrix_out_bin  = hcat([p[3] for p in profits_outsample_bin]...)

        pol_out_avg_bin = expected_revenues(pol_matrix_out_bin)
        flex_out_avg_bin = expected_revenues(flex_matrix_out_bin)
        reg_out_avg_bin = expected_revenues(reg_matrix_out_bin)

    end 
 
    if ind
        pol_matrix_out_ind  = hcat([p[1] for p in profits_outsample_ind]...)
        flex_matrix_out_ind = hcat([p[2] for p in profits_outsample_ind]...)
        reg_matrix_out_ind  = hcat([p[3] for p in profits_outsample_ind]...)
 
        pol_out_avg_ind = expected_revenues(pol_matrix_out_ind)
        flex_out_avg_ind = expected_revenues(flex_matrix_out_ind)
        reg_out_avg_ind = expected_revenues(reg_matrix_out_ind)
    end
 
    # --- 4. Return all required variables ---
 
    profits_results = (
        pol_in  = pol_matrix_in,
        flex_in = flex_matrix_in,
        reg_in  = reg_matrix_in,
        pol_out  = pol_matrix_out,
        flex_out = flex_matrix_out,
        reg_out  = reg_matrix_out,
    )
 
    bids = (
        B_pol_in        = B_bids_pol_in,
        B_pol_up_in     = B_bids_pol_up_in,
        B_pol_down_in   = B_bids_pol_down_in,
        B_flex_up_in    = B_bids_flex_up_in,
        B_flex_down_in  = B_bids_flex_down_in,
        B_pol_out       = B_bids_pol_out,
        B_pol_up_out    = B_bids_pol_up_out,
        B_pol_down_out  = B_bids_pol_down_out,
        B_flex_up_out   = B_bids_flex_up_out,
        B_flex_down_out = B_bids_flex_down_out,
        DA_bids_insample = DA_bids_insample,
        DA_bids_outsample = DA_bids_outsample,
        B_bids_all_insample = B_bids_all_insample,
        B_bids_all_outsample = B_bids_all_outsample,
    )
 
    aggregates = (
        pol_out_avg = pol_out_avg,
        flex_out_avg = flex_out_avg,
        reg_out_avg  = reg_out_avg,
        pol_in_sums  = pol_in_sums,
        pol_out_sums = pol_out_sums,
        flex_in_sums = flex_in_sums,
        flex_out_sums = flex_out_sums,
        reg_in_sums = reg_in_sums,
        reg_out_sums = reg_out_sums,
        mean_pol_bids_in         = mean_pol_bids_in,
        mean_pol_bids_out        = mean_pol_bids_out,
        mean_pol_bids_up_in      = mean_pol_bids_up_in,
        mean_pol_bids_up_out     = mean_pol_bids_up_out,
        mean_pol_bids_down_in    = mean_pol_bids_down_in,
        mean_pol_bids_down_out   = mean_pol_bids_down_out,
        mean_flex_bids_up_in     = mean_flex_bids_up_in,
        mean_flex_bids_up_out    = mean_flex_bids_up_out,
        mean_flex_bids_down_in   = mean_flex_bids_down_in,
        mean_flex_bids_down_out  = mean_flex_bids_down_out,
        mean_pol_DA_bids_out   = mean_pol_DA_bids_out,
        mean_flex_DA_bids_out  = mean_flex_DA_bids_out,
        mean_reg_DA_bids_out   = mean_reg_DA_bids_out,
    )

    
    split_info = (
        in_sample_prob = [s[1] for s in splits_info],
        out_sample_prob = [s[2] for s in splits_info],
        in_sample_Pw = [s[3] for s in splits_info],
        out_sample_Pw = [s[4] for s in splits_info],
    )
   
    if ind
        profits_results = merge(profits_results, (
            pol_out_ind  = pol_matrix_out_ind,
            flex_out_ind = flex_matrix_out_ind,
            reg_out_ind  = reg_matrix_out_ind,
        ))
 
        aggregates = merge(aggregates, (
            pol_out_avg_ind  = pol_out_avg_ind,
            flex_out_avg_ind = flex_out_avg_ind,
            reg_out_avg_ind  = reg_out_avg_ind,
        ))
    end
    if bin
        profits_results = merge(profits_results, (
            pol_out_bin  = pol_matrix_out_bin,
            flex_out_bin = flex_matrix_out_bin,
            reg_out_bin  = reg_matrix_out_bin,
        ))
 
        aggregates = merge(aggregates, (
            pol_out_avg_bin  = pol_out_avg_bin,
            flex_out_avg_bin = flex_out_avg_bin,
            reg_out_avg_bin  = reg_out_avg_bin,
        )) 

    end 

 
    history = (
        lambdas = lambdas,
        history_insample = history_insample,
        system_imbalance = system_imbalance,
    )
 
    return (profits = profits_results, bids = bids, aggregates = aggregates, history = history, split_info = split_info)
end



function run_kfold_simulation_chosen_folds( model::String, Pw, nS_all::Any, prob::Any,
    q_up_in::Any,  q_down_in::Any, k::Any, fold_size::Int,
    rf_up::Any, rf_down::Any, NI_up::Any, NI_down::Any, Cf_Rup::Any, Cf_Rdown::Any, R_bids, Lambda_R, hour;
     history = nothing, max_iter = 100, penalty_up = 100, penalty_down = 100, alpha = 1.0, tol = 0.01,
      penalty_t = 100, ind::Bool = false, bin::Bool = false)
    # --- 1. Initialize result vectors ---
    profits_outsample = []
    profits_insample = []
    B_bids_insample = []
    B_bids_outsample = []
    DA_bids_insample = []
    DA_bids_outsample = []
    history_insample = [] # Correctly initialized
    lambdas = []
    system_imbalance = []
    profits_outsample_ind = []
    profits_outsample_bin = []
    splits_info = []
    B_bids_all_insample = []
    B_bids_all_outsample = []
 
    # --- 2. K-Fold Loop ---
    for (i, fold) in enumerate(k)
        println("Fold: ", fold)
        println("Using history from fold index: ", i)
   
        # 2.1. Determine indices and initialize values
        inS = ((fold - 1) * fold_size + 1):(fold * fold_size)
        in_sample_Pw, out_sample_Pw, in_sample_prob, out_sample_prob, onS = sample_split(Pw[hour], nS_all, prob, fold_size, fold)
 
        # 2.2. Model Training (Conditional on 'model' type)
        local h_end, h_end_LP, DA_bids_cen
       
        if model == "TA_bin_Julia"
            @unpack k0_up, k0_down, ppB0_up, ppB0_down, z0_up, z0_down, Imb_abs_0, s0, s_up0, s_down0 = extract_initial_values_from_benchmark(history[i], fold_size)
            h = TA_bin_Julia(in_sample_Pw, length(inS), in_sample_prob, q_up_in, q_down_in, rf_up, rf_down, NI_up, NI_down, k0_up, k0_down, z0_up, z0_down, ppB0_up, ppB0_down, s_up0, s_down0, max_iter = max_iter, tol = 0.01, tol_m = 0.01, penalty_up = penalty_up, penalty_down = penalty_down)
            h_end_bin = h[end]
            h_LP = TA_bin_LP_Julia(h_end_bin, in_sample_Pw, length(inS), in_sample_prob, q_up_in, q_down_in, rf_up, rf_down, k0_up, k0_down, z0_up, z0_down, ppB0_up, ppB0_down, max_iter = max_iter, tol = 0.01)
            h_end_LP = h_LP[end]
 
            DA_bids_cen = [h_end_LP.ppDA, h_end_LP.pfDA, h_end_LP.pgDA]
            push!(history_insample, (h_end_bin, h_end_LP))
       
        elseif model == "TA_bin_GDCA"
            #@unpack k0_up, k0_down, ppB0_up, ppB0_down, z0_up, z0_down, Imb_abs_0, s0, s_up0, s_down0 = extract_initial_values_from_benchmark(history[i], fold_size)
            k0_up, k0_down, ppB0_up, ppB0_down, z0_up, z0_down, Imb_abs_0, s0, s_up0, s_down0 = initialize_values(inS)
            h = TA_bin_GDCA2(in_sample_Pw, length(inS), in_sample_prob, q_up_in, q_down_in, rf_up, rf_down, NI_up, NI_down, k0_up, k0_down, z0_up, z0_down, ppB0_up, ppB0_down, s_up0, s_down0, penalty_t = penalty_t, max_iter = max_iter, tol = tol, alpha = alpha)
            h_end_bin = h[end]
            h_LP = TA_bin_LP_GDCA(h_end_bin, in_sample_Pw, length(inS), in_sample_prob, q_up_in, q_down_in, rf_up, rf_down, NI_up, NI_down, k0_up, k0_down, z0_up, z0_down, ppB0_up, ppB0_down, max_iter = max_iter, penalty_t = penalty_t, tol = tol, alpha=alpha)
            h_end = h_LP[end]
 
            DA_bids_cen = [h_end.ppDA, h_end.pfDA, h_end.pgDA]
            push!(history_insample, (h_end_bin, h_end))
 
        elseif model == "Benchmark"
            h = Benchmark(in_sample_Pw, length(inS), in_sample_prob, rf_up, rf_down, NI_up, NI_down)
            h_end = h[end]
            DA_bids_cen = [h_end.ppDA, h_end.pfDA, h_end.pgDA]
            push!(history_insample, h_end)
        else
            error("Model '$model' not implemented in the k-fold function.")
        end
       
        # Ensure DA_bids_cen is defined for the next steps
        if !@isdefined(DA_bids_cen)
            error("Model training failed to produce DA bids for fold $fold.")
        end
 
        # 2.3. Market Simulation
       
        # DA Market Price
        Lambda_DA, f = merit_order_price_DA(DA_bids_cen)
 
        # In-sample profits and bids
        Lambda_B_in,  B_bids_pol_cen_in, B_bids_pol_up_cen_in, B_bids_pol_down_cen_in, act_bids_up_cen_in, act_bids_down_cen_in = merit_order_price_B(DA_bids_cen, in_sample_Pw, length(inS), R_bids[1], R_bids[2], NI_up, NI_down)
        #display(B_bids_pol_cen_in)
        #display(B_bids_pol_up_cen_in)
        #display(B_bids_pol_down_cen_in)
        B_bids_cen_in = [B_bids_pol_cen_in, B_bids_pol_up_cen_in, B_bids_pol_down_cen_in, act_bids_up_cen_in, act_bids_down_cen_in]
        pol_profit_in, flex_profit_in, reg_profit_in, SI_in = profits(length(inS), in_sample_prob, DA_bids_cen, B_bids_cen_in, R_bids, Lambda_DA, Lambda_B_in, Lambda_R, Cf_Rup, Cf_Rdown; Q_up = q_up_in, Q_down = q_down_in)
 
        # Out-of-sample profits and bids
        Lambda_B_out,  B_bids_pol_cen_out, B_bids_pol_up_cen_out, B_bids_pol_down_cen_out, act_bids_up_cen_out, act_bids_down_cen_out = merit_order_price_B(DA_bids_cen, out_sample_Pw, onS, R_bids[1], R_bids[2], NI_up, NI_down)
        #display("Out-of-sample Bids:")
        #display(B_bids_pol_cen_out)
        #display(B_bids_pol_up_cen_out)
        #display(B_bids_pol_down_cen_out)
        
        B_bids_cen_out = [B_bids_pol_cen_out, B_bids_pol_up_cen_out, B_bids_pol_down_cen_out, act_bids_up_cen_out, act_bids_down_cen_out]
        pol_profit_out, flex_profit_out, reg_profit_out, SI_out = profits(onS, out_sample_prob, DA_bids_cen, B_bids_cen_out, R_bids, Lambda_DA, Lambda_B_out, Lambda_R, Cf_Rup, Cf_Rdown; Q_up = q_up_in, Q_down = q_down_in)
 
        # Average Balancing Bids
        exp_pol_bids_in, exp_pol_bids_up_in, exp_pol_bids_down_in, exp_flex_bids_up_in, exp_flex_bids_down_in = average_bids_over_scenarios(in_sample_prob, B_bids_cen_in)
        exp_pol_bids_out, exp_pol_bids_up_out, exp_pol_bids_down_out, exp_flex_bids_up_out, exp_flex_bids_down_out = average_bids_over_scenarios(out_sample_prob, B_bids_cen_out)
        #display("eXPECTED out-sample Bids:")
        #display(exp_pol_bids_out)
        #display(exp_pol_bids_up_out)
        #display(exp_pol_bids_down_out)
        # 2.4. Collect Results
        push!(profits_outsample, (pol_profit_out, flex_profit_out, reg_profit_out))
        push!(profits_insample, (pol_profit_in, flex_profit_in, reg_profit_in))
        push!(B_bids_insample, (exp_pol_bids_in, exp_pol_bids_up_in, exp_pol_bids_down_in, exp_flex_bids_up_in, exp_flex_bids_down_in))
        push!(B_bids_outsample, (exp_pol_bids_out, exp_pol_bids_up_out, exp_pol_bids_down_out, exp_flex_bids_up_out, exp_flex_bids_down_out))
        push!(DA_bids_insample, DA_bids_cen)
        push!(DA_bids_outsample, DA_bids_cen)
        push!(lambdas, (Lambda_DA, Lambda_B_in, Lambda_B_out))
        push!(system_imbalance, (SI_in, SI_out))
        push!(splits_info, (in_sample_prob, out_sample_prob, in_sample_Pw, out_sample_Pw))
        push!(B_bids_all_insample, (B_bids_pol_up_cen_in, B_bids_pol_down_cen_in))
        push!(B_bids_all_outsample, (B_bids_pol_up_cen_out, B_bids_pol_down_cen_out))

     
        if bin 
            # Binary Optimization
            DA_bids_cen_bin = [h_end_bin.ppDA, h_end_bin.pfDA, h_end_bin.pgDA]
           
            # DA Market Price
            Lambda_DA_bin, _ = merit_order_price_DA(DA_bids_cen_bin)
 
            # Out-of-sample profits and bids
            Lambda_B_out_bin,  B_bids_pol_cen_out_bin, B_bids_pol_up_cen_out_bin, B_bids_pol_down_cen_out_bin, act_bids_up_cen_out_bin, act_bids_down_cen_out_bin = merit_order_price_B(DA_bids_cen, out_sample_Pw, onS, R_bids[1], R_bids[2], NI_up, NI_down)
            B_bids_cen_out_bin = [B_bids_pol_cen_out_bin, B_bids_pol_up_cen_out_bin, B_bids_pol_down_cen_out_bin, act_bids_up_cen_out_bin, act_bids_down_cen_out_bin]
            pol_profit_out_bin, flex_profit_out_bin, reg_profit_out_bin, _ = profits(onS, out_sample_prob, DA_bids_cen_bin, B_bids_cen_out_bin, R_bids, Lambda_DA_bin, Lambda_B_out_bin, Lambda_R, Cf_Rup, Cf_Rdown; Q_up = q_up_in, Q_down = q_down_in)
            push!(profits_outsample_bin, (pol_profit_out_bin, flex_profit_out_bin, reg_profit_out_bin))
        end 
 
        if ind
            #----------------------- Binary ---------------------------#
            
            #----------------------- Individual Opt ---------------------------#
            DA_bids_ind = [Float64[] for _ in 1:3]
 
            # --- Polluters ---
            for p in 1:nP
                k0_up_p, k0_down_p, z0_up_p, z0_down_p, ppB0_up_p, ppB0_down_p = k0_up[:, p], k0_down[:, p], z0_up, z0_down, ppB0_up[:, p], ppB0_down[:, p]
                h_new_LP = drop_polluter(h_end, p)
                if model =="TA_bin_GDCA"
                    h_ind = polluter_bin_LP(in_sample_Pw, length(inS), in_sample_prob, p, h_new_LP, q_up_in, q_down_in, NI_up, NI_down, k0_up_p, k0_down_p, z0_up_p, z0_down_p, ppB0_up_p, ppB0_down_p, penalty_t = penalty_t, max_iter = max_iter, tol = 0.01, alpha = 1.0)
                elseif model == "Benchmark"
                    h_ind = polluter_benchmark(p, h_new_LP, in_sample_Pw, length(inS), in_sample_prob)
                end
                append!(DA_bids_ind[1], h_ind[end].ppDA)
            end
            # --- Flexible ---
            for f in 1:nF
                if model =="TA_bin_GDCA"
                    h_ind = flex_bin_LP(in_sample_Pw, length(inS), in_sample_prob, f, h_end, rf_up, rf_down)
                elseif model == "Benchmark"
                    h_ind = flex_benchmark(f, h_end, in_sample_Pw, length(inS), in_sample_prob, rf_up, rf_down)
                end
                append!(DA_bids_ind[2], h_ind[end].pfDA)
            end
            # --- Regular ---
            for r in 1:nR
                h_ind = reg_opt(r, h_end)
                append!(DA_bids_ind[3], h_ind[end].pgDA)
            end
           
            # DA Market Price
            Lambda_DA_ind, _ = merit_order_price_DA(DA_bids_ind)
 
            # Out-of-sample profits and bids
            Lambda_B_out_ind,  B_bids_pol_out_ind, B_bids_pol_up_out_ind, B_bids_pol_down_out_ind, act_bids_up_out_ind, act_bids_down_out_ind = merit_order_price_B(DA_bids_ind, out_sample_Pw, onS, R_bids[1], R_bids[2], NI_up, NI_down)
            B_bids_out_ind = [B_bids_pol_out_ind, B_bids_pol_up_out_ind, B_bids_pol_down_out_ind, act_bids_up_out_ind, act_bids_down_out_ind]
            pol_profit_out_ind, flex_profit_out_ind, reg_profit_out_ind, _ = profits(onS, out_sample_prob, DA_bids_ind, B_bids_out_ind, R_bids, Lambda_DA_ind, Lambda_B_out_ind, Lambda_R, Cf_Rup, Cf_Rdown; Q_up = q_up_in, Q_down = q_down_in)
            push!(profits_outsample_ind, (pol_profit_out_ind, flex_profit_out_ind, reg_profit_out_ind))
        end
    end
 
    # --- 3. Aggregate and Post-Process Results ---
    pol_matrix_in  = hcat([p[1] for p in profits_insample]...)
    flex_matrix_in = hcat([p[2] for p in profits_insample]...)
    reg_matrix_in  = hcat([p[3] for p in profits_insample]...)
 
    pol_matrix_out  = hcat([p[1] for p in profits_outsample]...)
    flex_matrix_out = hcat([p[2] for p in profits_outsample]...)
    reg_matrix_out  = hcat([p[3] for p in profits_outsample]...)
 
    B_bids_pol_in = hcat([p[1] for p in B_bids_insample]...)
    B_bids_pol_up_in = hcat([p[2] for p in B_bids_insample]...)
    B_bids_pol_down_in = hcat([p[3] for p in B_bids_insample]...)
    B_bids_flex_up_in = hcat([p[4] for p in B_bids_insample]...)
    B_bids_flex_down_in = hcat([p[5] for p in B_bids_insample]...)
 
    B_bids_pol_out = hcat([p[1] for p in B_bids_outsample]...)
    B_bids_pol_up_out = hcat([p[2] for p in B_bids_outsample]...)
    B_bids_pol_down_out = hcat([p[3] for p in B_bids_outsample]...)
    B_bids_flex_up_out = hcat([p[4] for p in B_bids_outsample]...)
    B_bids_flex_down_out = hcat([p[5] for p in B_bids_outsample]...)
 
    DA_bids_pol_out= hcat([p[1] for p in  DA_bids_outsample]...)   # 10 × k
    DA_bids_flex_out = hcat([p[2] for p in DA_bids_outsample]...)   # 10 × k
    DA_bids_reg_out = hcat([p[3] for p in DA_bids_outsample]...)   # 10 × k
 
    pol_out_avg = expected_revenues(pol_matrix_out)
    flex_out_avg = expected_revenues(flex_matrix_out)
    reg_out_avg = expected_revenues(reg_matrix_out)
 
    pol_in_sums = sum_over_polluters(pol_matrix_in)
    pol_out_sums = sum_over_polluters(pol_matrix_out)
    flex_in_sums = sum_over_polluters(flex_matrix_in)
    flex_out_sums = sum_over_polluters(flex_matrix_out)
    reg_in_sums = sum_over_polluters(reg_matrix_in)
    reg_out_sums = sum_over_polluters(reg_matrix_out)
 
    mean_pol_bids_in = mean(B_bids_pol_in, dims = 2)
    mean_pol_bids_out = mean(B_bids_pol_out, dims = 2)
    mean_pol_bids_up_in = mean(B_bids_pol_up_in, dims = 2)
    mean_pol_bids_up_out = mean(B_bids_pol_up_out, dims = 2)
    mean_pol_bids_down_in = mean(B_bids_pol_down_in, dims = 2)
    mean_pol_bids_down_out = mean(B_bids_pol_down_out, dims = 2)
    mean_flex_bids_up_in = mean(B_bids_flex_up_in, dims = 2)
    mean_flex_bids_up_out = mean(B_bids_flex_up_out, dims = 2)
    mean_flex_bids_down_in = mean(B_bids_flex_down_in, dims = 2)
    mean_flex_bids_down_out = mean(B_bids_flex_down_out, dims = 2)
 
    mean_pol_DA_bids_out = mean(DA_bids_pol_out, dims = 2)
    mean_flex_DA_bids_out = mean(DA_bids_flex_out, dims = 2)
    mean_reg_DA_bids_out = mean(DA_bids_reg_out, dims = 2)

    if bin 
        pol_matrix_out_bin  = hcat([p[1] for p in profits_outsample_bin]...)
        flex_matrix_out_bin = hcat([p[2] for p in profits_outsample_bin]...)
        reg_matrix_out_bin  = hcat([p[3] for p in profits_outsample_bin]...)

        pol_out_avg_bin = expected_revenues(pol_matrix_out_bin)
        flex_out_avg_bin = expected_revenues(flex_matrix_out_bin)
        reg_out_avg_bin = expected_revenues(reg_matrix_out_bin)

    end 
 
    if ind
        pol_matrix_out_ind  = hcat([p[1] for p in profits_outsample_ind]...)
        flex_matrix_out_ind = hcat([p[2] for p in profits_outsample_ind]...)
        reg_matrix_out_ind  = hcat([p[3] for p in profits_outsample_ind]...)
 
        pol_out_avg_ind = expected_revenues(pol_matrix_out_ind)
        flex_out_avg_ind = expected_revenues(flex_matrix_out_ind)
        reg_out_avg_ind = expected_revenues(reg_matrix_out_ind)
    end
 
    # --- 4. Return all required variables ---
 
    profits_results = (
        pol_in  = pol_matrix_in,
        flex_in = flex_matrix_in,
        reg_in  = reg_matrix_in,
        pol_out  = pol_matrix_out,
        flex_out = flex_matrix_out,
        reg_out  = reg_matrix_out,
    )
 
    bids = (
        B_pol_in        = B_bids_pol_in,
        B_pol_up_in     = B_bids_pol_up_in,
        B_pol_down_in   = B_bids_pol_down_in,
        B_flex_up_in    = B_bids_flex_up_in,
        B_flex_down_in  = B_bids_flex_down_in,
        B_pol_out       = B_bids_pol_out,
        B_pol_up_out    = B_bids_pol_up_out,
        B_pol_down_out  = B_bids_pol_down_out,
        B_flex_up_out   = B_bids_flex_up_out,
        B_flex_down_out = B_bids_flex_down_out,
        DA_bids_insample = DA_bids_insample,
        DA_bids_outsample = DA_bids_outsample,
        B_bids_all_insample = B_bids_all_insample,
        B_bids_all_outsample = B_bids_all_outsample,
    )
 
    aggregates = (
        pol_out_avg = pol_out_avg,
        flex_out_avg = flex_out_avg,
        reg_out_avg  = reg_out_avg,
        pol_in_sums  = pol_in_sums,
        pol_out_sums = pol_out_sums,
        flex_in_sums = flex_in_sums,
        flex_out_sums = flex_out_sums,
        reg_in_sums = reg_in_sums,
        reg_out_sums = reg_out_sums,
        mean_pol_bids_in         = mean_pol_bids_in,
        mean_pol_bids_out        = mean_pol_bids_out,
        mean_pol_bids_up_in      = mean_pol_bids_up_in,
        mean_pol_bids_up_out     = mean_pol_bids_up_out,
        mean_pol_bids_down_in    = mean_pol_bids_down_in,
        mean_pol_bids_down_out   = mean_pol_bids_down_out,
        mean_flex_bids_up_in     = mean_flex_bids_up_in,
        mean_flex_bids_up_out    = mean_flex_bids_up_out,
        mean_flex_bids_down_in   = mean_flex_bids_down_in,
        mean_flex_bids_down_out  = mean_flex_bids_down_out,
        mean_pol_DA_bids_out   = mean_pol_DA_bids_out,
        mean_flex_DA_bids_out  = mean_flex_DA_bids_out,
        mean_reg_DA_bids_out   = mean_reg_DA_bids_out,
    )

    
    split_info = (
        in_sample_prob = [s[1] for s in splits_info],
        out_sample_prob = [s[2] for s in splits_info],
        in_sample_Pw = [s[3] for s in splits_info],
        out_sample_Pw = [s[4] for s in splits_info],
    )
   
    if ind
        profits_results = merge(profits_results, (
            pol_out_ind  = pol_matrix_out_ind,
            flex_out_ind = flex_matrix_out_ind,
            reg_out_ind  = reg_matrix_out_ind,
        ))
 
        aggregates = merge(aggregates, (
            pol_out_avg_ind  = pol_out_avg_ind,
            flex_out_avg_ind = flex_out_avg_ind,
            reg_out_avg_ind  = reg_out_avg_ind,
        ))
    end
    if bin
        profits_results = merge(profits_results, (
            pol_out_bin  = pol_matrix_out_bin,
            flex_out_bin = flex_matrix_out_bin,
            reg_out_bin  = reg_matrix_out_bin,
        ))
 
        aggregates = merge(aggregates, (
            pol_out_avg_bin  = pol_out_avg_bin,
            flex_out_avg_bin = flex_out_avg_bin,
            reg_out_avg_bin  = reg_out_avg_bin,
        )) 

    end 

 
    history = (
        lambdas = lambdas,
        history_insample = history_insample,
        system_imbalance = system_imbalance,
    )
 
    return (profits = profits_results, bids = bids, aggregates = aggregates, history = history, split_info = split_info)
end



function namedtuple_to_dict(x)
    if x isa NamedTuple
        return Dict(k => namedtuple_to_dict(v) for (k,v) in pairs(x))
    elseif x isa AbstractArray
        return [namedtuple_to_dict(v) for v in x]
    else
        return x
    end
end

####################################################################################

                        # MERIT ORDER CURVES

function merit_order_price_DA(DA_bids)
    # --- Collect bids and costs ---

    bids = vcat(
        [(b, Cp[i], :polluter, "P$(i)") for (i, b) in enumerate(DA_bids[1])],
        [(b, Cf[i], :flexible, "F$(i)") for (i, b) in enumerate(DA_bids[2])],
        [(b, Cr[i], :regular,  "R$(i)") for (i, b) in enumerate(DA_bids[3])]
    )

    # --- Sort by cost (ascending merit order) ---
    sorted_bids = sort(bids, by = x -> x[2])

    # --- Compute total demand ---
    total_demand = sum(DA_bids[1]) + sum(DA_bids[2]) + sum(DA_bids[3])

    # --- Initialize plot ---

    colors = Dict(:polluter => :red, :flexible => :blue, :regular => :green)

    fig = plot(
        xlabel = "Aggregated generation (MW)",
        ylabel = "Cost (€/MWh)",
        legend = :bottomright,
        title = "Merit Order Curve for Day-Ahead Market",
        lw = 2
    )

    # --- Draw each generator block ---
    cumulative = 0.0
    last_cost = nothing
    counts_by_cost = Dict{Float64, Int}()



    for (bid, cost, type, name) in sorted_bids
        color = colors[type]
    

        if !isnothing(last_cost) && cost == last_cost && cumulative > 0
            plot!([cumulative, cumulative], [cost - 0.2, cost + 0.2],
                  color = color, lw = 1, label = "", linealpha=1) # A subtle, slightly thinner line
        end
        
        # vertical line (cost step) if not first generator
        if !isnothing(last_cost)
            plot!([cumulative, cumulative], [last_cost, cost],
                  color = :gray, lw = 1, label = "")
        end
        
        plot!([cumulative, cumulative + bid], [cost, cost],
              lw = 2, color = color, label = "")
        # add label at the midpoint of the block
        counts_by_cost[cost] = get(counts_by_cost, cost, 0) + 1

        if bid > 0
            annotate!(cumulative + counts_by_cost[cost]*20, cost + 0.5, text(name, 8, color))
        end
        cumulative += bid
        last_cost = cost
    end

    # --- Find and mark market clearing price ---
    market_price = 0.0
    cumulative = 0.0
     for (bid, cost, _, _) in sorted_bids
        cumulative += bid
        if round(cumulative, digits = 3) >= round(total_demand, digits = 3)
            #println("cum", round(cumulative, digits = 3))
            #println("dem", round(total_demand, digits = 3))

            market_price = cost
            break
        end
    end
    vline!([total_demand], label = "Total demand", ls = :dash, color = :black)
    hline!([market_price], label = "Market price", ls = :dash, color = :gray)

    return market_price, fig
end


function merit_order_price_B(DA_bids, Pw, nS, rf_up, rf_down, NI_up, NI_down; tol=1e-6)
    # nP: number of polluters, nF: number of flexible generators
    nP = size(Pw, 2)
    nF = length(rf_up)
    
    # Initialize output arrays
    prices = zeros(nS)
    activated_bids_up = zeros(nS, nF)
    activated_bids_down = zeros(nS, nF)
    ppB = zeros(nS, nP)
    ppB_up_contr = zeros(nS, nP)
    ppB_down_contr = zeros(nS, nP)

    # Map flexible names to indices for consistent data structure access
    flex_names = ["F$(i)" for i in 1:nF]
    name2idx = Dict(name => i for (i, name) in enumerate(flex_names))

    for w in 1:nS
        # 1. Calculate individual and total imbalance (Realization - DA Bid)
        for p in 1:nP
            ppB[w,p] = Pw[w,p] - DA_bids[1][p]
        end
        total_req = sum(ppB[w, :]) 
        abs_req = abs(total_req)

        # 2. Assign contribution to imbalance (only if it worsens the system state)
        for p in 1:nP
            imb = ppB[w, p]
            if total_req > tol && imb > 0
                ppB_down_contr[w, p] = imb # System has surplus, polluter adds to it
            elseif total_req < -tol && imb < 0
                ppB_up_contr[w, p] = imb   # System has deficit, polluter adds to it
            end
        end

        # 3. Determine regulation type and identify available market bids
        if total_req > tol
            # Surplus: Down regulation needed (reduction of production)
            selected_bids = [(rf_down[i], Cf[i], flex_names[i]) for i in 1:nF if rf_down[i] > 0]
            NI_limit = NI_down
            regulation_type = "down"
        elseif total_req < -tol
            # Deficit: Up regulation needed (increase of production)
            selected_bids = [(rf_up[i], Cf[i], flex_names[i]) for i in 1:nF if rf_up[i] > 0]
            NI_limit = NI_up
            regulation_type = "up"
        else
            selected_bids = []
            NI_limit = 0.0
            regulation_type = "none"
        end

        # 4. Sort bids by cost to form the Merit Order Curve
        sorted_bids = sort(selected_bids, by = x -> x[2])
        activated_plotting_bids = []

        if abs_req <= NI_limit + tol
            # SCENARIO A: Within market reserve limits (Standard Merit Order)
            cumulative = 0.0
            for (bid_cap, cost, name) in sorted_bids
                f = name2idx[name]
                remaining = abs_req - cumulative
                amount = min(bid_cap, max(0.0, remaining))
                
                if amount > tol
                    if regulation_type == "down"
                        activated_bids_down[w, f] = amount
                    else
                        activated_bids_up[w, f] = amount
                    end
                    push!(activated_plotting_bids, (amount, cost, name))
                    cumulative += amount
                end
                if cumulative >= abs_req - tol break end
            end
        else
            # SCENARIO B: Imbalance exceeds market bids (|total_req| > NI)
            # Step 1: Activate all basic reserve bids
            current_allocations = Dict(name => bid for (bid, cost, name) in sorted_bids)
            remaining_to_balance = abs_req - sum(values(current_allocations))
            
            # Step 2: Distribute "extra" imbalance while respecting capacity constraints
            # For Down: cannot reduce more than DA_bids[2] (scheduled production)
            # For Up: no P_max provided, so we assume high technical limit
            active_generators = [name for (bid, cost, name) in sorted_bids]
            
            while remaining_to_balance > tol && !isempty(active_generators)
                share = remaining_to_balance / length(active_generators)
                to_remove = String[]
                
                for name in active_generators
                    f = name2idx[name]
                    # Capacity constraint: Down regulation is limited by DA contract
                    limit = (regulation_type == "down") ? DA_bids[2][f] : Inf 
                    
                    current = current_allocations[name]
                    can_add = limit - current
                    
                    added = min(share, max(0.0, can_add))
                    current_allocations[name] += added
                    remaining_to_balance -= added
                    
                    # If generator hits technical limit, remove from further distribution
                    if current_allocations[name] >= limit - tol
                        push!(to_remove, name)
                    end
                end
                filter!(x -> !(x in to_remove), active_generators)
                
                if isempty(active_generators) && remaining_to_balance > tol
                    @warn "Scenario $w: System remains unbalanced by $remaining_to_balance MW - all limits reached."
                    break
                end
            end
            
            # Map results back to output matrices
            for (bid_cap, cost, name) in sorted_bids
                f = name2idx[name]
                final_amt = current_allocations[name]
                if regulation_type == "down"
                    activated_bids_down[w, f] = final_amt
                else
                    activated_bids_up[w, f] = final_amt
                end
                push!(activated_plotting_bids, (final_amt, cost, name))
            end
        end

        # 5. Determine balancing price based on marginal activated bid
        if regulation_type != "none" && !isempty(activated_plotting_bids)
            prices[w] = activated_plotting_bids[end][2] 
        else
            prices[w] = 0.0
        end
    end

    return prices, ppB, ppB_up_contr, ppB_down_contr, activated_bids_up, activated_bids_down
end


function merit_order_price_reserve(rf_up, rf_down, Cf_Rup, Cf_Rdown; dir = "up regulation", 
    figsize = (1000, 800), fontsize = 19)
    Flexibles = vcat(
        [("SCGT", i) for i in 1:4],
        [("CCGT", i) for i in 1:4]
    )    
    bids_up   = [(rf_up[i],   Cf_Rup[i],   Flexibles[i]) for i in 1:nF]
    bids_down = [(rf_down[i], Cf_Rdown[i], Flexibles[i]) for i in 1:nF]

    # --- Sort by cost (ascending merit order) ---
    sorted_bids_up = sort(bids_up, by = x -> x[2])
    sorted_bids_down = sort(bids_down, by = x -> x[2])

    # --- Compute total demand ---
    reserve_up = sum(rf_up) # This shold actually be equal to RI + NI up
    reserve_down = sum(rf_down) # This shold actually be equal to RI + NI up
    
    # --- Find and mark market clearing price UP regulation ---
    market_price_up = 0.0
    cumulative_up = 0.0
    for (bid, cost, _) in sorted_bids_up
        cumulative_up += bid
        if cumulative_up >= reserve_up
            market_price_up = cost
            break
        end
    end

    # --- Find and mark market clearing price DOWN regulation---
    market_price_down = 0.0
    cumulative_down = 0.0
    for (bid, cost, _) in sorted_bids_down
        cumulative_down += bid
        if cumulative_down >= reserve_down
            market_price_down = cost
            break
        end
    end
    
    if dir == "up regulation"
        sorted_bids = sorted_bids_up
        total_demand = reserve_up
        market_price = market_price_up
        price = L"\lambda_r^{\uparrow *}"
    else
        sorted_bids = sorted_bids_down
        total_demand = reserve_down
        market_price = market_price_down
        price = L"\lambda_r^{\downarrow *}"
    end

    # --- Initialize plot ---
    fig = plot(size = figsize,
        xlabel = "Aggregated $(dir) reserve capacity (MW)",
        ylabel = "Cost (€/MW)",
        legend = :bottomright,
        #title = "Merit Order Curve for Reserve Market",
        title = "",
        lw = 2, 
        ylim = (0, market_price*1.5), 
        xlim = (-10, total_demand*1.1),

        guidefont=fontsize,
        tickfont=fontsize,
        xtickfont=fontsize,
        ytickfont=fontsize,
        titlefont=fontsize,
        legendfont=fontsize-2,
        ylabelfont=fontsize,
        xlabelfont=fontsize,
        bottom_margin=10Plots.mm,
        left_margin=10Plots.mm,
        top_margin=5Plots.mm,
        right_margin=10Plots.mm,
        framestyle = :box,
        grid = true,
        gridalpha = 0.1
    )
    # --- Draw each generator block ---
    cumulative = 0.0
    last_cost = nothing
   
    counts_by_cost = Dict{Float64, Int}()

    for (bid, cost, (tech, idx)) in sorted_bids
        color = tech == "SCGT" ? "#3f8fc4" : "#08306b"   # blue / red

        if !isnothing(last_cost) && cost == last_cost && cumulative > 0 && bid > 0
            plot!([cumulative, cumulative], [cost - 0.15, cost + 0.15],
                color =  color, lw = 2, label = "", linealpha=1) # A subtle, slightly thinner line
        end

        # vertical line (cost step) if not first generator
        if !isnothing(last_cost) && bid >0
            plot!([cumulative, cumulative], [last_cost, cost],
                  color = color, lw = 1, label = "")
        end

        # Colors
        if bid > 0
            plot!([cumulative, cumulative + bid], [cost, cost],
                lw = 4, color = color, label = "")
        end

        # add label at the midpoint of the block
        counts_by_cost[cost] = get(counts_by_cost, cost, 0) + 1

        if trunc(Int, bid) > eps()
            annotate!(cumulative + bid/2, cost + cost*0.08, text(string(idx), fontsize, color = color))
        end
        cumulative += bid
        last_cost = cost
    end

    plot!([NaN], [NaN], lw = 4, color = "#3f8fc4", label = "SCGT")
    plot!([NaN], [NaN], lw = 4, color = "#08306b", label = "CCGT")

    vline!([total_demand], label = "Total reserve capacity", ls = :dash, color = "#2E6F40")
    annotate!(cumulative, market_price*1.4, text("$(trunc(Int, cumulative)) MW", fontsize, color = "#388E3C"))
    annotate!(cumulative*1.06, market_price + market_price*0.08, text(price, fontsize, color = :black))

    hline!([market_price], label = "Reserve price", ls = :dash, color = :black)

    return market_price_up, market_price_down, fig
end
####################################################################################

                        # K-FOLD HANDLING

function expected_revenues(a)
    nPolluters, nFolds = size(a)
    nMarkets = length(a[1,1])
    avg_mat = zeros(Float64, nPolluters, nMarkets)
    for i in 1:nPolluters           # each polluter
        for m in 1:nMarkets         # each market
            values = [a[i, f][m] for f in 1:nFolds]
            avg_mat[i, m] = mean(values)
        end
    end

    return avg_mat
end

function sum_over_polluters(mat)
    participants, nFolds = size(mat)
    result = zeros(nFolds, size(mat[1])[1])  # folds × the nb of markets, which will differ for reg
    for k in 1:nFolds
        for i in 1:participants
            result[k, :] += mat[i, k]   # vector addition
        end
    end
    return result
end

function df_in_vs_out(in_mat, out_mat, markets)
    df = DataFrame()
    nFolds, nMarkets = size(in_mat)

    # Add each market with _in and _out columns
    for j in 1:nMarkets
        df[!, "$(markets[j])_in"] = in_mat[:, j]
        df[!, "$(markets[j])_out"] = out_mat[:, j]
    end

    # Compute mean row
    mean_row = Dict()
    for j in 1:nMarkets
        mean_row["$(markets[j])_in"] = mean(in_mat[:, j])
        mean_row["$(markets[j])_out"] = mean(out_mat[:, j])
    end
    
    # Append mean row
    df = vcat(df, DataFrame(mean_row))

    return df, mean_row
end

function initialize_values(inS)
    println("initializing values")
    ################################## Centralized Optimization ##################################
    k0_up   = fill(1.0/nP, length(inS), nP)
    k0_down = fill(1.0/nP, length(inS), nP)

    ppB0_up   = fill(1.0, length(inS), nP)
    ppB0_down = fill(1.0, length(inS), nP)

    ################################## Individual Optimization ##################################
    k0_up_singular   = fill(1.0, length(inS))
    k0_down_singular = fill(1.0, length(inS))  

    ppB0_up_singular   = fill(10.0, length(inS))
    ppB0_down_singular = fill(10.0, length(inS))

    ################################## General ##################################
    z0_up   = fill(1.0, length(inS))
    z0_down = fill(1.0, length(inS))

    Imb_abs_0 = fill(1.0, length(inS))
    s0        = fill(1.0, length(inS))
    s_up0     = fill(1.0, length(inS))
    s_down0   = fill(1.0, length(inS))
    
    return k0_up, k0_down, ppB0_up, ppB0_down, z0_up, z0_down, Imb_abs_0, s0, s_up0, s_down0
end

##profits summed per fold per generator type 
function kfold_profit_table(profits;
        which = (:pol, :flex, :reg),
        label::AbstractString = "Total profit",
        base_suffix::AbstractString = "out",
    )

    # --- helpers ---
    # Extract suffix from e.g. :pol_out_ind -> "out_ind"
    suffix_of(agent::Symbol, field::Symbol) = begin
        pref = string(agent) * "_"
        s = string(field)
        startswith(s, pref) || return nothing
        return s[length(pref)+1:end]
    end

    # Safely get Matrix{Vector}-like field and flatten all units for a fold
    function collect_fold_vec(profits, field::Symbol, fold::Int)
        mat = getproperty(profits, field)  # expected Matrix{Vector}
        v = Float64[]
        for unit in 1:size(mat, 1)
            append!(v, flatten1(mat[unit, fold]))
        end
        return v
    end

    # --- discover relevant fields & suffixes ---
    props = propertynames(profits)

    # Map: suffix => Vector{Symbol fields} (across all chosen agents)
    fields_by_suffix = Dict{String, Vector{Symbol}}()

    for agent in which
        for p in props
            suf = suffix_of(agent, p)
            suf === nothing && continue
            push!(get!(fields_by_suffix, suf, Symbol[]), p)
        end
    end

    isempty(fields_by_suffix) && error("No matching profit fields found for agents: $(which)")

    # Determine k (#folds) from base suffix if available, else from any suffix
    k = if haskey(fields_by_suffix, base_suffix)
        size(getproperty(profits, first(fields_by_suffix[base_suffix])), 2)
    else
        any_suf = first(keys(fields_by_suffix))
        size(getproperty(profits, first(fields_by_suffix[any_suf])), 2)
    end

    # Keep a stable, readable suffix order: base first, then the rest sorted
    suffixes = sort(collect(keys(fields_by_suffix)))
    if base_suffix in suffixes
        suffixes = vcat([base_suffix], filter(!=(base_suffix), suffixes))
    end

    # --- compute per-fold totals/stdevs per suffix ---
    totals = Dict{String, Vector{Float64}}()
    stds   = Dict{String, Vector{Float64}}()

    for suf in suffixes
        totals[suf] = Float64[]
        stds[suf]   = Float64[]
        flds = fields_by_suffix[suf]

        for fold in 1:k
            vec_all = Float64[]
            for fld in flds
                append!(vec_all, collect_fold_vec(profits, fld, fold))
            end
            push!(totals[suf], sum(vec_all))
            push!(stds[suf],   std(vec_all))
        end
    end

    # --- build DataFrame ---
    folds = string.(1:k)
    push!(folds, "Mean")

    which_str = join(string.(which), ", ")

    df = DataFrame(:Fold => folds)

    # Add sum + std columns for each suffix
    for suf in suffixes
        col_sum = Symbol("$label ($which_str) — $suf")
        col_std = Symbol("Std ($which_str) — $suf")

        df[!, col_sum] = vcat(totals[suf], mean(totals[suf]))
        df[!, col_std] = vcat(stds[suf],   mean(stds[suf]))
    end

    # Add differences vs base_suffix (only if base exists)
    if haskey(totals, base_suffix)
        base = totals[base_suffix]
        for suf in suffixes
            suf == base_suffix && continue
            dif = totals[suf] .- base
            df[!, Symbol("Diff ($suf - $base_suffix)")] = vcat(dif, mean(dif))
        end
    end

    return df
end


#profits averaged per fold for each generator 
function kfold_profit_by_unit_table(profits;
        which = (:pol, :flex, :reg),
        label::AbstractString = "Profit",
        base_suffix::AbstractString = "out",
        agent_aliases = Dict(
            :pol  => ["pol", "polluter"],
            :flex => ["flex", "flexible"],
            :reg  => ["reg", "regular"],
        )
    )

    props = propertynames(profits)
    prop_strs = string.(props)

    # suffix_of: wykrywa suffix po dowolnym aliasie agenta
    function suffix_of(agent::Symbol, field::Symbol)
        s = string(field)
        for a in get(agent_aliases, agent, [string(agent)])
            pref = a * "_"
            if startswith(s, pref)
                return s[length(pref)+1:end]  # np "out", "out_ind"
            end
        end
        return nothing
    end

    # Map: (agent, suffix) => Vector{Symbol pól}
    fields_by_agent_suffix = Dict{Tuple{Symbol,String}, Vector{Symbol}}()

    for agent in which
        for p in props
            suf = suffix_of(agent, p)
            suf === nothing && continue
            push!(get!(fields_by_agent_suffix, (agent, suf), Symbol[]), p)
        end
    end

    # lista suffixów
    suffixes = sort(unique(last(k) for k in keys(fields_by_agent_suffix)))
    if base_suffix in suffixes
        suffixes = vcat([base_suffix], filter(!=(base_suffix), suffixes))
    end

    # ustal k (#folds) i #units per agent na bazie pierwszego pola
    function mat_of(field::Symbol)
        getproperty(profits, field)  # expected Matrix{Vector} lub podobne
    end

    # rozmiary per agent mogą być różne, więc robimy per agent osobno
    rows = Vector{NamedTuple}()

    for agent in which
        # znajdź jakieś pole dla tego agenta
        agent_keys = filter(k -> first(k) == agent, collect(keys(fields_by_agent_suffix)))
        isempty(agent_keys) && continue

        # wybierz bazowy suffix jeśli jest, inaczej pierwszy
        suf0 = (agent, base_suffix) in agent_keys ? base_suffix : last(first(agent_keys))
        field0 = first(fields_by_agent_suffix[(agent, suf0)])
        m0 = mat_of(field0)
        n_units = size(m0, 1)
        kfolds  = size(m0, 2)

        for unit in 1:n_units
            row = Dict{Symbol,Any}()
            row[:Agent] = String(agent)
            row[:Unit]  = unit

            for suf in suffixes
                flds = get(fields_by_agent_suffix, (agent, suf), Symbol[])
                isempty(flds) && continue

                # per fold: suma profitów dla tej jednostki, sumując po wszystkich polach danego sufiksu
                per_fold = Float64[]
                for fold in 1:kfolds
                    total = 0.0
                    for fld in flds
                        m = mat_of(fld)
                        total += sum(flatten1(m[unit, fold]))
                    end
                    push!(per_fold, total)
                end

                row[Symbol("$label Mean — $suf")] = mean(per_fold)
                row[Symbol("$label Std — $suf")]  = std(per_fold)
            end

            push!(rows, (; row...))
        end
    end

    df = DataFrame(rows)

    return df
end



function analyze_detailed_profits(profits_data, gen_type, stage, names; models=["", "ind", "bin"])
    df = DataFrame()
    
    for model_suffix in models
        suffix_part = isempty(model_suffix) ? "" : "_$model_suffix"
        key_name = Symbol("$(gen_type)_$(stage)$(suffix_part)")
        
        # Pobranie macierzy (np. 12 generatorów x 10 foldów)
        data_matrix = getfield(profits_data, key_name)
        
        # 1. Sumujemy 4 składowe wewnątrz każdego wektora [Matrix{Float64}]
        # Mapowanie sumy na każdym elemencie macierzy zachowuje jej wymiary (12x10)
        total_profits_matrix = map(sum, data_matrix)
        
        # 2. Średnia i std wzdłuż foldów (dims=2 uśrednia kolumny dla każdego rzędu)
        means = mean(total_profits_matrix, dims=2)[:]
        stds = std(total_profits_matrix, dims=2)[:]
        
        # Sprawdzenie czy liczba nazw zgadza się z rzędami danych
        if length(names) != length(means)
            error("Błąd: Przekazałeś $(length(names)) nazw, a w danych '$key_name' jest $(length(means)) generatorów.")
        end
        
        model_label = isempty(model_suffix) ? "base" : model_suffix
        
        # Budowanie kolumn dla konkretnego modelu
        temp_df = DataFrame(
            :Name => names,
            :Type => fill(uppercase(gen_type), length(names)),
            Symbol("mean_$(model_label)") => means,
            Symbol("std_$(model_label)") => stds
        )
        
        if isempty(df)
            df = temp_df
        else
            # Łączymy wyniki modeli obok siebie
            df = leftjoin(df, temp_df, on=[:Name, :Type])
        end
    end
    
    return df
end


function profits_generators_table(profits_data, types_to_include, stage, names_dict; models=["", "ind", "bin"])
    final_df = DataFrame()
    
    for gen_type in types_to_include
        # Pobieramy nazwy specyficzne dla danego typu (np. pol -> deviator_names)
        current_names = get(names_dict, gen_type, ["Unknown_$i" for i in 1:12])
        
        # Używamy stworzonej wcześniej funkcji do analizy pojedynczego typu
        type_df = analyze_detailed_profits(profits_data, gen_type, stage, current_names; models=models)
        
        # Łączymy wyniki pionowo
        final_df = vcat(final_df, type_df)
    end
    
    return final_df
end



flatten1(x) = x isa AbstractVector{<:Real} ? x :
              x isa AbstractVector{<:AbstractVector{<:Real}} ? vcat(x...) :
              error("Unexpected shape: $(typeof(x))")


function plot_side_by_side_stacked(
    avg_bench, avg_model;
    polluters=nothing,
    markets = ["Day-ahead","Balancing","Penalty"],
    # totals across folds (optional):
    bench_minmax::Union{Nothing,AbstractMatrix}=nothing,  # (nP,2): min,max
    bench_median::Union{Nothing,AbstractVector}=nothing,  # (nP,)
    model_minmax::Union{Nothing,AbstractMatrix}=nothing,
    model_median::Union{Nothing,AbstractVector}=nothing,
    fig_size::Tuple{Int,Int}=(900,560))
    nP = size(avg_bench, 1)
    bar_width = 0.32
    gap = 0.04
    x = 1:nP
    xL = x .- bar_width/2 .- gap/2
    xR = x .+ bar_width/2 .+ gap/2

    # Colors
    blue_dark  = RGBA(0.3, 0.5, 0.8, 1.0)
    blue_med   = RGBA(0.5, 0.7, 0.9, 1.0)
    blue_light = RGBA(0.7, 0.85, 1.0, 1.0)
    orange_dark  = RGBA(1.0, 0.5, 0.0, 1.0)
    orange_med   = RGBA(1.0, 0.7, 0.2, 0.8)
    orange_light = RGBA(1.0, 0.85, 0.6, 0.6)
    blue_colors   = [blue_dark, blue_med, blue_light]
    orange_colors = [orange_dark, orange_med, orange_light]

    # Typography / style
    default(fontfamily="Aptos", guidefont=font(14), tickfont=font(12),
            legendfont=font(12), titlefont=font(18))

    # Omit Benchmark "Penalty" if all zeros
    bench_cols = collect(1:size(avg_bench,2))
    if size(avg_bench,2) >= 3 && all(iszero, vec(avg_bench[:,3]))
        bench_cols = [1,2]
    end

    # y-lims from shown data
    shown_bench = hcat([avg_bench[:,c] for c in bench_cols]...)
    ymin = minimum(vcat(vec(shown_bench), vec(avg_model)))
    ymax = maximum(vcat(vec(shown_bench), vec(avg_model)))
    pad  = 0.08*(ymax - ymin + eps())
    ylim_tuple = (ymin - pad, ymax + pad)

    p = plot(fig_size=fig_size)  # create figure first to apply global opts
    p = plot(size=fig_size,
             xlabel="Polluter", ylabel="Profit [€]",
             title="Profit decomposition per polluter: Model vs Benchmark",
             legend=:bottomright, framestyle=:box,
             background_color=:transparent, background_color_outside=:transparent,
             xlim=(0.5, nP+0.5), ylim=ylim_tuple,  left_margin = 10Plots.mm)

    # Benchmark bars (left)
    for (j,c) in enumerate(bench_cols)
        bar!(xL, avg_bench[:,c]; label="Benchmark – $(markets[c])",
             color=blue_colors[j], bar_width=bar_width, bar_position=:stack, lw=0)
    end
    # Model bars (right)
    for c in 1:length(markets)
        bar!(xR, avg_model[:,c]; label="Model – $(markets[c])",
             color=orange_colors[c], bar_width=bar_width, bar_position=:stack, lw=0)
    end

    # Zero line
    hline!([0.0]; lc=:gray30, ls=:dash, lw=1, label=false)

    # Whiskers helpers
    bench_mean_tot = sum(avg_bench, dims=2)[:]
    model_mean_tot = sum(avg_model, dims=2)[:]

    function add_whiskers!(xs, ymean, mm::AbstractMatrix; col=:gray25)
        lower = ymean .- mm[:,1]
        upper = mm[:,2] .- ymean
        scatter!(xs, ymean; yerror=(lower, upper), ms=0, lc=col, lw=2, label=false)
        cap = bar_width*0.28
        for i in eachindex(xs)
            plot!([xs[i]-cap, xs[i]+cap], [mm[i,1], mm[i,1]]; lc=col, lw=2, label=false)
            plot!([xs[i]-cap, xs[i]+cap], [mm[i,2], mm[i,2]]; lc=col, lw=2, label=false)
        end
    end
    function add_median_ticks!(xs, med; col=:black)
        cap = bar_width*0.22
        for i in eachindex(xs)
            plot!([xs[i]-cap, xs[i]+cap], [med[i], med[i]]; lc=col, lw=2, label=false)
        end
    end

    if bench_minmax !== nothing; add_whiskers!(xL, bench_mean_tot, bench_minmax); end
    if bench_median !== nothing; add_median_ticks!(xL, bench_median); end
    if model_minmax !== nothing; add_whiskers!(xR, model_mean_tot, model_minmax); end
    if model_median !== nothing; add_median_ticks!(xR, model_median); end

    xticks!(x, polluters === nothing ? string.(1:nP) : polluters)
    return p
end
            
function profits_side_by_side(avg_bench, avg_model; polluters=nothing, markets = ["Day-ahead" "Balancing" "Penalty"])
    nP = size(avg_bench, 1)
    bar_width = 0.32
    gap = 0.04
    x = 1:nP

    # Colors
    # Blues for Benchmark: darker to lighter
    blue_dark  = RGBA(0.3, 0.5, 0.8, 1.0)
    blue_med   = RGBA(0.5, 0.7, 0.9, 1.0)
    blue_light = RGBA(0.7, 0.85, 1.0, 1.0)

    # Oranges for Model: darker to lighter
    orange_dark   = RGBA(1.0, 0.5, 0.0, 1.0)
    orange_med    = RGBA(1.0, 0.7, 0.2, 0.8)
    orange_light  = RGBA(1.0, 0.85, 0.6, 0.6)

    blue_colors = [blue_dark, blue_med, blue_light]
    orange_colors = [orange_dark, orange_med, orange_light]

    # Default fonts
    default(
        fontfamily = "Aptos",
        guidefont = font(8),
        tickfont = font(8),
        legendfont = font(8),
        titlefont = font(8)
    )

    # --- Benchmark ---
    p = bar(
        x .- bar_width/2 .- gap/2,
        avg_bench[:,1],
        label = "Benchmark – $(markets[1])",
        color = blue_colors[1],
        bar_width = bar_width,
        bar_position = :stack,
        xlabel = "Polluter",
        ylabel = "Profit [€]",
        title = "Profit decomposition per polluter: Model vs Benchmark",
        legend = :bottomright,
        framestyle = :box,
        ylim = (minimum([avg_bench; avg_model])*1.1, maximum([avg_bench; avg_model])*1.1)
    )

    for m in 2:length(markets)
        bar!(
            x .- bar_width/2 .- gap/2,
            avg_bench[:,m],
            label = "Benchmark – $(markets[m])",
            color = blue_colors[m],
            bar_width = bar_width,
            bar_position = :stack
        )
    end

    # --- Model ---
    for m in 1:length(markets)
        bar!(
            x .+ bar_width/2 .+ gap/2,
            avg_model[:,m],
            label = "Model – $(markets[m])",
            color = orange_colors[m],
            bar_width = bar_width,
            bar_position = :stack
        )
    end

    # X-axis labels
    xticks!(x, polluters === nothing ? string.(1:nP) : polluters)

    return p
end

function average_bids_over_scenarios(prob, B_bids)
    # B_bids is defined as follows :  B_bids_cen_out = [B_bids_pol_cen_out, B_bids_pol_up_out, B_bids_pol_down_out, act_bids_up_out, act_bids_down_out]
    mean_pol_bids= zeros(nP)
    mean_pol_bids_up= zeros(nP)
    mean_pol_bids_down= zeros(nP)
    
    mean_flex_bids_up= zeros(nF)
    mean_flex_bids_down= zeros(nF)
    
    #println(B_bids[1])
    for p in 1:nP
        mean_pol_bids[p] = sum(prob.*B_bids[1][:, p])
        mean_pol_bids_up[p] = sum(prob.*B_bids[2][:, p])
        mean_pol_bids_down[p] = sum(prob.*B_bids[3][:, p])
    end
    #println(mean_pol_bids[1])
    for f in 1:nF
        mean_flex_bids_down[f] = sum(prob.*B_bids[5][:, f])
        mean_flex_bids_up[f] = sum(prob.*B_bids[4][:, f])
    end
    return mean_pol_bids, mean_pol_bids_up, mean_pol_bids_down, mean_flex_bids_up, mean_flex_bids_down  
end


function compute_contributions(bids::AbstractVector, NI_down, NI_up)
    total_imbalance = sum(bids)
    println(total_imbalance)
    if total_imbalance == 0
        return zeros(length(bids))  # no imbalance
    end

    # Determine which bids contribute in the same direction as total
    if total_imbalance > 0
        contributing = bids .> 0
        NI = NI_down
    else
        contributing = bids .< 0
        NI = NI_up
    end
    
    # Compute contributions
    total_contributors = sum(bids[contributing])
    contributions = zeros(length(bids))

    # Avoid division by zero if all bids are opposite direction
    if total_contributors != 0
        contributions[contributing] .= bids[contributing] ./ total_contributors * NI
    end

    return contributions
end


function merge_results(dfs::Vector{DataFrame}, suffixes::Vector{String}; ids=["Scenario", "Gen"])
    @assert length(dfs) == length(suffixes) "Each dataframe must have a suffix"

    # normalize ids to symbols
    ids = Symbol.(ids)

    # make deep copies so we don't mutate originals
    dfs = [copy(df) for df in dfs]

    # normalize all column names in all dfs to symbols
    for df in dfs
        rename!(df, Dict(col => Symbol(col) for col in names(df) if !(col isa Symbol)))
    end

    # rename overlapping columns in each DF (except ids)
    seen = Set(ids)
    for (i, df) in enumerate(dfs)
        suffix = suffixes[i]
        # ensure names are symbols
        cols = Symbol.(names(df))
        overlap = intersect(setdiff(cols, ids), seen)
        if !isempty(overlap)
            rename!(df, Dict(col => Symbol(string(col, suffix)) for col in overlap))
        end
        union!(seen, setdiff(Symbol.(names(df)), ids))
    end

    # join all together
    df_merged = reduce((a, b) -> innerjoin(a, b; on=ids), dfs)

    # build desired column order
    vars = setdiff(union(vcat([Symbol.(names(df)) for df in dfs]...)), ids)
    grouped = Dict{String, Vector{Symbol}}()
    for v in vars
        base = replace(string(v), r"(_.*)$" => "")
        push!(get!(grouped, base, Symbol[]), v)
    end

    ordered = vcat(ids, reduce(vcat, values(grouped)))
    return df_merged[:, ordered]
end



###### Scenario handling ####################


function preprocess_scenario_files(file_list::Vector{String};
                                   season::String = "Spring",
                                   seed::Int = 1234,
                                file_of_generator::Union{Nothing,Vector{Int}} = nothing)
    """
        preprocess_scenario_files(file_list; season="Spring", seed=1234)

    Loads one or more scenario files. If only one file is supplied, it will be
    replicated internally for all generators in generate_Pw().

    Returns:
        centers::Vector{Matrix}
        Prob::Vector{Matrix}
        hours_unique::Vector
        nS_all::Int
    """

    nF = length(file_list)
    centers = Vector{Matrix{Float64}}(undef, nF)
    probabilities    = Vector{Matrix{Float64}}(undef, nF)

    hours_unique = nothing
    nS_all = nothing

    for f in 1:nF
        df = CSV.read(file_list[f], DataFrame)

        # Filter by season
        Sc = df[df.season .== season, :]

        # Initialize shared metadata on first pass
        if f == 1
            hours_unique = sort(unique(Sc.hour))
            nS_all = maximum(Sc.scenario)
        end

        n_hours = length(hours_unique)

        # Random scenario permutation
        Random.seed!(seed)
        perm = randperm(nS_all)

        C = Array{Float64}(undef, nS_all, n_hours)
        P = Array{Float64}(undef, nS_all, n_hours)

        for (j, h) in enumerate(hours_unique)
            sub = Sc[Sc.hour .== h, :]
            sub = sort(sub, :scenario)
            sub = sub[perm, :]

            C[:, j] = sub.value
            P[:, j] = sub.probability
        end

        centers[f] = C
        probabilities[f]    = P
    end

    n_hours = size(centers[1], 2)

    Pw = Vector{Matrix{Float64}}(undef, n_hours)
    #Prob = Vector{Matrix{Float64}}(undef, n_hours)

    # If only one scenario file, force all generators → file 1
    if file_of_generator == nothing
        if nF == 1
            file_of_generator = fill(1, nP)
        else
            error("Multiple files provided → file_of_generator required.")
        end
    end
    
    nOnshore = 6
    first_offshore = 7

    for t in 1:n_hours
        Pw[t] = zeros(nS_all, nP)
        #Prob[t] = zeros(nS_all, nP)

        for p in 1:nOnshore
            f = file_of_generator[p]    # choose correct scenario file
            col = centers[f][:, t]      # pick that file for this hour
            Pw[t][:, p] .= round.(col .* Pmaxp[p], digits=2)
            #Prob[t][:, p] .= round.(probabilities[f][:, t], digits = 4)
        end
        for p in first_offshore:nP 
            f = file_of_generator[p]    # choose correct scenario file
            col = centers[f][:, t]      # pick that file for this hour
            Pw[t][:, p] .= round.(1.6 .* col .* Pmaxp[p], digits=2)
            #Prob[t][:, p] .= round.(probabilities[f][:, t], digits = 4)
        end
    end
    return Pw, centers, probabilities, nS_all
end

function plot_wind_power_scenarios_and_boxplots(Pw_hour, nS_all)
    # Combine offshore and onshore into 6 effective generators by scaling offshore
    Pw_hour_effective = hcat(
        Pw_hour[:,1],          # Onshore σ=0.05s
        Pw_hour[:,3],          # Onshore σ=2.5
        Pw_hour[:,5],          # Onshore σ=1.125
        Pw_hour[:,7],      # Offshore σ=0.05
        Pw_hour[:,9],      # Offshore σ=2.5
        Pw_hour[:,11]      # Offshore σ=1.125
    )

    sigmas_efficient = [0.05, 0.125, 0.25,
                        0.05, 0.125,  0.25]
    # Colors
    sigma_colors = Dict(
        0.05 => "#EE7518",    # Good
        0.125 => "#C85700",   # Average
        0.25 => "#EA985A",      # Bad
    )

    # Sigma for 12 generators (onshore/offshore)
    sigmas = [0.05, 0.05, 0.125, 0.125, 0.25, 0.25, 
            0.05, 0.05, 0.125, 0.125, 0.25, 0.25, ]

    # Legend labels (shared)
    legend_labels = [
        L"$\gamma=0.05$",
        L"$\gamma=0.125$",
        L"$\gamma=0.25$",
    ]

    # --- Scenario lines plot (6 effective generators) ---
    scenarios = 1:nS_all
    p1 = plot(size=(600,600),
            xlabel=L"Scenario index $\omega$",
            ylabel=L"Real Time Wind Power $P_{g, \omega}^{RT}$ (MW)",
            #title=L"$P_{g, \omega}^{RT}$ accross scenarios",
            grid=true,
            gridalpha=0.1,
            #legend=false,             # disable legend here
            guidefont=font(18),
            bottom_margin=10Plots.mm,
            left_margin=10Plots.mm,
            xtickfont=font(16, rotation = 45),
            top_margin=5Plots.mm,
            ytickfont=font(16),
            rotation = 45,
            titlefont=font(18),
            legend=:topright,
            framestyle = :box,
            legendfont=font(16))
            
    colors_array = collect(values(sigma_colors))

    # Example mapping: first 3 onshore, next 3 offshore
    for g in 1:6
        label = if g <=3 legend_labels[g] else "" end
        plot!(p1, scenarios, Pw_hour_effective[:,g],
            color=sigma_colors[sigmas_efficient[g]],
            lw=2,label=label)
    end

    # Add diagonal arrows for onshore/offshore (optional)
    quiver!(p1, [20], [maximum(Pw_hour_effective[:,1])+6], quiver=([10],[ -4]),
            color=:black, lw=2, arrow=true)
    annotate!(p1, 16, maximum(Pw_hour_effective[:,1])+7, text("Onshore", 14, :black))

    quiver!(p1, [38], [minimum(Pw_hour_effective[:,4])-9], quiver=([-10],[4]),
            color=:black, lw=2, arrow=true)
    annotate!(p1, 44, minimum(Pw_hour_effective[:,4])-10, text("Offshore", 14, :black))

    # --- Boxplot for 12 generators ---
    gen_labels = ["D$i" for i in 1:12]  # simplified
    p2 = plot(size=(600,600),
            xlabel="Deviator",
            #ylabel="Wind Power (MW)",
            #title=L"Distribution of $P_{g, \omega}^{RT}$ per deviator",
            grid=true,
            gridalpha=0.1,
            guidefont=font(16),
            xtickfont=font(15),
            rotation = 45,
            ytickfont=font(15),
            titlefont=font(16),
            xticks=(1:12, gen_labels),
            bottom_margin=10Plots.mm,
            right_margin=10Plots.mm,
            top_margin=5Plots.mm,
            legend=:bottomright,
            framestyle = :box,
            legendfont=font(15))

    for g in 1:12
        label = if g in (1, 3, 5)
            legend_labels[mod1(findfirst(==(g), (1,3,5)), length(legend_labels))] 
        else
            ""
        end    
        boxplot!(p2, Pw_hour[:, g],
        positions=[g],
        color=sigma_colors[sigmas[g]],
        boxwidth=0.5, label=label)
    end

    return plot(p1, p2, layout=(1,2), size=(1200,600))
end

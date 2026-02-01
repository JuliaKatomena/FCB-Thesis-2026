####################################################################################

                        # RESERVE MARKET OPT
function reserve_RI_NI(NI_up, NI_down, Cf_Rup, Cf_Rdown)
    results_RI_NI = []

    ########### Model ############
    model = Model(HiGHS.Optimizer)
    set_silent(model)
    @variable(model, 0 <= rf_up[f =1:nF] <= Rf_up[f])
    @variable(model, 0 <= rf_down[f =1:nF] <= Rf_down[f])

    ############ Constraints ############  
    down_constraint = @constraint(model, sum(rf_down[f] for f in 1:nF) == RI_down + NI_down)  
    up_constraint = @constraint(model, sum(rf_up[f] for f in 1:nF) == RI_up + NI_up)

    ############ Objective ############
    @objective(model, Min, 
        sum(Cf_Rup[f]*rf_up[f] + Cf_Rdown[f]*rf_down[f] for f in 1:nF))

    optimize!(model)
    
    ########## Model solved ############
    if termination_status(model) != MOI.OPTIMAL
        #println("Model did not solve to optimality.")
        return nothing
    end

    push!(results_RI_NI, (
        rf_up = deepcopy(value.(rf_up)),
        rf_down = deepcopy(value.(rf_down)),
        obj = objective_value(model),
        Lambda_RI_NI_up = dual(up_constraint),
        Lambda_RI_NI_down = dual(down_constraint),
        RI_NI_up = RI_up + NI_up,
        RI_NI_down = RI_down + NI_down,
    ))

    return results_RI_NI
end
            
function reserve_RI(Cf_Rup, Cf_Rdown)
    results_RI = []

    ########### Model ############
    model = Model(HiGHS.Optimizer)
    set_silent(model)
    @variable(model, 0 <= rf_up[f =1:nF] <= Rf_up[f])
    @variable(model, 0 <= rf_down[f =1:nF] <= Rf_down[f])

    ############ Constraints ############  
    down_constraint = @constraint(model, sum(rf_down[f] for f in 1:nF) == RI_down)  
    up_constraint = @constraint(model, sum(rf_up[f] for f in 1:nF) == RI_up)

    ############ Objective ############
    @objective(model, Min, 
        sum(Cf_Rup[f]*rf_up[f] + Cf_Rdown[f]*rf_down[f] for f in 1:nF))

    optimize!(model)
    
    ########## Model solved ############
    if termination_status(model) != MOI.OPTIMAL
        #println("Model did not solve to optimality.")
        return nothing
    end

    push!(results_RI, (
        rf_up = deepcopy(value.(rf_up)),
        rf_down = deepcopy(value.(rf_down)),
        obj = objective_value(model),
        Lambda_RI_up = dual(up_constraint),
        Lambda_RI_down = dual(down_constraint),
        RI_up = RI_up,
        RI_down = RI_down 
    ))

    return results_RI
end         


####################################################################################

                                    # TA BIN and BIN RELAXED

################# Polluter's Individual Optimization 
function polluter_bin_LP(Pw, nS, prob, pol, central_sol, Q_up, Q_down,
        NI_up, NI_down,
        k_up_prev, k_down_prev, z_up_prev, z_down_prev,
        ppB_up_prev, ppB_down_prev; max_iter=100,
        tol=0.01, penalty_t=100, M = 400, alpha = 1.0)
   
    # make sure we’re working on private copies
    k_up_prev   = copy(k_up_prev)
    k_down_prev = copy(k_down_prev)
    z_up_prev   = copy(z_up_prev)
    z_down_prev = copy(z_down_prev)
    ppB_up_prev = copy(ppB_up_prev)
    ppB_down_prev = copy(ppB_down_prev)
   
    ############## Penalty update (safer) ##############
    max_penalty = 1e5          # hard cap to avoid blow-up
    mult_factor = 1.05         # multiplicative increase (use >1) OR set add_sigma>0 to use additive
    iteration_history = []
 
    for iter in 1:max_iter
       
        ########### Model ############
        model = Model(HiGHS.Optimizer)
        set_silent(model)
 
        @variable(model, 0 <= ppDA<= Pmaxp[pol])
 
        @variable(model, ppB[1:nS])
        @variable(model, 0 <= ppB_up[1:nS] <= Pmaxp[pol])
        @variable(model, 0 <= ppB_down[1:nS] <= Pmaxp[pol])
 
        @variable(model, 0 <= k_up[1:nS] <= 1)
        @variable(model, 0 <= k_down[1:nS] <= 1)
 
        @variable(model, 0 <= z_up[1:nS] <= NI_up)
        @variable(model, 0 <= z_down[1:nS] <= NI_down)
       
        @variable(model, δ_up[1:nS] >= 0)
        @variable(model, δ_down[1:nS] >= 0)
 
        @variable(model, t[1:nS] >= 0)
   
        @variable(model, Imb[1:nS])
        @variable(model, Imb_abs[1:nS] >= 0)

 
        ############ Constraints ############
     
        ##########Polluter ###############
        @constraint(model, [w in 1:nS],
            ppDA + ppB[w] - Pw[w, pol] == 0)
   
        @constraint(model, [w in 1:nS],
            ppB[w] == ppB_down[w] - ppB_up[w])
 
        @constraint(model, [w in 1:nS],ppB_down[w] <= M*(1 - central_sol.y[w, pol]))
        @constraint(model, [w in 1:nS],ppB_up[w] <= M*central_sol.y[w,pol])
 
        @constraint(model, [w=1:nS], z_up[w] - ppB_up[w] - sum(central_sol.ppB_up[w, :]) == 0)
        @constraint(model, [w=1:nS], z_down[w] - ppB_down[w] - sum(central_sol.ppB_down[w, :]) == 0)
       
        @expression(model, linearized_diff_up[w=1:nS],
            k_up_prev[w] * z_up[w] + z_up_prev[w]*k_up[w] - z_up_prev[w]*k_up_prev[w] - (central_sol.s_up[w]*ppB_up[w]))
 
        @expression(model, linearized_diff_down[w=1:nS],
            k_down_prev[w] * z_down[w] + z_down_prev[w]*k_down[w] - z_down_prev[w]*k_down_prev[w] - (central_sol.s_down[w]*ppB_down[w]))
 
        @constraint(model, [w=1:nS], δ_up[w] >= linearized_diff_up[w])
        @constraint(model, [w=1:nS], δ_up[w] >= -linearized_diff_up[w])
 
        @constraint(model, [w=1:nS], δ_down[w] >= linearized_diff_down[w])
        @constraint(model, [w=1:nS], δ_down[w] >= -linearized_diff_down[w])
 
        @constraint(model, [w=1:nS], k_up[w] + sum(central_sol.k_up[w, :]) == central_sol.s_up[w])
        @constraint(model, [w=1:nS], k_down[w] + sum(central_sol.k_down[w, :]) == central_sol.s_down[w])
 
        t_constraint = @constraint(model, [w=1:nS], δ_up[w] + δ_down[w] <= t[w])

 
        @constraint(model, [w=1:nS], Imb[w] == ppB[w] + sum(central_sol.ppB[w, :]))  
        @constraint(model, [w=1:nS], Imb[w] == z_down[w] - z_up[w])
        @constraint(model, [w=1:nS], Imb_abs[w] >= Imb[w])
        @constraint(model, [w=1:nS], Imb_abs[w] >= -Imb[w])
 
        @constraint(model, [w=1:nS], Imb_abs[w] - Imb[w] <= 2*NI_up*(central_sol.s_up[w]))
        @constraint(model, [w=1:nS], Imb_abs[w] + Imb[w] <= 2*NI_down*(central_sol.s_down[w]))
 
        ############ Objective ############
 
        #Penalty for the bilinear terms
        @expression(model, bilinear_cost, sum(prob[w] * penalty_t *t[w] for w in 1:nS))
 
        #Cost of the pollution - with one s
        @expression(model, penalty_cost, sum(prob[w] *(Q_down*k_down[w] + Q_up*k_up[w])  for w in 1:nS))
 
        @objective(model, Max,
            ppDA * (central_sol.Lambda_DA - Cp[pol])
            +sum(prob[w]*(central_sol.Lambda_B[w] -  Cp[pol])*ppB[w] for w in 1:nS)
            - bilinear_cost - penalty_cost )
 
        optimize!(model)
       
        ########## Model solved ############
        if termination_status(model) != MOI.OPTIMAL
            #println("Model did not solve to optimality.")
            return nothing
        end
        Lambda_t = dual.(t_constraint)
        delta_1 = 0.01
 
        zk_val_up = value.(z_up) .* value.(k_up)
        zk_val_down = value.(z_down) .* value.(k_down)
 
        pct_diff_val_up = 100 .* (zk_val_up .- value.(ppB_up)) ./ (value.(ppB_up).*value.(central_sol.s_up) .+ eps())
        pct_diff_val_down = 100 .* (zk_val_down .- value.(ppB_down)) ./ (value.(ppB_down).*value.(central_sol.s_down) .+ eps())
 
        push!(iteration_history, (
            iter = iter,
            Gen = pol,
            ppB = deepcopy(value.(ppB)),
            z_up = deepcopy(value.(z_up)),
            k_up = deepcopy(value.(k_up)),
            ppB_up = deepcopy(value.(ppB_up)),
            diff_up = deepcopy(pct_diff_val_up),
            z_down = deepcopy(value.(z_down)),
            k_down = deepcopy(value.(k_down)),
            ppB_down = deepcopy(value.(ppB_down)),
            diff_down = deepcopy(pct_diff_val_down),
            ppDA = deepcopy(value(ppDA)),
            Imb = deepcopy(value.(Imb)),
            Imb_abs = deepcopy(value.(Imb_abs)),
            obj = objective_value(model),
            bilinear_cost = value(bilinear_cost),
            penalty_cost = value(penalty_cost),
            obj_sans_imbalance = objective_value(model)
        ))
 
        ######### Convergence check ############
        diff_matrix_up = zk_val_up .- value.(ppB_up).*value.(central_sol.s_up)
        diff_matrix_down = zk_val_down .- value.(ppB_down).*value.(central_sol.s_down)
 
        step_norm = sqrt(
            norm(value.(k_up) - k_up_prev)^2 +
            norm(value.(z_up) - z_up_prev)^2 +
            norm(value.(ppB_up) - ppB_up_prev)^2 +
            norm(value.(k_down) - k_down_prev)^2 +
            norm(value.(z_down) - z_down_prev)^2 +
            norm(value.(ppB_down) - ppB_down_prev)^2
        )
 
        if step_norm < tol && all(abs.(diff_matrix_up) .< tol) && all(abs.(diff_matrix_down) .< tol) && all(value.(t) .< tol)
            if any(abs.(diff_matrix_up) .> tol)
                #println("⚠️ Some values in diff_matrix_up exceed tolerance $tol")
            end
            if any(abs.(diff_matrix_up) .> tol)
                #println("⚠️ Some values in diff_matrix_down exceed tolerance $tol")
            end
            #println("TA Ind Converged after  $iter iterations")
            return iteration_history
        end
 
        # ---- UP ----
        inv_step = 1.0 / max(step_norm, eps())
        rk = min(inv_step, norm(Lambda_t, 1) + delta_1)
 
        # replace your if-condition with this (example for UP):
        if (penalty_t < rk)
            penalty_t = min(max_penalty, penalty_t * mult_factor)
        end
 
        if iter % 100 == 0  # print every 100 iterations
            #println("Iteration $iter summary:")
            #println("  UP diag: step_norm=$(step_norm_up), inv_step=$(inv_step_up), norm(Lambda_up,1)=$(norm(Lambda_up,1)), penalty_up=$(penalty_up), violation_up=$(violation_up)")
            #println("  DOWN diag: step_norm=$(step_norm_down), inv_step=$(inv_step_down), norm(Lambda_down,1)=$(norm(Lambda_down,1)), penalty_down=$(penalty_down), violation_down=$(violation_down)")
        end
 
        ############# Update k, z, ppB and y ############
        k_up_prev    .= (1-alpha).*k_up_prev    .+ alpha.*value.(k_up)
        z_up_prev    .= (1-alpha).*z_up_prev    .+ alpha.*value.(z_up)
        ppB_up_prev  .= (1-alpha).*ppB_up_prev  .+ alpha.*value.(ppB_up)
 
        k_down_prev    .= (1-alpha).*k_down_prev    .+ alpha.*value.(k_down)
        z_down_prev    .= (1-alpha).*z_down_prev    .+ alpha.*value.(z_down)
        ppB_down_prev  .= (1-alpha).*ppB_down_prev  .+ alpha.*value.(ppB_down)
    end
 
    #println("TA Did not converge after $max_iter iterations")
    return iteration_history
end

function polluter_bin_LP_k(Pw, nS, prob, pol, central_sol, Q_up, Q_down,
        NI_up, NI_down,
        k_up_prev, k_down_prev, z_up_prev, z_down_prev,
        ppB_up_prev, ppB_down_prev; max_iter=100,
        tol=0.01, penalty_t=100, M = 400, alpha = 1.0)
   
    # make sure we’re working on private copies
    k_up_prev   = copy(k_up_prev)
    k_down_prev = copy(k_down_prev)
    z_up_prev   = copy(z_up_prev)
    z_down_prev = copy(z_down_prev)
    ppB_up_prev = copy(ppB_up_prev)
    ppB_down_prev = copy(ppB_down_prev)
   
    ############## Penalty update (safer) ##############
    max_penalty = 1e5          # hard cap to avoid blow-up
    mult_factor = 1.05         # multiplicative increase (use >1) OR set add_sigma>0 to use additive
    iteration_history = []
 
    for iter in 1:max_iter
       
        ########### Model ############
        model = Model(HiGHS.Optimizer)
        set_silent(model)
 
        @variable(model, 0 <= ppDA<= Pmaxp[pol])
 
        @variable(model, ppB[1:nS])
        @variable(model, 0 <= ppB_up[1:nS] <= Pmaxp[pol])
        @variable(model, 0 <= ppB_down[1:nS] <= Pmaxp[pol])
 
        @variable(model, 0 <= k_up[1:nS] <= 1)
        @variable(model, 0 <= k_down[1:nS] <= 1)
 
        @variable(model, 0 <= z_up[1:nS] <= NI_up)
        @variable(model, 0 <= z_down[1:nS] <= NI_down)
       
        @variable(model, δ_up[1:nS] >= 0)
        @variable(model, δ_down[1:nS] >= 0)
 
        @variable(model, t[1:nS] >= 0)
   
        @variable(model, Imb[1:nS])
        @variable(model, Imb_abs[1:nS] >= 0)

       
        @variable(model, k_max[1:nS] >= 0)
 
        ############ Constraints ############
     
        ##########Polluter ###############
        @constraint(model, [w in 1:nS],
            ppDA + ppB[w] - Pw[w, pol] == 0)
   
        @constraint(model, [w in 1:nS],
            ppB[w] == ppB_down[w] - ppB_up[w])
 
        @constraint(model, [w in 1:nS],ppB_down[w] <= M*(1 - central_sol.y[w, pol]))
        @constraint(model, [w in 1:nS],ppB_up[w] <= M*central_sol.y[w,pol])
 
        @constraint(model, [w=1:nS], z_up[w] - ppB_up[w] - sum(central_sol.ppB_up[w, :]) == 0)
        @constraint(model, [w=1:nS], z_down[w] - ppB_down[w] - sum(central_sol.ppB_down[w, :]) == 0)
       
        @expression(model, linearized_diff_up[w=1:nS],
            k_up_prev[w] * z_up[w] + z_up_prev[w]*k_up[w] - z_up_prev[w]*k_up_prev[w] - (central_sol.s_up[w]*ppB_up[w]))
 
        @expression(model, linearized_diff_down[w=1:nS],
            k_down_prev[w] * z_down[w] + z_down_prev[w]*k_down[w] - z_down_prev[w]*k_down_prev[w] - (central_sol.s_down[w]*ppB_down[w]))
 
        @constraint(model, [w=1:nS], δ_up[w] >= linearized_diff_up[w])
        @constraint(model, [w=1:nS], δ_up[w] >= -linearized_diff_up[w])
 
        @constraint(model, [w=1:nS], δ_down[w] >= linearized_diff_down[w])
        @constraint(model, [w=1:nS], δ_down[w] >= -linearized_diff_down[w])
 
        @constraint(model, [w=1:nS], k_up[w] + sum(central_sol.k_up[w, :]) == central_sol.s_up[w])
        @constraint(model, [w=1:nS], k_down[w] + sum(central_sol.k_down[w, :]) == central_sol.s_down[w])
 
        t_constraint = @constraint(model, [w=1:nS], δ_up[w] + δ_down[w] <= t[w])

        @constraint(model, [w=1:nS], k_max[w] >= k_up[w])
        @constraint(model, [w=1:nS], k_max[w] >= k_down[w])
 
        @constraint(model, [w=1:nS], Imb[w] == ppB[w] + sum(central_sol.ppB[w, :]))  
        @constraint(model, [w=1:nS], Imb[w] == z_down[w] - z_up[w])
        @constraint(model, [w=1:nS], Imb_abs[w] >= Imb[w])
        @constraint(model, [w=1:nS], Imb_abs[w] >= -Imb[w])
 
        @constraint(model, [w=1:nS], Imb_abs[w] - Imb[w] <= 2*NI_up*(central_sol.s_up[w]))
        @constraint(model, [w=1:nS], Imb_abs[w] + Imb[w] <= 2*NI_down*(central_sol.s_down[w]))
 
        ############ Objective ############
 
        #Penalty for the bilinear terms
        @expression(model, bilinear_cost, sum(prob[w] * penalty_t *t[w] for w in 1:nS))
 
        #Cost of the pollution - with one s
        @expression(model, penalty_cost, sum(prob[w] *(Q_down*k_down[w] + Q_up*k_up[w])  for w in 1:nS))

        #Equalization penalty 
        #@expression(model, spread_penalty, sum(prob[w]*(dev_up[w] + dev_down[w]) for w in 1:nS))
        @expression(model, spread_penalty, 200*sum(k_max[w] for w in 1:nS))
 
        @objective(model, Max,
            ppDA * (central_sol.Lambda_DA - Cp[pol])
            +sum(prob[w]*(central_sol.Lambda_B[w] -  Cp[pol])*ppB[w] for w in 1:nS)
            - bilinear_cost - penalty_cost - spread_penalty)
 
        optimize!(model)
       
        ########## Model solved ############
        if termination_status(model) != MOI.OPTIMAL
            #println("Model did not solve to optimality.")
            return nothing
        end
        Lambda_t = dual.(t_constraint)
        delta_1 = 0.01
 
        zk_val_up = value.(z_up) .* value.(k_up)
        zk_val_down = value.(z_down) .* value.(k_down)
 
        pct_diff_val_up = 100 .* (zk_val_up .- value.(ppB_up)) ./ (value.(ppB_up).*value.(central_sol.s_up) .+ eps())
        pct_diff_val_down = 100 .* (zk_val_down .- value.(ppB_down)) ./ (value.(ppB_down).*value.(central_sol.s_down) .+ eps())
 
        push!(iteration_history, (
            iter = iter,
            Gen = pol,
            ppB = deepcopy(value.(ppB)),
            z_up = deepcopy(value.(z_up)),
            k_up = deepcopy(value.(k_up)),
            ppB_up = deepcopy(value.(ppB_up)),
            diff_up = deepcopy(pct_diff_val_up),
            z_down = deepcopy(value.(z_down)),
            k_down = deepcopy(value.(k_down)),
            ppB_down = deepcopy(value.(ppB_down)),
            diff_down = deepcopy(pct_diff_val_down),
            ppDA = deepcopy(value(ppDA)),
            Imb = deepcopy(value.(Imb)),
            Imb_abs = deepcopy(value.(Imb_abs)),
            obj = objective_value(model),
            bilinear_cost = value(bilinear_cost),
            penalty_cost = value(penalty_cost),
            spread_penalty = value(spread_penalty),
            obj_sans_imbalance = objective_value(model) - value(spread_penalty)
        ))
 
        ######### Convergence check ############
        diff_matrix_up = zk_val_up .- value.(ppB_up).*value.(central_sol.s_up)
        diff_matrix_down = zk_val_down .- value.(ppB_down).*value.(central_sol.s_down)
 
        step_norm = sqrt(
            norm(value.(k_up) - k_up_prev)^2 +
            norm(value.(z_up) - z_up_prev)^2 +
            norm(value.(ppB_up) - ppB_up_prev)^2 +
            norm(value.(k_down) - k_down_prev)^2 +
            norm(value.(z_down) - z_down_prev)^2 +
            norm(value.(ppB_down) - ppB_down_prev)^2
        )
 
        if step_norm < tol && all(abs.(diff_matrix_up) .< tol) && all(abs.(diff_matrix_down) .< tol) && all(value.(t) .< tol)
            if any(abs.(diff_matrix_up) .> tol)
                #println("⚠️ Some values in diff_matrix_up exceed tolerance $tol")
            end
            if any(abs.(diff_matrix_up) .> tol)
                #println("⚠️ Some values in diff_matrix_down exceed tolerance $tol")
            end
            #println("TA Ind Converged after  $iter iterations")
            return iteration_history
        end
 
        # ---- UP ----
        inv_step = 1.0 / max(step_norm, eps())
        rk = min(inv_step, norm(Lambda_t, 1) + delta_1)
 
        # replace your if-condition with this (example for UP):
        if (penalty_t < rk)
            penalty_t = min(max_penalty, penalty_t * mult_factor)
        end
 
        if iter % 100 == 0  # print every 100 iterations
            #println("Iteration $iter summary:")
            #println("  UP diag: step_norm=$(step_norm_up), inv_step=$(inv_step_up), norm(Lambda_up,1)=$(norm(Lambda_up,1)), penalty_up=$(penalty_up), violation_up=$(violation_up)")
            #println("  DOWN diag: step_norm=$(step_norm_down), inv_step=$(inv_step_down), norm(Lambda_down,1)=$(norm(Lambda_down,1)), penalty_down=$(penalty_down), violation_down=$(violation_down)")
        end
 
        ############# Update k, z, ppB and y ############
        k_up_prev    .= (1-alpha).*k_up_prev    .+ alpha.*value.(k_up)
        z_up_prev    .= (1-alpha).*z_up_prev    .+ alpha.*value.(z_up)
        ppB_up_prev  .= (1-alpha).*ppB_up_prev  .+ alpha.*value.(ppB_up)
 
        k_down_prev    .= (1-alpha).*k_down_prev    .+ alpha.*value.(k_down)
        z_down_prev    .= (1-alpha).*z_down_prev    .+ alpha.*value.(z_down)
        ppB_down_prev  .= (1-alpha).*ppB_down_prev  .+ alpha.*value.(ppB_down)
    end
 
    #println("TA Did not converge after $max_iter iterations")
    return iteration_history
end
 
 

function polluter_bin_LP_no_step(Pw, nS, prob, pol, central_sol, Q_up, Q_down,
        NI_up, NI_down,
        k_up_prev, k_down_prev, z_up_prev, z_down_prev,
        ppB_up_prev, ppB_down_prev; max_iter=100,
        tol=0.01, penalty_t=100, M = 400, alpha = 1.0)
   
    # make sure we’re working on private copies
    k_up_prev   = copy(k_up_prev)
    k_down_prev = copy(k_down_prev)
    z_up_prev   = copy(z_up_prev)
    z_down_prev = copy(z_down_prev)
    ppB_up_prev = copy(ppB_up_prev)
    ppB_down_prev = copy(ppB_down_prev)
   
    ############## Penalty update (safer) ##############
    max_penalty = 1e5          # hard cap to avoid blow-up
    mult_factor = 1.05         # multiplicative increase (use >1) OR set add_sigma>0 to use additive
    iteration_history = []
 
    for iter in 1:max_iter
       
        ########### Model ############
        model = Model(HiGHS.Optimizer)
        set_silent(model)
 
        @variable(model, 0 <= ppDA<= Pmaxp[pol])
 
        @variable(model, ppB[1:nS])
        @variable(model, 0 <= ppB_up[1:nS] <= Pmaxp[pol])
        @variable(model, 0 <= ppB_down[1:nS] <= Pmaxp[pol])
 
        @variable(model, 0 <= k_up[1:nS] <= 1)
        @variable(model, 0 <= k_down[1:nS] <= 1)
 
        @variable(model, 0 <= z_up[1:nS] <= NI_up)
        @variable(model, 0 <= z_down[1:nS] <= NI_down)
       
        @variable(model, δ_up[1:nS] >= 0)
        @variable(model, δ_down[1:nS] >= 0)
 
        @variable(model, t[1:nS] >= 0)
   
        @variable(model, Imb[1:nS])
        @variable(model, Imb_abs[1:nS] >= 0)
 
        ############ Constraints ############
     
        ##########Polluter ###############
        @constraint(model, [w in 1:nS],
            ppDA + ppB[w] - Pw[w, pol] == 0)
   
        @constraint(model, [w in 1:nS],
            ppB[w] == ppB_down[w] - ppB_up[w])
 
        @constraint(model, [w in 1:nS],ppB_down[w] <= M*(1 - central_sol.y[w, pol]))
        @constraint(model, [w in 1:nS],ppB_up[w] <= M*central_sol.y[w,pol])
 
        @constraint(model, [w=1:nS], z_up[w] - ppB_up[w] - sum(central_sol.ppB_up[w, :]) == 0)
        @constraint(model, [w=1:nS], z_down[w] - ppB_down[w] - sum(central_sol.ppB_down[w, :]) == 0)
       
        @expression(model, linearized_diff_up[w=1:nS],
            k_up_prev[w] * z_up[w] + z_up_prev[w]*k_up[w] - z_up_prev[w]*k_up_prev[w] - (central_sol.s_up[w]*ppB_up[w]))
 
        @expression(model, linearized_diff_down[w=1:nS],
            k_down_prev[w] * z_down[w] + z_down_prev[w]*k_down[w] - z_down_prev[w]*k_down_prev[w] - (central_sol.s_down[w]*ppB_down[w]))
 
        @constraint(model, [w=1:nS], δ_up[w] >= linearized_diff_up[w])
        @constraint(model, [w=1:nS], δ_up[w] >= -linearized_diff_up[w])
 
        @constraint(model, [w=1:nS], δ_down[w] >= linearized_diff_down[w])
        @constraint(model, [w=1:nS], δ_down[w] >= -linearized_diff_down[w])
 
        @constraint(model, [w=1:nS], k_up[w] + sum(central_sol.k_up[w, :]) == central_sol.s_up[w])
        @constraint(model, [w=1:nS], k_down[w] + sum(central_sol.k_down[w, :]) == central_sol.s_down[w])
 
        t_constraint = @constraint(model, [w=1:nS], δ_up[w] + δ_down[w] <= t[w])
 
        @constraint(model, [w=1:nS], Imb[w] == ppB[w] + sum(central_sol.ppB[w, :]))  
        @constraint(model, [w=1:nS], Imb[w] == z_down[w] - z_up[w])
        @constraint(model, [w=1:nS], Imb_abs[w] >= Imb[w])
        @constraint(model, [w=1:nS], Imb_abs[w] >= -Imb[w])
 
        @constraint(model, [w=1:nS], Imb_abs[w] - Imb[w] <= 2*NI_up*(central_sol.s_up[w]))
        @constraint(model, [w=1:nS], Imb_abs[w] + Imb[w] <= 2*NI_down*(central_sol.s_down[w]))
 
        ############ Objective ############
 
        #Penalty for the bilinear terms
        @expression(model, bilinear_cost, sum(prob[w] * penalty_t *t[w] for w in 1:nS))
 
        #Cost of the pollution - with one s
        @expression(model, penalty_cost, sum(prob[w] *(Q_down*k_down[w] + Q_up*k_up[w])  for w in 1:nS))
 
        @objective(model, Max,
            ppDA * (central_sol.Lambda_DA - Cp[pol])
            +sum(prob[w]*(central_sol.Lambda_B[w] -  Cp[pol])*ppB[w] for w in 1:nS)
            - bilinear_cost - penalty_cost)
 
        optimize!(model)
       
        ########## Model solved ############
        if termination_status(model) != MOI.OPTIMAL
            #println("Model did not solve to optimality.")
            return nothing
        end
        Lambda_t = dual.(t_constraint)
        delta_1 = 0.01
 
        zk_val_up = value.(z_up) .* value.(k_up)
        zk_val_down = value.(z_down) .* value.(k_down)
 
        pct_diff_val_up = 100 .* (zk_val_up .- value.(ppB_up)) ./ (value.(ppB_up).*value.(central_sol.s_up) .+ eps())
        pct_diff_val_down = 100 .* (zk_val_down .- value.(ppB_down)) ./ (value.(ppB_down).*value.(central_sol.s_down) .+ eps())
 
        push!(iteration_history, (
            iter = iter,
            Gen = pol,
            ppB = deepcopy(value.(ppB)),
            z_up = deepcopy(value.(z_up)),
            k_up = deepcopy(value.(k_up)),
            ppB_up = deepcopy(value.(ppB_up)),
            diff_up = deepcopy(pct_diff_val_up),
            z_down = deepcopy(value.(z_down)),
            k_down = deepcopy(value.(k_down)),
            ppB_down = deepcopy(value.(ppB_down)),
            diff_down = deepcopy(pct_diff_val_down),
            ppDA = deepcopy(value(ppDA)),
            Imb = deepcopy(value.(Imb)),
            Imb_abs = deepcopy(value.(Imb_abs)),
            obj = objective_value(model),
            bilinear_cost = value(bilinear_cost),
            penalty_cost = value(penalty_cost),
            obj_sans_imbalance = objective_value(model)
        ))
 
        ######### Convergence check ############
        diff_matrix_up = zk_val_up .- value.(ppB_up).*value.(central_sol.s_up)
        diff_matrix_down = zk_val_down .- value.(ppB_down).*value.(central_sol.s_down)
 
        step_norm = sqrt(
            norm(value.(k_up) - k_up_prev)^2 +
            norm(value.(z_up) - z_up_prev)^2 +
            norm(value.(ppB_up) - ppB_up_prev)^2 +
            norm(value.(k_down) - k_down_prev)^2 +
            norm(value.(z_down) - z_down_prev)^2 +
            norm(value.(ppB_down) - ppB_down_prev)^2
        )
        ##step_norm < tol &&
        if  all(abs.(diff_matrix_up) .< tol) && all(abs.(diff_matrix_down) .< tol) && all(value.(t) .< tol)
            if any(abs.(diff_matrix_up) .> tol)
                #println("⚠️ Some values in diff_matrix_up exceed tolerance $tol")
            end
            if any(abs.(diff_matrix_up) .> tol)
                #println("⚠️ Some values in diff_matrix_down exceed tolerance $tol")
            end
            #println("TA Ind Converged after  $iter iterations")
            return iteration_history
        end
 
        # ---- UP ----
        inv_step = 1.0 / max(step_norm, eps())
        rk = min(inv_step, norm(Lambda_t, 1) + delta_1)
 
        # replace your if-condition with this (example for UP):
        if (penalty_t < rk)
            penalty_t = min(max_penalty, penalty_t * mult_factor)
        end
 

 
        ############# Update k, z, ppB and y ############
        k_up_prev    .= (1-alpha).*k_up_prev    .+ alpha.*value.(k_up)
        z_up_prev    .= (1-alpha).*z_up_prev    .+ alpha.*value.(z_up)
        ppB_up_prev  .= (1-alpha).*ppB_up_prev  .+ alpha.*value.(ppB_up)
 
        k_down_prev    .= (1-alpha).*k_down_prev    .+ alpha.*value.(k_down)
        z_down_prev    .= (1-alpha).*z_down_prev    .+ alpha.*value.(z_down)
        ppB_down_prev  .= (1-alpha).*ppB_down_prev  .+ alpha.*value.(ppB_down)
    end
 
    #println("TA Did not converge after $max_iter iterations")
    return iteration_history
end


# We need to distinguish the binary opt for flex as it includes the binary variables in some of their constraints
function flex_bin_LP(Pw, nS, prob, flex, central_sol, rf_up, rf_down)
 
    ########### Model ############
    model = Model(HiGHS.Optimizer)
    set_silent(model)
    @variable(model, 0 <= pfDA)
 
    @variable(model, 0 <= pfB_up[1:nS])
    @variable(model, 0 <= pfB_down[1:nS])
 
    ############ Constraints ############        
 
    ########## Flex ###############
 
    @constraint(model, rf_down[flex] <= pfDA)
    @constraint(model, pfDA <= Pmaxf[flex] - rf_up[flex])
    @constraint(model, [w=1:nS], pfB_up[w] <= rf_up[flex]*central_sol.s_up[w])
    @constraint(model, [w=1:nS], pfB_down[w] <= rf_down[flex]*central_sol.s_down[w])
 
    ############ Objective ############
 
    @objective(model, Max,
        pfDA * (central_sol.Lambda_DA - Cf[flex])
        + sum(prob[w]*(central_sol.Lambda_B[w] -  Cf[flex])*(pfB_up[w]-pfB_down[w]) for w in 1:nS ))
 
    optimize!(model)
   
    ########## Model solved ############
    if termination_status(model) != MOI.OPTIMAL
        #println("Model did not solve to optimality.")
        return nothing
    end
 
    return [(
        Gen = flex,
        pfDA = JuMP.value(pfDA),
        pfB_up = value.(pfB_up),
        pfB_down = value.(pfB_down),
        obj = objective_value(model)
    )]
end

################# Centralized Optimization with several flexible gen
function TA_bin_GDCA(Pw, nS, prob, Q_up, Q_down,  rf_up, rf_down,
    NI_up, NI_down,
    k_up_prev, k_down_prev, z_up_prev, 
    z_down_prev, ppB_up_prev, ppB_down_prev,
    s_up_prev, s_down_prev;
    max_iter=2000, tol=0.01, penalty_t=500,
    M=300, alpha = 0.9)

    # make sure we’re working on private copies
    k_up_prev   = copy(k_up_prev)
    k_down_prev = copy(k_down_prev)
    z_up_prev   = copy(z_up_prev)
    z_down_prev = copy(z_down_prev)
    ppB_up_prev = copy(ppB_up_prev)
    ppB_down_prev = copy(ppB_down_prev)
    s_up_prev   = copy(s_up_prev)
    s_down_prev = copy(s_down_prev)
    
    iteration_history = []

    ############## Penalty update (safer) ##############
    # parameters you can tune:
    max_penalty = 1e5          # hard cap to avoid blow-up
    mult_factor = 1.5         # multiplicative increase (use >1) OR set add_sigma>0 to use additive
    iteration_history = []

    for iter in 1:max_iter
        ########### Model ############
        model = Model(HiGHS.Optimizer)
        set_silent(model)


        @variable(model, 0 <= ppDA[p=1:nP] <= Pmaxp[p])
        @variable(model, 0 <= pfDA[f=1:nF] <= Pmaxf[f])
        @variable(model, 0 <= pgDA[r=1:nR] <= Pmaxg[r])

        @variable(model, ppB[1:nS, p=1:nP])
        @variable(model, 0 <= ppB_up[1:nS, p=1:nP] <= Pmaxp[p])
        @variable(model, 0 <= ppB_down[1:nS, p=1:nP] <= Pmaxp[p])

        @variable(model, 0 <= pfB_up[1:nS, f=1:nF] <= rf_up[f])
        @variable(model, 0 <= pfB_down[1:nS, f=1:nF] <= rf_down[f])

        @variable(model, 0 <= k_up[1:nS, p=1:nP] <= 1)
        @variable(model, 0 <= k_down[1:nS, p=1:nP] <= 1)

        @variable(model, 0 <= z_up[1:nS] <= NI_up)
        @variable(model, 0 <= z_down[1:nS] <= NI_down)
        
        @variable(model, δ_up[1:nS, 1:nP] >= 0)
        @variable(model, δ_down[1:nS, 1:nP] >= 0)
        @variable(model, t[1:nS, 1:nP] >= 0)
    
        @variable(model, Imb[1:nS])
        @variable(model, Imb_abs[1:nS] >= 0)
        @variable(model, s_up[1:nS], Bin)
        @variable(model, s_down[1:nS], Bin)

        @variable(model, y[1:nS, 1:nP], Bin)

        ############ Constraints ############
     
        ##########Polluter ###############
        @constraint(model, [w in 1:nS, p in 1:nP],
            ppDA[p] + ppB[w, p] - Pw[w, p] == 0)
        
        @constraint(model, [w in 1:nS, p in 1:nP],
            ppB[w, p] == ppB_down[w, p] - ppB_up[w, p])

        @constraint(model, [w in 1:nS, p in 1:nP], ppB_down[w, p] <= M*(1 - y[w,p]))
        @constraint(model, [w in 1:nS, p in 1:nP], ppB_up[w, p] <= M*y[w,p])
    
        @constraint(model, [w=1:nS], z_up[w] - sum(ppB_up[w, p] for p in 1:nP) == 0)
        @constraint(model, [w=1:nS], z_down[w] - sum(ppB_down[w, p] for p in 1:nP) == 0)
        
        @expression(model, linearized_diff_up[w=1:nS, p=1:nP],
            k_up_prev[w,p] * z_up[w] + z_up_prev[w]*k_up[w,p] - z_up_prev[w]*k_up_prev[w,p] - (s_up_prev[w]*ppB_up[w,p] + s_up[w]*ppB_up_prev[w,p] - s_up_prev[w]*ppB_up_prev[w, p]))

        @expression(model, linearized_diff_down[w=1:nS, p=1:nP],
            k_down_prev[w,p] * z_down[w] + z_down_prev[w]*k_down[w,p] - z_down_prev[w]*k_down_prev[w,p] - (s_down_prev[w]*ppB_down[w,p] + s_down[w]*ppB_down_prev[w,p] - s_down_prev[w]*ppB_down_prev[w, p]))

        @constraint(model, [w=1:nS, p=1:nP], δ_up[w,p] >= linearized_diff_up[w,p])
        @constraint(model, [w=1:nS, p=1:nP], δ_up[w,p] >= -linearized_diff_up[w,p])

        @constraint(model, [w=1:nS, p=1:nP], δ_down[w,p] >= linearized_diff_down[w,p])
        @constraint(model, [w=1:nS, p=1:nP], δ_down[w,p] >= -linearized_diff_down[w,p])

        @constraint(model, [w=1:nS], sum(k_up[w, p] for p in 1:nP) == s_up[w])
        @constraint(model, [w=1:nS], sum(k_down[w, p] for p in 1:nP) == s_down[w])
 
        t_constraint = @constraint(model, [w=1:nS, p=1:nP], δ_down[w,p] + δ_up[w,p] <= t[w,p])

        #@constraint(model, [w=1:nS], Imb[w] == sum(ppB[w, p] for p in 1:nP)) 
        @constraint(model, [w=1:nS], Imb[w] == z_down[w] - z_up[w]) 
        @constraint(model, [w=1:nS], Imb_abs[w] >= Imb[w])
        @constraint(model, [w=1:nS], Imb_abs[w] >= -Imb[w])
        @constraint(model, [w=1:nS], s_up[w] + s_down[w] <= 1)

        @constraint(model, [w=1:nS], Imb_abs[w] - Imb[w] <= 2*NI_up*(s_up[w]))
        @constraint(model, [w=1:nS], Imb_abs[w] + Imb[w] <= 2*NI_down*(s_down[w]))

        ########## Flex ###############

        @constraint(model,[f=1:nF], rf_down[f] <= pfDA[f])
        @constraint(model,[f=1:nF], pfDA[f] <= Pmaxf[f] - rf_up[f])
        
        @constraint(model, [w=1:nS, f=1:nF], pfB_down[w, f] <= rf_down[f]*s_down[w])
        @constraint(model, [w=1:nS, f=1:nF], pfB_up[w, f] <= rf_up[f]*s_up[w])

        ############### DA market clearing ###############

        demand_con = @constraint(model, 
            -D + sum(ppDA[p] for p in 1:nP) + sum(pfDA[f] for f in 1:nF) + sum(pgDA[r] for r in 1:nR) == 0)

        ############### Balancing market clearing ###############

        balancing_cons = @constraint(model, [w in 1:nS], 
            prob[w] * (sum(ppB[w, p] for p in 1:nP) - (sum(pfB_down[w,f] - pfB_up[w,f] for f in 1:nF))) == 0)

        ############ Objective ############

        #Penalty for the bilinear terms
        @expression(model, bilinear_cost, sum(prob[w] * penalty_t *t[w,p] for w in 1:nS, p in 1:nP))

        #Cost of the pollution - with one s 
        @expression(model, penalty_cost, sum(prob[w] *(Q_down*k_down[w,p] + Q_up*k_up[w,p])  for w in 1:nS, p in 1:nP))



        @objective(model, Min, 
            sum(pfDA[f] * Cf[f] for f in 1:nF) 
            + sum(ppDA[p] * Cp[p] for p in 1:nP) 
            + sum(pgDA[r] * Cr[r] for r in 1:nR)
            + sum(prob[w]*Cp[p] * ppB[w,p]  for w in 1:nS,   p in 1:nP) 
            + sum(prob[w]*Cf[f] * (pfB_up[w,f] - pfB_down[w,f]) for w in 1:nS, f in 1:nF)
            + bilinear_cost +  penalty_cost)

        optimize!(model)
        
        ########## Model solved ############
        if termination_status(model) != MOI.OPTIMAL
            #println("Model did not solve to optimality.")
            #println("primal", primal_status(model))
            #println("dual", dual_status(model))
            return nothing
        end

        Lambda_t = dual.(t_constraint)
        delta_1 = 0.01

        zk_val_up = value.(z_up) .* JuMP.value.(k_up)
        zk_val_down = value.(z_down) .* JuMP.value.(k_down)

        pct_diff_val_up = 100 .* (zk_val_up .- JuMP.value.(ppB_up)) ./ (JuMP.value.(ppB_up).*value.(s_up) .+ eps())
        pct_diff_val_down = 100 .* (zk_val_down .- JuMP.value.(ppB_down)) ./ (JuMP.value.(ppB_down).*value.(s_down) .+ eps())

        diff_matrix_up = zk_val_up .- JuMP.value.(ppB_up).*value.(s_up)
        diff_matrix_down = zk_val_down .- JuMP.value.(ppB_down).*value.(s_down)

        push!(iteration_history, (
            iter = iter,
            ppB = deepcopy(JuMP.value.(ppB)),
            z_up = deepcopy(value.(z_up)),
            k_up = deepcopy(JuMP.value.(k_up)),
            ppB_up = deepcopy(JuMP.value.(ppB_up)),
            diff_up_pct = deepcopy(pct_diff_val_up),
            z_down = deepcopy(value.(z_down)),
            k_down = deepcopy(JuMP.value.(k_down)),
            ppB_down = deepcopy(JuMP.value.(ppB_down)),
            diff_down_pct = deepcopy(pct_diff_val_down),
            ppDA = deepcopy(JuMP.value.(ppDA)),
            pfDA = deepcopy(JuMP.value.(pfDA)),
            pgDA = deepcopy(JuMP.value.(pgDA)),
            pfB_up = deepcopy(JuMP.value.(pfB_up)),
            pfB_down = deepcopy(JuMP.value.(pfB_down)),
            s_up = deepcopy(value.(s_up)),
            s_down = deepcopy(value.(s_down)),
            Imb = deepcopy(value.(Imb)),
            Imb_abs = deepcopy(value.(Imb_abs)),
            y = deepcopy(value.(y)),
            obj = objective_value(model),
            bilinear_cost = deepcopy(value.(bilinear_cost)),
            penalty_cost = deepcopy(value.(penalty_cost)),
            rf_up = deepcopy(rf_up),
            rf_down = deepcopy(rf_down), 
            diff_up = deepcopy(diff_matrix_up),
            diff_down = deepcopy(diff_matrix_down),
            t = deepcopy(JuMP.value.(t))
        ))

        ######### Convergence check ############
        
        
        step_norm = sqrt(
            norm(value.(k_up) - k_up_prev)^2 +
            #norm(value.(z_up) - z_up_prev)^2 +
            #norm(value.(ppB_up) - ppB_up_prev)^2 +
            norm(value.(s_up) - s_up_prev)^2 +
            norm(value.(k_down) - k_down_prev)^2 +
            ##norm(value.(z_down) - z_down_prev)^2 +
            #norm(value.(ppB_down) - ppB_down_prev)^2 +
            norm(value.(s_down) - s_down_prev)^2
        )
        # step_norm < tol &&
        if step_norm < tol && all(abs.(diff_matrix_up) .< tol) && all(abs.(diff_matrix_down) .< tol) && all(value.(t) .< tol)
            if any(abs.(diff_matrix_up) .> tol)
                #println("⚠️ Some values in diff_matrix_up exceed tolerance $tol")
            end
            if any(abs.(diff_matrix_up) .> tol)
                #println("⚠️ Some values in diff_matrix_down exceed tolerance $tol")
            end
            #   println("TA Bin Converged after  $iter iterations")
            return iteration_history
        end

        inv_step = 1.0 / max(step_norm, eps())
        rk = min(inv_step, norm(Lambda_t, 1) + delta_1)
        if (penalty_t < rk)
            #println("penalty updated!")
            penalty_t = min(max_penalty, penalty_t * mult_factor)
        end

        if iter % 100 == 0  # print every 100 iterations
            #println("Iteration $iter summary:")
            #println("  UP diag: step_norm=$(step_norm_up), inv_step=$(inv_step_up), norm(Lambda_up,1)=$(norm(Lambda_up,1)), penalty_up=$(penalty_up), violation_up=$(violation_up)")
            #println("  DOWN diag: step_norm=$(step_norm_down), inv_step=$(inv_step_down), norm(Lambda_down,1)=$(norm(Lambda_down,1)), penalty_down=$(penalty_down), violation_down=$(violation_down)")
        end

        ############# Update k, z, ppB and y ############
        k_up_prev    .= (1-alpha).*k_up_prev    .+ alpha.*value.(k_up)
        z_up_prev    .= (1-alpha).*z_up_prev    .+ alpha.*value.(z_up)
        ppB_up_prev  .= (1-alpha).*ppB_up_prev  .+ alpha.*value.(ppB_up)
        s_up_prev    .= (1-alpha).*s_up_prev    .+ alpha.*value.(s_up)

        k_down_prev    .= (1-alpha).*k_down_prev    .+ alpha.*value.(k_down)
        z_down_prev    .= (1-alpha).*z_down_prev    .+ alpha.*value.(z_down)
        ppB_down_prev  .= (1-alpha).*ppB_down_prev  .+ alpha.*value.(ppB_down)
        s_down_prev    .= (1-alpha).*s_down_prev    .+ alpha.*value.(s_down)

    end

    #println("TA Did not converge after $max_iter iterations") 
    return iteration_history
end
function TA_bin_GDCA_k(Pw, nS, prob, Q_up, Q_down,  rf_up, rf_down,
    NI_up, NI_down,
    k_up_prev, k_down_prev, z_up_prev, 
    z_down_prev, ppB_up_prev, ppB_down_prev,
    s_up_prev, s_down_prev;
    max_iter=2000, tol=0.01, penalty_t=500,
    M=300, alpha = 0.9)

    # make sure we’re working on private copies
    k_up_prev   = copy(k_up_prev)
    k_down_prev = copy(k_down_prev)
    z_up_prev   = copy(z_up_prev)
    z_down_prev = copy(z_down_prev)
    ppB_up_prev = copy(ppB_up_prev)
    ppB_down_prev = copy(ppB_down_prev)
    s_up_prev   = copy(s_up_prev)
    s_down_prev = copy(s_down_prev)
    
    iteration_history = []

    ############## Penalty update (safer) ##############
    # parameters you can tune:
    max_penalty = 1e5          # hard cap to avoid blow-up
    mult_factor = 1.5         # multiplicative increase (use >1) OR set add_sigma>0 to use additive
    iteration_history = []

    for iter in 1:max_iter
        ########### Model ############
        model = Model(HiGHS.Optimizer)
        set_silent(model)


        @variable(model, 0 <= ppDA[p=1:nP] <= Pmaxp[p])
        @variable(model, 0 <= pfDA[f=1:nF] <= Pmaxf[f])
        @variable(model, 0 <= pgDA[r=1:nR] <= Pmaxg[r])

        @variable(model, ppB[1:nS, p=1:nP])
        @variable(model, 0 <= ppB_up[1:nS, p=1:nP] <= Pmaxp[p])
        @variable(model, 0 <= ppB_down[1:nS, p=1:nP] <= Pmaxp[p])

        @variable(model, 0 <= pfB_up[1:nS, f=1:nF] <= rf_up[f])
        @variable(model, 0 <= pfB_down[1:nS, f=1:nF] <= rf_down[f])

        @variable(model, 0 <= k_up[1:nS, p=1:nP] <= 1)
        @variable(model, 0 <= k_down[1:nS, p=1:nP] <= 1)

        @variable(model, 0 <= z_up[1:nS] <= NI_up)
        @variable(model, 0 <= z_down[1:nS] <= NI_down)
        
        @variable(model, δ_up[1:nS, 1:nP] >= 0)
        @variable(model, δ_down[1:nS, 1:nP] >= 0)
        @variable(model, t[1:nS, 1:nP] >= 0)
    
        @variable(model, Imb[1:nS])
        @variable(model, Imb_abs[1:nS] >= 0)
        @variable(model, s_up[1:nS], Bin)
        @variable(model, s_down[1:nS], Bin)

        @variable(model, y[1:nS, 1:nP], Bin)

        #@variable(model, dev_up[1:nS, 1:nP] >= 0)
        #@variable(model, dev_down[1:nS, 1:nP] >= 0)
        @variable(model, k_max[1:nS] >= 0)

        ############ Constraints ############
     
        ##########Polluter ###############
        @constraint(model, [w in 1:nS, p in 1:nP],
            ppDA[p] + ppB[w, p] - Pw[w, p] == 0)
        
        @constraint(model, [w in 1:nS, p in 1:nP],
            ppB[w, p] == ppB_down[w, p] - ppB_up[w, p])

        @constraint(model, [w in 1:nS, p in 1:nP], ppB_down[w, p] <= M*(1 - y[w,p]))
        @constraint(model, [w in 1:nS, p in 1:nP], ppB_up[w, p] <= M*y[w,p])
    
        @constraint(model, [w=1:nS], z_up[w] - sum(ppB_up[w, p] for p in 1:nP) == 0)
        @constraint(model, [w=1:nS], z_down[w] - sum(ppB_down[w, p] for p in 1:nP) == 0)
        
        @expression(model, linearized_diff_up[w=1:nS, p=1:nP],
            k_up_prev[w,p] * z_up[w] + z_up_prev[w]*k_up[w,p] - z_up_prev[w]*k_up_prev[w,p] - (s_up_prev[w]*ppB_up[w,p] + s_up[w]*ppB_up_prev[w,p] - s_up_prev[w]*ppB_up_prev[w, p]))

        @expression(model, linearized_diff_down[w=1:nS, p=1:nP],
            k_down_prev[w,p] * z_down[w] + z_down_prev[w]*k_down[w,p] - z_down_prev[w]*k_down_prev[w,p] - (s_down_prev[w]*ppB_down[w,p] + s_down[w]*ppB_down_prev[w,p] - s_down_prev[w]*ppB_down_prev[w, p]))

        @constraint(model, [w=1:nS, p=1:nP], δ_up[w,p] >= linearized_diff_up[w,p])
        @constraint(model, [w=1:nS, p=1:nP], δ_up[w,p] >= -linearized_diff_up[w,p])

        @constraint(model, [w=1:nS, p=1:nP], δ_down[w,p] >= linearized_diff_down[w,p])
        @constraint(model, [w=1:nS, p=1:nP], δ_down[w,p] >= -linearized_diff_down[w,p])

        @constraint(model, [w=1:nS], sum(k_up[w, p] for p in 1:nP) == s_up[w])
        @constraint(model, [w=1:nS], sum(k_down[w, p] for p in 1:nP) == s_down[w])
 
        t_constraint = @constraint(model, [w=1:nS, p=1:nP], δ_down[w,p] + δ_up[w,p] <= t[w,p])

        @constraint(model, [w=1:nS, p=1:nP], k_max[w] >= k_up[w, p])
        @constraint(model, [w=1:nS, p=1:nP], k_max[w] >= k_down[w, p])
        #prod_scale = [Pw[w,p] / (sum(Pw[w,j] for j in 1:nP)) for w in 1:nS, p in 1:nP] 
        #@constraint(model, [w=1:nS, p=1:nP], dev_up[w,p] >= ppB_up[w,p] - z_up[w]*prod_scale[w,p])
        #@constraint(model, [w=1:nS, p=1:nP], dev_up[w,p] >= -(ppB_up[w,p] - z_up[w]*prod_scale[w,p]))

        #@constraint(model, [w=1:nS, p=1:nP], dev_down[w,p] >= ppB_down[w,p] - z_down[w]*prod_scale[w,p])
        #@constraint(model, [w=1:nS, p=1:nP], dev_down[w,p] >= -(ppB_down[w,p] - z_down[w]*prod_scale[w,p]))
        #@constraint(model, [w=1:nS], Imb[w] == sum(ppB[w, p] for p in 1:nP)) 
        @constraint(model, [w=1:nS], Imb[w] == z_down[w] - z_up[w]) 
        @constraint(model, [w=1:nS], Imb_abs[w] >= Imb[w])
        @constraint(model, [w=1:nS], Imb_abs[w] >= -Imb[w])
        @constraint(model, [w=1:nS], s_up[w] + s_down[w] <= 1)

        @constraint(model, [w=1:nS], Imb_abs[w] - Imb[w] <= 2*NI_up*(s_up[w]))
        @constraint(model, [w=1:nS], Imb_abs[w] + Imb[w] <= 2*NI_down*(s_down[w]))

        ########## Flex ###############

        @constraint(model,[f=1:nF], rf_down[f] <= pfDA[f])
        @constraint(model,[f=1:nF], pfDA[f] <= Pmaxf[f] - rf_up[f])
        
        @constraint(model, [w=1:nS, f=1:nF], pfB_down[w, f] <= rf_down[f]*s_down[w])
        @constraint(model, [w=1:nS, f=1:nF], pfB_up[w, f] <= rf_up[f]*s_up[w])

        ############### DA market clearing ###############

        demand_con = @constraint(model, 
            -D + sum(ppDA[p] for p in 1:nP) + sum(pfDA[f] for f in 1:nF) + sum(pgDA[r] for r in 1:nR) == 0)

        ############### Balancing market clearing ###############

        balancing_cons = @constraint(model, [w in 1:nS], 
            prob[w] * (sum(ppB[w, p] for p in 1:nP) - (sum(pfB_down[w,f] - pfB_up[w,f] for f in 1:nF))) == 0)

        ############ Objective ############

        #Penalty for the bilinear terms
        @expression(model, bilinear_cost, sum(prob[w] * penalty_t *t[w,p] for w in 1:nS, p in 1:nP))

        #Cost of the pollution - with one s 
        @expression(model, penalty_cost, sum(prob[w] *(Q_down*k_down[w,p] + Q_up*k_up[w,p])  for w in 1:nS, p in 1:nP))

        #Equality constraints costs
       # @expression(model, spread_penalty, sum(prob[w] * (dev_up[w,p] + dev_down[w,p]) for w in 1:nS, p in 1:nP))

        #equality constraints costs
        @expression(model, spread_penalty, 200*sum(k_max[w] for w in 1:nS))

        @objective(model, Min, 
            sum(pfDA[f] * Cf[f] for f in 1:nF) 
            + sum(ppDA[p] * Cp[p] for p in 1:nP) 
            + sum(pgDA[r] * Cr[r] for r in 1:nR)
            + sum(prob[w]*Cp[p] * ppB[w,p]  for w in 1:nS,   p in 1:nP) 
            + sum(prob[w]*Cf[f] * (pfB_up[w,f] - pfB_down[w,f]) for w in 1:nS, f in 1:nF)
            + bilinear_cost +  penalty_cost + spread_penalty)

        optimize!(model)
        
        ########## Model solved ############
        if termination_status(model) != MOI.OPTIMAL
            #println("Model did not solve to optimality.")
            #println("primal", primal_status(model))
            #println("dual", dual_status(model))
            return nothing
        end

        Lambda_t = dual.(t_constraint)
        delta_1 = 0.01

        zk_val_up = value.(z_up) .* JuMP.value.(k_up)
        zk_val_down = value.(z_down) .* JuMP.value.(k_down)

        pct_diff_val_up = 100 .* (zk_val_up .- JuMP.value.(ppB_up)) ./ (JuMP.value.(ppB_up).*value.(s_up) .+ eps())
        pct_diff_val_down = 100 .* (zk_val_down .- JuMP.value.(ppB_down)) ./ (JuMP.value.(ppB_down).*value.(s_down) .+ eps())

        diff_matrix_up = zk_val_up .- JuMP.value.(ppB_up).*value.(s_up)
        diff_matrix_down = zk_val_down .- JuMP.value.(ppB_down).*value.(s_down)

        push!(iteration_history, (
            iter = iter,
            ppB = deepcopy(JuMP.value.(ppB)),
            z_up = deepcopy(value.(z_up)),
            k_up = deepcopy(JuMP.value.(k_up)),
            ppB_up = deepcopy(JuMP.value.(ppB_up)),
            diff_up_pct = deepcopy(pct_diff_val_up),
            z_down = deepcopy(value.(z_down)),
            k_down = deepcopy(JuMP.value.(k_down)),
            ppB_down = deepcopy(JuMP.value.(ppB_down)),
            diff_down_pct = deepcopy(pct_diff_val_down),
            ppDA = deepcopy(JuMP.value.(ppDA)),
            pfDA = deepcopy(JuMP.value.(pfDA)),
            pgDA = deepcopy(JuMP.value.(pgDA)),
            pfB_up = deepcopy(JuMP.value.(pfB_up)),
            pfB_down = deepcopy(JuMP.value.(pfB_down)),
            s_up = deepcopy(value.(s_up)),
            s_down = deepcopy(value.(s_down)),
            Imb = deepcopy(value.(Imb)),
            Imb_abs = deepcopy(value.(Imb_abs)),
            y = deepcopy(value.(y)),
            obj = objective_value(model),
            bilinear_cost = deepcopy(value.(bilinear_cost)),
            penalty_cost = deepcopy(value.(penalty_cost)),
            spread_penalty = deepcopy(value.(spread_penalty)),
            obj_sans_imbalance = objective_value(model) - value.(spread_penalty),
            rf_up = deepcopy(rf_up),
            rf_down = deepcopy(rf_down), 
            diff_up = deepcopy(diff_matrix_up),
            diff_down = deepcopy(diff_matrix_down),
            t = deepcopy(JuMP.value.(t))
        ))

        ######### Convergence check ############
        
        
        step_norm = sqrt(
            norm(value.(k_up) - k_up_prev)^2 +
            #norm(value.(z_up) - z_up_prev)^2 +
            #norm(value.(ppB_up) - ppB_up_prev)^2 +
            norm(value.(s_up) - s_up_prev)^2 +
            norm(value.(k_down) - k_down_prev)^2 +
            ##norm(value.(z_down) - z_down_prev)^2 +
            #norm(value.(ppB_down) - ppB_down_prev)^2 +
            norm(value.(s_down) - s_down_prev)^2
        )
        # step_norm < tol &&
        if step_norm < tol && all(abs.(diff_matrix_up) .< tol) && all(abs.(diff_matrix_down) .< tol) && all(value.(t) .< tol)
            if any(abs.(diff_matrix_up) .> tol)
                #println("⚠️ Some values in diff_matrix_up exceed tolerance $tol")
            end
            if any(abs.(diff_matrix_up) .> tol)
                #println("⚠️ Some values in diff_matrix_down exceed tolerance $tol")
            end
            #   println("TA Bin Converged after  $iter iterations")
            return iteration_history
        end

        inv_step = 1.0 / max(step_norm, eps())
        rk = min(inv_step, norm(Lambda_t, 1) + delta_1)
        if (penalty_t < rk)
            #println("penalty updated!")
            penalty_t = min(max_penalty, penalty_t * mult_factor)
        end

        if iter % 100 == 0  # print every 100 iterations
            #println("Iteration $iter summary:")
            #println("  UP diag: step_norm=$(step_norm_up), inv_step=$(inv_step_up), norm(Lambda_up,1)=$(norm(Lambda_up,1)), penalty_up=$(penalty_up), violation_up=$(violation_up)")
            #println("  DOWN diag: step_norm=$(step_norm_down), inv_step=$(inv_step_down), norm(Lambda_down,1)=$(norm(Lambda_down,1)), penalty_down=$(penalty_down), violation_down=$(violation_down)")
        end

        ############# Update k, z, ppB and y ############
        k_up_prev    .= (1-alpha).*k_up_prev    .+ alpha.*value.(k_up)
        z_up_prev    .= (1-alpha).*z_up_prev    .+ alpha.*value.(z_up)
        ppB_up_prev  .= (1-alpha).*ppB_up_prev  .+ alpha.*value.(ppB_up)
        s_up_prev    .= (1-alpha).*s_up_prev    .+ alpha.*value.(s_up)

        k_down_prev    .= (1-alpha).*k_down_prev    .+ alpha.*value.(k_down)
        z_down_prev    .= (1-alpha).*z_down_prev    .+ alpha.*value.(z_down)
        ppB_down_prev  .= (1-alpha).*ppB_down_prev  .+ alpha.*value.(ppB_down)
        s_down_prev    .= (1-alpha).*s_down_prev    .+ alpha.*value.(s_down)

    end

    #println("TA Did not converge after $max_iter iterations") 
    return iteration_history
end

function TA_bin_GDCA_no_step(Pw, nS, prob, Q_up, Q_down,  rf_up, rf_down,
    NI_up, NI_down,
    k_up_prev, k_down_prev, z_up_prev, 
    z_down_prev, ppB_up_prev, ppB_down_prev,
    s_up_prev, s_down_prev;
    max_iter=2000, tol=0.01, penalty_t=500,
    M=300, alpha = 0.9)

    # make sure we’re working on private copies
    k_up_prev   = copy(k_up_prev)
    k_down_prev = copy(k_down_prev)
    z_up_prev   = copy(z_up_prev)
    z_down_prev = copy(z_down_prev)
    ppB_up_prev = copy(ppB_up_prev)
    ppB_down_prev = copy(ppB_down_prev)
    s_up_prev   = copy(s_up_prev)
    s_down_prev = copy(s_down_prev)
    
    iteration_history = []

    ############## Penalty update (safer) ##############
    # parameters you can tune:
    max_penalty = 1e5          # hard cap to avoid blow-up
    mult_factor = 1.5         # multiplicative increase (use >1) OR set add_sigma>0 to use additive
    iteration_history = []

    for iter in 1:max_iter
        ########### Model ############
        model = Model(HiGHS.Optimizer)
        set_silent(model)

        @variable(model, 0 <= ppDA[p=1:nP] <= Pmaxp[p])
        @variable(model, 0 <= pfDA[f=1:nF] <= Pmaxf[f])
        @variable(model, 0 <= pgDA[r=1:nR] <= Pmaxg[r])

        @variable(model, ppB[1:nS, p=1:nP])
        @variable(model, 0 <= ppB_up[1:nS, p=1:nP] <= Pmaxp[p])
        @variable(model, 0 <= ppB_down[1:nS, p=1:nP] <= Pmaxp[p])

        @variable(model, 0 <= pfB_up[1:nS, f=1:nF] <= rf_up[f])
        @variable(model, 0 <= pfB_down[1:nS, f=1:nF] <= rf_down[f])

        @variable(model, 0 <= k_up[1:nS, p=1:nP] <= 1)
        @variable(model, 0 <= k_down[1:nS, p=1:nP] <= 1)

        @variable(model, 0 <= z_up[1:nS] <= NI_up)
        @variable(model, 0 <= z_down[1:nS] <= NI_down)
        
        @variable(model, δ_up[1:nS, 1:nP] >= 0)
        @variable(model, δ_down[1:nS, 1:nP] >= 0)
        @variable(model, t[1:nS, 1:nP] >= 0)
    
        @variable(model, Imb[1:nS])
        @variable(model, Imb_abs[1:nS] >= 0)
        @variable(model, s_up[1:nS], Bin)
        @variable(model, s_down[1:nS], Bin)

        @variable(model, y[1:nS, 1:nP], Bin)

        ############ Constraints ############
     
        ##########Polluter ###############
        @constraint(model, [w in 1:nS, p in 1:nP],
            ppDA[p] + ppB[w, p] - Pw[w, p] == 0)
        
        @constraint(model, [w in 1:nS, p in 1:nP],
            ppB[w, p] == ppB_down[w, p] - ppB_up[w, p])

        @constraint(model, [w in 1:nS, p in 1:nP], ppB_down[w, p] <= M*(1 - y[w,p]))
        @constraint(model, [w in 1:nS, p in 1:nP], ppB_up[w, p] <= M*y[w,p])
    
        @constraint(model, [w=1:nS], z_up[w] - sum(ppB_up[w, p] for p in 1:nP) == 0)
        @constraint(model, [w=1:nS], z_down[w] - sum(ppB_down[w, p] for p in 1:nP) == 0)
        
        @expression(model, linearized_diff_up[w=1:nS, p=1:nP],
            k_up_prev[w,p] * z_up[w] + z_up_prev[w]*k_up[w,p] - z_up_prev[w]*k_up_prev[w,p] - (s_up_prev[w]*ppB_up[w,p] + s_up[w]*ppB_up_prev[w,p] - s_up_prev[w]*ppB_up_prev[w, p]))

        @expression(model, linearized_diff_down[w=1:nS, p=1:nP],
            k_down_prev[w,p] * z_down[w] + z_down_prev[w]*k_down[w,p] - z_down_prev[w]*k_down_prev[w,p] - (s_down_prev[w]*ppB_down[w,p] + s_down[w]*ppB_down_prev[w,p] - s_down_prev[w]*ppB_down_prev[w, p]))

        @constraint(model, [w=1:nS, p=1:nP], δ_up[w,p] >= linearized_diff_up[w,p])
        @constraint(model, [w=1:nS, p=1:nP], δ_up[w,p] >= -linearized_diff_up[w,p])

        @constraint(model, [w=1:nS, p=1:nP], δ_down[w,p] >= linearized_diff_down[w,p])
        @constraint(model, [w=1:nS, p=1:nP], δ_down[w,p] >= -linearized_diff_down[w,p])

        @constraint(model, [w=1:nS], sum(k_up[w, p] for p in 1:nP) == s_up[w])
        @constraint(model, [w=1:nS], sum(k_down[w, p] for p in 1:nP) == s_down[w])
 
        t_constraint = @constraint(model, [w=1:nS, p=1:nP], δ_down[w,p] + δ_up[w,p] <= t[w,p])

        #@constraint(model, [w=1:nS], Imb[w] == sum(ppB[w, p] for p in 1:nP)) 
        @constraint(model, [w=1:nS], Imb[w] == z_down[w] - z_up[w]) 
        @constraint(model, [w=1:nS], Imb_abs[w] >= Imb[w])
        @constraint(model, [w=1:nS], Imb_abs[w] >= -Imb[w])
        @constraint(model, [w=1:nS], s_up[w] + s_down[w] <= 1)

        @constraint(model, [w=1:nS], Imb_abs[w] - Imb[w] <= 2*NI_up*(s_up[w]))
        @constraint(model, [w=1:nS], Imb_abs[w] + Imb[w] <= 2*NI_down*(s_down[w]))

        ########## Flex ###############

        @constraint(model,[f=1:nF], rf_down[f] <= pfDA[f])
        @constraint(model,[f=1:nF], pfDA[f] <= Pmaxf[f] - rf_up[f])
        
        @constraint(model, [w=1:nS, f=1:nF], pfB_down[w, f] <= rf_down[f]*s_down[w])
        @constraint(model, [w=1:nS, f=1:nF], pfB_up[w, f] <= rf_up[f]*s_up[w])

        ############### DA market clearing ###############

        demand_con = @constraint(model, 
            -D + sum(ppDA[p] for p in 1:nP) + sum(pfDA[f] for f in 1:nF) + sum(pgDA[r] for r in 1:nR) == 0)

        ############### Balancing market clearing ###############

        balancing_cons = @constraint(model, [w in 1:nS], 
            prob[w] * (sum(ppB[w, p] for p in 1:nP) - (sum(pfB_down[w,f] - pfB_up[w,f] for f in 1:nF))) == 0)

        ############ Objective ############

        #Penalty for the bilinear terms
        @expression(model, bilinear_cost, sum(prob[w] * penalty_t *t[w,p] for w in 1:nS, p in 1:nP))

        #Cost of the pollution - with one s 
        @expression(model, penalty_cost, sum(prob[w] *(Q_down*k_down[w,p] + Q_up*k_up[w,p])  for w in 1:nS, p in 1:nP))

        @objective(model, Min, 
            sum(pfDA[f] * Cf[f] for f in 1:nF) 
            + sum(ppDA[p] * Cp[p] for p in 1:nP) 
            + sum(pgDA[r] * Cr[r] for r in 1:nR)
            + sum(prob[w]*Cp[p] * ppB[w,p]  for w in 1:nS,   p in 1:nP) 
            + sum(prob[w]*Cf[f] * (pfB_up[w,f] - pfB_down[w,f]) for w in 1:nS, f in 1:nF)
            + bilinear_cost +  penalty_cost)

        optimize!(model)
        
        ########## Model solved ############
        if termination_status(model) != MOI.OPTIMAL
            #println("Model did not solve to optimality.")
            #println("primal", primal_status(model))
            #println("dual", dual_status(model))
            return nothing
        end

        Lambda_t = dual.(t_constraint)
        delta_1 = 0.01

        zk_val_up = value.(z_up) .* JuMP.value.(k_up)
        zk_val_down = value.(z_down) .* JuMP.value.(k_down)

        pct_diff_val_up = 100 .* (zk_val_up .- JuMP.value.(ppB_up)) ./ (JuMP.value.(ppB_up).*value.(s_up) .+ eps())
        pct_diff_val_down = 100 .* (zk_val_down .- JuMP.value.(ppB_down)) ./ (JuMP.value.(ppB_down).*value.(s_down) .+ eps())

        diff_matrix_up = zk_val_up .- JuMP.value.(ppB_up).*value.(s_up)
        diff_matrix_down = zk_val_down .- JuMP.value.(ppB_down).*value.(s_down)

        push!(iteration_history, (
            iter = iter,
            ppB = deepcopy(JuMP.value.(ppB)),
            z_up = deepcopy(value.(z_up)),
            k_up = deepcopy(JuMP.value.(k_up)),
            ppB_up = deepcopy(JuMP.value.(ppB_up)),
            diff_up_pct = deepcopy(pct_diff_val_up),
            z_down = deepcopy(value.(z_down)),
            k_down = deepcopy(JuMP.value.(k_down)),
            ppB_down = deepcopy(JuMP.value.(ppB_down)),
            diff_down_pct = deepcopy(pct_diff_val_down),
            ppDA = deepcopy(JuMP.value.(ppDA)),
            pfDA = deepcopy(JuMP.value.(pfDA)),
            pgDA = deepcopy(JuMP.value.(pgDA)),
            pfB_up = deepcopy(JuMP.value.(pfB_up)),
            pfB_down = deepcopy(JuMP.value.(pfB_down)),
            s_up = deepcopy(value.(s_up)),
            s_down = deepcopy(value.(s_down)),
            Imb = deepcopy(value.(Imb)),
            Imb_abs = deepcopy(value.(Imb_abs)),
            y = deepcopy(value.(y)),
            obj = objective_value(model),
            bilinear_cost = deepcopy(value.(bilinear_cost)),
            penalty_cost = deepcopy(value.(penalty_cost)),
            rf_up = deepcopy(rf_up),
            rf_down = deepcopy(rf_down), 
            diff_up = deepcopy(diff_matrix_up),
            diff_down = deepcopy(diff_matrix_down),
            t = deepcopy(JuMP.value.(t))
        ))

        ######### Convergence check ############
        
        
        step_norm = sqrt(
            norm(value.(k_up) - k_up_prev)^2 +
            #norm(value.(z_up) - z_up_prev)^2 +
            #norm(value.(ppB_up) - ppB_up_prev)^2 +
            norm(value.(s_up) - s_up_prev)^2 +
            norm(value.(k_down) - k_down_prev)^2 +
            ##norm(value.(z_down) - z_down_prev)^2 +
            #norm(value.(ppB_down) - ppB_down_prev)^2 +
            norm(value.(s_down) - s_down_prev)^2
        )
        # step_norm < tol &&
        if all(abs.(diff_matrix_up) .< tol) && all(abs.(diff_matrix_down) .< tol) && all(value.(t) .< tol)
            if any(abs.(diff_matrix_up) .> tol)
                #println("⚠️ Some values in diff_matrix_up exceed tolerance $tol")
            end
            if any(abs.(diff_matrix_up) .> tol)
                #println("⚠️ Some values in diff_matrix_down exceed tolerance $tol")
            end
            #   println("TA Bin Converged after  $iter iterations")
            return iteration_history
        end

        inv_step = 1.0 / max(step_norm, eps())
        rk = min(inv_step, norm(Lambda_t, 1) + delta_1)
        if (penalty_t < rk)
            #println("penalty updated!")
            penalty_t = min(max_penalty, penalty_t * mult_factor)
        end

        if iter % 100 == 0  # print every 100 iterations
            #println("Iteration $iter summary:")
            #println("  UP diag: step_norm=$(step_norm_up), inv_step=$(inv_step_up), norm(Lambda_up,1)=$(norm(Lambda_up,1)), penalty_up=$(penalty_up), violation_up=$(violation_up)")
            #println("  DOWN diag: step_norm=$(step_norm_down), inv_step=$(inv_step_down), norm(Lambda_down,1)=$(norm(Lambda_down,1)), penalty_down=$(penalty_down), violation_down=$(violation_down)")
        end

        ############# Update k, z, ppB and y ############
        k_up_prev    .= (1-alpha).*k_up_prev    .+ alpha.*value.(k_up)
        z_up_prev    .= (1-alpha).*z_up_prev    .+ alpha.*value.(z_up)
        ppB_up_prev  .= (1-alpha).*ppB_up_prev  .+ alpha.*value.(ppB_up)
        s_up_prev    .= (1-alpha).*s_up_prev    .+ alpha.*value.(s_up)

        k_down_prev    .= (1-alpha).*k_down_prev    .+ alpha.*value.(k_down)
        z_down_prev    .= (1-alpha).*z_down_prev    .+ alpha.*value.(z_down)
        ppB_down_prev  .= (1-alpha).*ppB_down_prev  .+ alpha.*value.(ppB_down)
        s_down_prev    .= (1-alpha).*s_down_prev    .+ alpha.*value.(s_down)

    end

    #println("TA Did not converge after $max_iter iterations") 
    return iteration_history
end

function TA_bin_LP_GDCA(history, Pw, nS, prob, Q_up, Q_down, 
    rf_up, rf_down, 
    NI_up, NI_down, 
    k_up_prev, k_down_prev, 
    z_up_prev, z_down_prev, ppB_up_prev, ppB_down_prev; 
     max_iter=2000, tol=0.01, penalty_t=1000,
     M=300, alpha= 1.0)
        
    # make sure we’re working on private copies
    k_up_prev   = copy(k_up_prev)
    k_down_prev = copy(k_down_prev)
    z_up_prev   = copy(z_up_prev)
    z_down_prev = copy(z_down_prev)
    ppB_up_prev = copy(ppB_up_prev)
    ppB_down_prev = copy(ppB_down_prev)
    iteration_history = []
    
    ############## Penalty update (safer) ##############
    # parameters you can tune:
    max_penalty = 1e5          # hard cap to avoid blow-up
    mult_factor = 1.05         # multiplicative increase (use >1) OR set add_sigma>0 to use additive
    iteration_history = []

    for iter in 1:max_iter
        s_up   = history.s_up
        s_down = history.s_down
        y       = history.y

        ########### Model ############
        model = Model(HiGHS.Optimizer)
        set_silent(model)

        @variable(model, 0 <= ppDA[p=1:nP] <= Pmaxp[p])
        @variable(model, 0 <= pfDA[f =1:nF] <= Pmaxf[f])
        @variable(model, 0 <= pgDA[r = 1:nR] <= Pmaxg[r])

        @variable(model, ppB[1:nS, p=1:nP])
        @variable(model, 0 <= ppB_up[1:nS, p=1:nP] <= Pmaxp[p])
        @variable(model, 0 <= ppB_down[1:nS, p=1:nP] <= Pmaxp[p])

        @variable(model, 0 <= pfB_up[1:nS, f=1:nF] <= rf_up[f])
        @variable(model, 0 <= pfB_down[1:nS, f=1:nF] <= rf_down[f])

        @variable(model, 0 <= k_up[1:nS, p=1:nP] <= 1)
        @variable(model, 0 <= k_down[1:nS, p=1:nP] <= 1)

        @variable(model, 0 <= z_up[1:nS] <= NI_up)
        @variable(model, 0 <= z_down[1:nS] <= NI_down)
        
        @variable(model, δ_up[1:nS, 1:nP] >= 0)
        @variable(model, δ_down[1:nS, 1:nP] >= 0)
        
        @variable(model, t[1:nS, 1:nP] >= 0)

        @variable(model, Imb[1:nS])
        @variable(model, Imb_abs[1:nS] >= 0)
        

        ############ Constraints ############
     
        ##########Polluter ###############
        @constraint(model, [w in 1:nS, p in 1:nP],
            ppDA[p] + ppB[w, p] - Pw[w, p] == 0)
        
        @constraint(model, [w in 1:nS, p in 1:nP],
            ppB[w, p] == ppB_down[w, p] - ppB_up[w, p])

        @constraint(model, [w in 1:nS, p in 1:nP], ppB_down[w, p] <= M*(1 - y[w,p]))
        @constraint(model, [w in 1:nS, p in 1:nP], ppB_up[w, p] <= M*y[w,p])
    
        @constraint(model, [w=1:nS], z_up[w] - sum(ppB_up[w, p] for p in 1:nP) == 0)
        @constraint(model, [w=1:nS], z_down[w] - sum(ppB_down[w, p] for p in 1:nP) == 0)
        
        @expression(model, linearized_diff_up[w=1:nS, p=1:nP],
            k_up_prev[w,p] * z_up[w] + z_up_prev[w]*k_up[w,p] - z_up_prev[w]*k_up_prev[w,p] - (s_up[w]*ppB_up[w,p] ))

        @expression(model, linearized_diff_down[w=1:nS, p=1:nP],
            k_down_prev[w,p] * z_down[w] + z_down_prev[w]*k_down[w,p] - z_down_prev[w]*k_down_prev[w,p] - (s_down[w]*ppB_down[w,p] ))

        @constraint(model, [w=1:nS, p=1:nP], δ_up[w,p] >= linearized_diff_up[w,p])
        @constraint(model, [w=1:nS, p=1:nP], δ_up[w,p] >= -linearized_diff_up[w,p])

        @constraint(model, [w=1:nS, p=1:nP], δ_down[w,p] >= linearized_diff_down[w,p])
        @constraint(model, [w=1:nS, p=1:nP], δ_down[w,p] >= -linearized_diff_down[w,p])

        @constraint(model, [w=1:nS], sum(k_up[w, p] for p in 1:nP) == s_up[w])
        @constraint(model, [w=1:nS], sum(k_down[w, p] for p in 1:nP) == s_down[w])
 
        t_constraint = @constraint(model, [w=1:nS, p=1:nP], δ_down[w,p] + δ_up[w,p] <= t[w,p])

        @constraint(model, [w=1:nS], Imb[w] == sum(ppB[w, p] for p in 1:nP)) 
        @constraint(model, [w=1:nS], Imb[w] == z_down[w] - z_up[w]) 
        @constraint(model, [w=1:nS], Imb_abs[w] >= Imb[w])
        @constraint(model, [w=1:nS], Imb_abs[w] >= -Imb[w])

        @constraint(model, [w=1:nS], Imb_abs[w] - Imb[w] <= 2*NI_up*(s_up[w]))
        @constraint(model, [w=1:nS], Imb_abs[w] + Imb[w] <= 2*NI_down*(s_down[w]))

        ########## Flex ###############

        @constraint(model,[f=1:nF], rf_down[f] <= pfDA[f])
        @constraint(model,[f=1:nF], pfDA[f] <= Pmaxf[f] - rf_up[f])
        
        @constraint(model, [w=1:nS, f=1:nF], pfB_down[w, f] <= rf_down[f]*s_down[w])
        @constraint(model, [w=1:nS, f=1:nF], pfB_up[w, f] <= rf_up[f]*s_up[w])

        ############### DA market clearing ###############

        demand_con = @constraint(model, 
            -D + sum(ppDA[p] for p in 1:nP) + sum(pfDA[f] for f in 1:nF) + sum(pgDA[r] for r in 1:nR) == 0)

        ############### Balancing market clearing ###############

        balancing_cons = @constraint(model, [w in 1:nS], 
            prob[w] * (sum(ppB[w, p] for p in 1:nP) - sum(pfB_down[w,f] - pfB_up[w,f] for f in 1:nF)) == 0)

        ############ Objective ############

        #Penalty for the bilinear terms
        @expression(model, bilinear_cost, sum(prob[w] * penalty_t *t[w,p] for w in 1:nS, p in 1:nP))

        #Cost of the pollution - with one s 
        @expression(model, penalty_cost, sum(prob[w] *((Q_down)*k_down[w,p] + (Q_up)*k_up[w,p])  for w in 1:nS, p in 1:nP))

    

        @expression(model, ppDA_cost, sum(ppDA[p] * Cp[p] for p in 1:nP))
        @expression(model, pfDA_cost, sum(pfDA[f] * Cf[f] for f in 1:nF))
        @expression(model, pgDA_cost, sum(pgDA[r] * Cr[r] for r in 1:nR))

        @expression(model, ppDA_single_cost[p=1:nP], ppDA[p] * Cp[p])

        @expression(model, ppB_single_cost[w=1:nS, p=1:nP], prob[w] * Cp[p] * ppB[w,p])

        @expression(model, pfB_single_cost[w=1:nS, f=1:nF], prob[w] * Cf[f] * (pfB_up[w, f] - pfB_down[w, f]))

        @objective(model, Min,
            pfDA_cost + pgDA_cost + ppDA_cost
            + sum(ppB_single_cost[w, p] for w in 1:nS, p in 1:nP) 
            + sum(pfB_single_cost[w, f] for w in 1:nS, f in 1:nF)
            +  bilinear_cost +  penalty_cost )

        optimize!(model)
        
        ########## Model solved ############
        if termination_status(model) != MOI.OPTIMAL
            #println("Model did not solve to optimality.")
            return nothing
        end

        zk_val_up = value.(z_up) .* JuMP.value.(k_up)
        zk_val_down = value.(z_down) .* JuMP.value.(k_down)

        pct_diff_val_up = 100 .* (zk_val_up .- JuMP.value.(ppB_up)) ./ (JuMP.value.(ppB_up).*value.(s_up) .+ eps())
        pct_diff_val_down = 100 .* (zk_val_down .- JuMP.value.(ppB_down)) ./ (JuMP.value.(ppB_down).*value.(s_down) .+ eps())

        Lambda_t = dual.(t_constraint)
        delta_1 = 0.01

        diff_matrix_up = zk_val_up .- JuMP.value.(ppB_up).*s_up
        diff_matrix_down = zk_val_down .- JuMP.value.(ppB_down).*s_down

        push!(iteration_history, (
            iter = iter,
            ppB = deepcopy(JuMP.value.(ppB)),
            z_up = deepcopy(value.(z_up)),
            k_up = deepcopy(JuMP.value.(k_up)),
            ppB_up = deepcopy(JuMP.value.(ppB_up)),
            diff_up_pct = deepcopy(pct_diff_val_up),
            z_down = deepcopy(value.(z_down)),
            k_down = deepcopy(JuMP.value.(k_down)),
            ppB_down = deepcopy(JuMP.value.(ppB_down)),
            diff_down_pct = deepcopy(pct_diff_val_down),
            ppDA = deepcopy(JuMP.value.(ppDA)),
            pfDA = deepcopy(JuMP.value.(pfDA)),
            pgDA = deepcopy(JuMP.value.(pgDA)),
            pfB_up = deepcopy(JuMP.value.(pfB_up)),
            pfB_down = deepcopy(JuMP.value.(pfB_down)),
            s_up = deepcopy(value.(s_up)),
            s_down = deepcopy(value.(s_down)),
            Imb = deepcopy(value.(Imb)),
            Imb_abs = deepcopy(value.(Imb_abs)),
            y = deepcopy(value.(y)),
            obj = objective_value(model),
            ppDA_cost = deepcopy(value.(ppDA_cost)),
            ppDA_single_cost = deepcopy(value.(ppDA_single_cost)),
            pfDA_cost = deepcopy(JuMP.value.(pfDA_cost)),
            pgDA_cost = deepcopy(value.(pgDA_cost)),
            #ppB_cost = deepcopy(value.(ppB_cost)),
            ppB_single_cost = deepcopy(JuMP.value.(ppB_single_cost)),
            #pf1B_cost = deepcopy(value.(pf1B_cost)),
            pfB_single_cost = deepcopy(JuMP.value.(pfB_single_cost)),
            bilinear_cost = deepcopy(value.(bilinear_cost)),
            penalty_cost = deepcopy(value.(penalty_cost)),
            Lambda_DA = dual(demand_con),
            Lambda_B = dual.(balancing_cons),
            t = deepcopy(value.(t)),
            rf_up = deepcopy(rf_up),
            rf_down = deepcopy(rf_down), 
            diff_up = deepcopy(diff_matrix_up),
            diff_down = deepcopy(diff_matrix_down)
        ))


        ######### Convergence check ############
        

        
        step_norm = sqrt(
            norm(value.(k_up) - k_up_prev)^2 +
            #norm(value.(z_up) - z_up_prev)^2 +
            #norm(value.(ppB_up) - ppB_up_prev)^2 +
            norm(value.(k_down) - k_down_prev)^2 
            #norm(value.(z_down) - z_down_prev)^2 +
            #norm(value.(ppB_down) - ppB_down_prev)^2
        )

        if step_norm < tol && all(abs.(diff_matrix_up) .< tol) && all(abs.(diff_matrix_down) .< tol) && all(value.(t) .< tol)
            if any(abs.(diff_matrix_up) .> tol)
                #println("⚠️ Some values in diff_matrix_up exceed tolerance $tol")
            end
            if any(abs.(diff_matrix_up) .> tol)
                #println("⚠️ Some values in diff_matrix_down exceed tolerance $tol")
            end
            #println("TA bin LP Converged after  $iter iterations")
            return iteration_history
        end

        inv_step = 1.0 / max(step_norm, eps())
        rk = min(inv_step, norm(Lambda_t, 1) + delta_1)
        #violation_up = maximum(abs.(diff_matrix_up))
        if (penalty_t < rk)
            #println("Penalty update in iter $iter")
            penalty_t = min(max_penalty, penalty_t * mult_factor)
        end

        if iter % 1000 == 0  # print every 100 iterations
            #println("Iteration $iter summary:")
        end

        ############# Update k, z, ppB and y ############
        k_up_prev    .= (1-alpha).*k_up_prev    .+ alpha.*value.(k_up)
        z_up_prev    .= (1-alpha).*z_up_prev    .+ alpha.*value.(z_up)
        ppB_up_prev  .= (1-alpha).*ppB_up_prev  .+ alpha.*value.(ppB_up)
   
        k_down_prev    .= (1-alpha).*k_down_prev    .+ alpha.*value.(k_down)
        z_down_prev    .= (1-alpha).*z_down_prev    .+ alpha.*value.(z_down)
        ppB_down_prev  .= (1-alpha).*ppB_down_prev  .+ alpha.*value.(ppB_down)
    
    end

    #println("TA Did not converge after $max_iter iterations") 
    return iteration_history
end

function TA_bin_LP_GDCA_k(history, Pw, nS, prob, Q_up, Q_down, 
    rf_up, rf_down, 
    NI_up, NI_down, 
    k_up_prev, k_down_prev, 
    z_up_prev, z_down_prev, ppB_up_prev, ppB_down_prev; 
     max_iter=2000, tol=0.01, penalty_t=1000,
     M=300, alpha= 1.0)
        
    # make sure we’re working on private copies
    k_up_prev   = copy(k_up_prev)
    k_down_prev = copy(k_down_prev)
    z_up_prev   = copy(z_up_prev)
    z_down_prev = copy(z_down_prev)
    ppB_up_prev = copy(ppB_up_prev)
    ppB_down_prev = copy(ppB_down_prev)
    iteration_history = []
    
    ############## Penalty update (safer) ##############
    # parameters you can tune:
    max_penalty = 1e5          # hard cap to avoid blow-up
    mult_factor = 1.05         # multiplicative increase (use >1) OR set add_sigma>0 to use additive
    iteration_history = []

    for iter in 1:max_iter
        s_up   = history.s_up
        s_down = history.s_down
        y       = history.y

        ########### Model ############
        model = Model(HiGHS.Optimizer)
        set_silent(model)

        @variable(model, 0 <= ppDA[p=1:nP] <= Pmaxp[p])
        @variable(model, 0 <= pfDA[f =1:nF] <= Pmaxf[f])
        @variable(model, 0 <= pgDA[r = 1:nR] <= Pmaxg[r])

        @variable(model, ppB[1:nS, p=1:nP])
        @variable(model, 0 <= ppB_up[1:nS, p=1:nP] <= Pmaxp[p])
        @variable(model, 0 <= ppB_down[1:nS, p=1:nP] <= Pmaxp[p])

        @variable(model, 0 <= pfB_up[1:nS, f=1:nF] <= rf_up[f])
        @variable(model, 0 <= pfB_down[1:nS, f=1:nF] <= rf_down[f])

        @variable(model, 0 <= k_up[1:nS, p=1:nP] <= 1)
        @variable(model, 0 <= k_down[1:nS, p=1:nP] <= 1)

        @variable(model, 0 <= z_up[1:nS] <= NI_up)
        @variable(model, 0 <= z_down[1:nS] <= NI_down)
        
        @variable(model, δ_up[1:nS, 1:nP] >= 0)
        @variable(model, δ_down[1:nS, 1:nP] >= 0)
        
        @variable(model, t[1:nS, 1:nP] >= 0)

        @variable(model, Imb[1:nS])
        @variable(model, Imb_abs[1:nS] >= 0)
        @variable(model, k_max[1:nS] >= 0)

        ############ Constraints ############
     
        ##########Polluter ###############
        @constraint(model, [w in 1:nS, p in 1:nP],
            ppDA[p] + ppB[w, p] - Pw[w, p] == 0)
        
        @constraint(model, [w in 1:nS, p in 1:nP],
            ppB[w, p] == ppB_down[w, p] - ppB_up[w, p])

        @constraint(model, [w in 1:nS, p in 1:nP], ppB_down[w, p] <= M*(1 - y[w,p]))
        @constraint(model, [w in 1:nS, p in 1:nP], ppB_up[w, p] <= M*y[w,p])
    
        @constraint(model, [w=1:nS], z_up[w] - sum(ppB_up[w, p] for p in 1:nP) == 0)
        @constraint(model, [w=1:nS], z_down[w] - sum(ppB_down[w, p] for p in 1:nP) == 0)
        
        @expression(model, linearized_diff_up[w=1:nS, p=1:nP],
            k_up_prev[w,p] * z_up[w] + z_up_prev[w]*k_up[w,p] - z_up_prev[w]*k_up_prev[w,p] - (s_up[w]*ppB_up[w,p] ))

        @expression(model, linearized_diff_down[w=1:nS, p=1:nP],
            k_down_prev[w,p] * z_down[w] + z_down_prev[w]*k_down[w,p] - z_down_prev[w]*k_down_prev[w,p] - (s_down[w]*ppB_down[w,p] ))

        @constraint(model, [w=1:nS, p=1:nP], δ_up[w,p] >= linearized_diff_up[w,p])
        @constraint(model, [w=1:nS, p=1:nP], δ_up[w,p] >= -linearized_diff_up[w,p])

        @constraint(model, [w=1:nS, p=1:nP], δ_down[w,p] >= linearized_diff_down[w,p])
        @constraint(model, [w=1:nS, p=1:nP], δ_down[w,p] >= -linearized_diff_down[w,p])

        @constraint(model, [w=1:nS], sum(k_up[w, p] for p in 1:nP) == s_up[w])
        @constraint(model, [w=1:nS], sum(k_down[w, p] for p in 1:nP) == s_down[w])
 
        t_constraint = @constraint(model, [w=1:nS, p=1:nP], δ_down[w,p] + δ_up[w,p] <= t[w,p])

        @constraint(model, [w=1:nS, p=1:nP], k_max[w] >= k_up[w,p])
        @constraint(model, [w=1:nS, p=1:nP], k_max[w] >= k_down[w,p])


        @constraint(model, [w=1:nS], Imb[w] == sum(ppB[w, p] for p in 1:nP)) 
        @constraint(model, [w=1:nS], Imb[w] == z_down[w] - z_up[w]) 
        @constraint(model, [w=1:nS], Imb_abs[w] >= Imb[w])
        @constraint(model, [w=1:nS], Imb_abs[w] >= -Imb[w])

        @constraint(model, [w=1:nS], Imb_abs[w] - Imb[w] <= 2*NI_up*(s_up[w]))
        @constraint(model, [w=1:nS], Imb_abs[w] + Imb[w] <= 2*NI_down*(s_down[w]))

        ########## Flex ###############

        @constraint(model,[f=1:nF], rf_down[f] <= pfDA[f])
        @constraint(model,[f=1:nF], pfDA[f] <= Pmaxf[f] - rf_up[f])
        
        @constraint(model, [w=1:nS, f=1:nF], pfB_down[w, f] <= rf_down[f]*s_down[w])
        @constraint(model, [w=1:nS, f=1:nF], pfB_up[w, f] <= rf_up[f]*s_up[w])

        ############### DA market clearing ###############

        demand_con = @constraint(model, 
            -D + sum(ppDA[p] for p in 1:nP) + sum(pfDA[f] for f in 1:nF) + sum(pgDA[r] for r in 1:nR) == 0)

        ############### Balancing market clearing ###############

        balancing_cons = @constraint(model, [w in 1:nS], 
            prob[w] * (sum(ppB[w, p] for p in 1:nP) - sum(pfB_down[w,f] - pfB_up[w,f] for f in 1:nF)) == 0)

        ############ Objective ############

        #Penalty for the bilinear terms
        @expression(model, bilinear_cost, sum(prob[w] * penalty_t *t[w,p] for w in 1:nS, p in 1:nP))

        #Cost of the pollution - with one s 
        @expression(model, penalty_cost, sum(prob[w] *((Q_down)*k_down[w,p] + (Q_up)*k_up[w,p])  for w in 1:nS, p in 1:nP))

        #Equality constraints costs
        #@expression(model, spread_penalty, sum(prob[w] * (dev_up[w,p] + dev_down[w,p]) for w in 1:nS, p in 1:nP))
        @expression(model, spread_penalty, 200*sum(k_max[w] for w in 1:nS))

        @expression(model, ppDA_cost, sum(ppDA[p] * Cp[p] for p in 1:nP))
        @expression(model, pfDA_cost, sum(pfDA[f] * Cf[f] for f in 1:nF))
        @expression(model, pgDA_cost, sum(pgDA[r] * Cr[r] for r in 1:nR))

        @expression(model, ppDA_single_cost[p=1:nP], ppDA[p] * Cp[p])

        @expression(model, ppB_single_cost[w=1:nS, p=1:nP], prob[w] * Cp[p] * ppB[w,p])

        @expression(model, pfB_single_cost[w=1:nS, f=1:nF], prob[w] * Cf[f] * (pfB_up[w, f] - pfB_down[w, f]))

        @objective(model, Min,
            pfDA_cost + pgDA_cost + ppDA_cost
            + sum(ppB_single_cost[w, p] for w in 1:nS, p in 1:nP) 
            + sum(pfB_single_cost[w, f] for w in 1:nS, f in 1:nF)
            +  bilinear_cost +  penalty_cost + spread_penalty)

        optimize!(model)
        
        ########## Model solved ############
        if termination_status(model) != MOI.OPTIMAL
            #println("Model did not solve to optimality.")
            return nothing
        end

        zk_val_up = value.(z_up) .* JuMP.value.(k_up)
        zk_val_down = value.(z_down) .* JuMP.value.(k_down)

        pct_diff_val_up = 100 .* (zk_val_up .- JuMP.value.(ppB_up)) ./ (JuMP.value.(ppB_up).*value.(s_up) .+ eps())
        pct_diff_val_down = 100 .* (zk_val_down .- JuMP.value.(ppB_down)) ./ (JuMP.value.(ppB_down).*value.(s_down) .+ eps())

        Lambda_t = dual.(t_constraint)
        delta_1 = 0.01

        diff_matrix_up = zk_val_up .- JuMP.value.(ppB_up).*s_up
        diff_matrix_down = zk_val_down .- JuMP.value.(ppB_down).*s_down

        push!(iteration_history, (
            iter = iter,
            ppB = deepcopy(JuMP.value.(ppB)),
            z_up = deepcopy(value.(z_up)),
            k_up = deepcopy(JuMP.value.(k_up)),
            ppB_up = deepcopy(JuMP.value.(ppB_up)),
            diff_up_pct = deepcopy(pct_diff_val_up),
            z_down = deepcopy(value.(z_down)),
            k_down = deepcopy(JuMP.value.(k_down)),
            ppB_down = deepcopy(JuMP.value.(ppB_down)),
            diff_down_pct = deepcopy(pct_diff_val_down),
            ppDA = deepcopy(JuMP.value.(ppDA)),
            pfDA = deepcopy(JuMP.value.(pfDA)),
            pgDA = deepcopy(JuMP.value.(pgDA)),
            pfB_up = deepcopy(JuMP.value.(pfB_up)),
            pfB_down = deepcopy(JuMP.value.(pfB_down)),
            s_up = deepcopy(value.(s_up)),
            s_down = deepcopy(value.(s_down)),
            Imb = deepcopy(value.(Imb)),
            Imb_abs = deepcopy(value.(Imb_abs)),
            y = deepcopy(value.(y)),
            obj = objective_value(model),
            ppDA_cost = deepcopy(value.(ppDA_cost)),
            ppDA_single_cost = deepcopy(value.(ppDA_single_cost)),
            pfDA_cost = deepcopy(JuMP.value.(pfDA_cost)),
            pgDA_cost = deepcopy(value.(pgDA_cost)),
            #ppB_cost = deepcopy(value.(ppB_cost)),
            ppB_single_cost = deepcopy(JuMP.value.(ppB_single_cost)),
            #pf1B_cost = deepcopy(value.(pf1B_cost)),
            pfB_single_cost = deepcopy(JuMP.value.(pfB_single_cost)),
            bilinear_cost = deepcopy(value.(bilinear_cost)),
            penalty_cost = deepcopy(value.(penalty_cost)),
            spread_penalty = deepcopy(value.(spread_penalty)),
            obj_sans_imbalance = objective_value(model) - value.(spread_penalty),
            Lambda_DA = dual(demand_con),
            Lambda_B = dual.(balancing_cons),
            t = deepcopy(value.(t)),
            rf_up = deepcopy(rf_up),
            rf_down = deepcopy(rf_down), 
            diff_up = deepcopy(diff_matrix_up),
            diff_down = deepcopy(diff_matrix_down)
        ))


        ######### Convergence check ############
        

        
        step_norm = sqrt(
            norm(value.(k_up) - k_up_prev)^2 +
            #norm(value.(z_up) - z_up_prev)^2 +
            #norm(value.(ppB_up) - ppB_up_prev)^2 +
            norm(value.(k_down) - k_down_prev)^2 
            #norm(value.(z_down) - z_down_prev)^2 +
            #norm(value.(ppB_down) - ppB_down_prev)^2
        )

        if step_norm < tol && all(abs.(diff_matrix_up) .< tol) && all(abs.(diff_matrix_down) .< tol) && all(value.(t) .< tol)
            if any(abs.(diff_matrix_up) .> tol)
                #println("⚠️ Some values in diff_matrix_up exceed tolerance $tol")
            end
            if any(abs.(diff_matrix_up) .> tol)
                #println("⚠️ Some values in diff_matrix_down exceed tolerance $tol")
            end
            #println("TA bin LP Converged after  $iter iterations")
            return iteration_history
        end

        inv_step = 1.0 / max(step_norm, eps())
        rk = min(inv_step, norm(Lambda_t, 1) + delta_1)
        #violation_up = maximum(abs.(diff_matrix_up))
        if (penalty_t < rk)
            #println("Penalty update in iter $iter")
            penalty_t = min(max_penalty, penalty_t * mult_factor)
        end

        if iter % 1000 == 0  # print every 100 iterations
            #println("Iteration $iter summary:")
        end

        ############# Update k, z, ppB and y ############
        k_up_prev    .= (1-alpha).*k_up_prev    .+ alpha.*value.(k_up)
        z_up_prev    .= (1-alpha).*z_up_prev    .+ alpha.*value.(z_up)
        ppB_up_prev  .= (1-alpha).*ppB_up_prev  .+ alpha.*value.(ppB_up)
   
        k_down_prev    .= (1-alpha).*k_down_prev    .+ alpha.*value.(k_down)
        z_down_prev    .= (1-alpha).*z_down_prev    .+ alpha.*value.(z_down)
        ppB_down_prev  .= (1-alpha).*ppB_down_prev  .+ alpha.*value.(ppB_down)
    
    end

    #println("TA Did not converge after $max_iter iterations") 
    return iteration_history
end

function TA_bin_LP_GDCA_no_step(history, Pw, nS, prob, Q_up, Q_down, 
    rf_up, rf_down, 
    NI_up, NI_down, 
    k_up_prev, k_down_prev, 
    z_up_prev, z_down_prev, ppB_up_prev, ppB_down_prev; 
     max_iter=2000, tol=0.01, penalty_t=1000,
     M=300, alpha= 1.0)
        
    # make sure we’re working on private copies
    k_up_prev   = copy(k_up_prev)
    k_down_prev = copy(k_down_prev)
    z_up_prev   = copy(z_up_prev)
    z_down_prev = copy(z_down_prev)
    ppB_up_prev = copy(ppB_up_prev)
    ppB_down_prev = copy(ppB_down_prev)
    iteration_history = []
    
    ############## Penalty update (safer) ##############
    # parameters you can tune:
    max_penalty = 1e5          # hard cap to avoid blow-up
    mult_factor = 1.05         # multiplicative increase (use >1) OR set add_sigma>0 to use additive
    iteration_history = []

    for iter in 1:max_iter
        s_up   = history.s_up
        s_down = history.s_down
        y       = history.y

        ########### Model ############
        model = Model(HiGHS.Optimizer)
        set_silent(model)

        @variable(model, 0 <= ppDA[p=1:nP] <= Pmaxp[p])
        @variable(model, 0 <= pfDA[f =1:nF] <= Pmaxf[f])
        @variable(model, 0 <= pgDA[r = 1:nR] <= Pmaxg[r])

        @variable(model, ppB[1:nS, p=1:nP])
        @variable(model, 0 <= ppB_up[1:nS, p=1:nP] <= Pmaxp[p])
        @variable(model, 0 <= ppB_down[1:nS, p=1:nP] <= Pmaxp[p])

        @variable(model, 0 <= pfB_up[1:nS, f=1:nF] <= rf_up[f])
        @variable(model, 0 <= pfB_down[1:nS, f=1:nF] <= rf_down[f])

        @variable(model, 0 <= k_up[1:nS, p=1:nP] <= 1)
        @variable(model, 0 <= k_down[1:nS, p=1:nP] <= 1)

        @variable(model, 0 <= z_up[1:nS] <= NI_up)
        @variable(model, 0 <= z_down[1:nS] <= NI_down)
        
        @variable(model, δ_up[1:nS, 1:nP] >= 0)
        @variable(model, δ_down[1:nS, 1:nP] >= 0)
        
        @variable(model, t[1:nS, 1:nP] >= 0)

        @variable(model, Imb[1:nS])
        @variable(model, Imb_abs[1:nS] >= 0)

        ############ Constraints ############
     
        ##########Polluter ###############
        @constraint(model, [w in 1:nS, p in 1:nP],
            ppDA[p] + ppB[w, p] - Pw[w, p] == 0)
        
        @constraint(model, [w in 1:nS, p in 1:nP],
            ppB[w, p] == ppB_down[w, p] - ppB_up[w, p])

        @constraint(model, [w in 1:nS, p in 1:nP], ppB_down[w, p] <= M*(1 - y[w,p]))
        @constraint(model, [w in 1:nS, p in 1:nP], ppB_up[w, p] <= M*y[w,p])
    
        @constraint(model, [w=1:nS], z_up[w] - sum(ppB_up[w, p] for p in 1:nP) == 0)
        @constraint(model, [w=1:nS], z_down[w] - sum(ppB_down[w, p] for p in 1:nP) == 0)
        
        @expression(model, linearized_diff_up[w=1:nS, p=1:nP],
            k_up_prev[w,p] * z_up[w] + z_up_prev[w]*k_up[w,p] - z_up_prev[w]*k_up_prev[w,p] - (s_up[w]*ppB_up[w,p] ))

        @expression(model, linearized_diff_down[w=1:nS, p=1:nP],
            k_down_prev[w,p] * z_down[w] + z_down_prev[w]*k_down[w,p] - z_down_prev[w]*k_down_prev[w,p] - (s_down[w]*ppB_down[w,p] ))

        @constraint(model, [w=1:nS, p=1:nP], δ_up[w,p] >= linearized_diff_up[w,p])
        @constraint(model, [w=1:nS, p=1:nP], δ_up[w,p] >= -linearized_diff_up[w,p])

        @constraint(model, [w=1:nS, p=1:nP], δ_down[w,p] >= linearized_diff_down[w,p])
        @constraint(model, [w=1:nS, p=1:nP], δ_down[w,p] >= -linearized_diff_down[w,p])

        @constraint(model, [w=1:nS], sum(k_up[w, p] for p in 1:nP) == s_up[w])
        @constraint(model, [w=1:nS], sum(k_down[w, p] for p in 1:nP) == s_down[w])
 
        t_constraint = @constraint(model, [w=1:nS, p=1:nP], δ_down[w,p] + δ_up[w,p] <= t[w,p])

        @constraint(model, [w=1:nS], Imb[w] == sum(ppB[w, p] for p in 1:nP)) 
        @constraint(model, [w=1:nS], Imb[w] == z_down[w] - z_up[w]) 
        @constraint(model, [w=1:nS], Imb_abs[w] >= Imb[w])
        @constraint(model, [w=1:nS], Imb_abs[w] >= -Imb[w])

        @constraint(model, [w=1:nS], Imb_abs[w] - Imb[w] <= 2*NI_up*(s_up[w]))
        @constraint(model, [w=1:nS], Imb_abs[w] + Imb[w] <= 2*NI_down*(s_down[w]))

        ########## Flex ###############

        @constraint(model,[f=1:nF], rf_down[f] <= pfDA[f])
        @constraint(model,[f=1:nF], pfDA[f] <= Pmaxf[f] - rf_up[f])
        
        @constraint(model, [w=1:nS, f=1:nF], pfB_down[w, f] <= rf_down[f]*s_down[w])
        @constraint(model, [w=1:nS, f=1:nF], pfB_up[w, f] <= rf_up[f]*s_up[w])

        ############### DA market clearing ###############

        demand_con = @constraint(model, 
            -D + sum(ppDA[p] for p in 1:nP) + sum(pfDA[f] for f in 1:nF) + sum(pgDA[r] for r in 1:nR) == 0)

        ############### Balancing market clearing ###############

        balancing_cons = @constraint(model, [w in 1:nS], 
            prob[w] * (sum(ppB[w, p] for p in 1:nP) - sum(pfB_down[w,f] - pfB_up[w,f] for f in 1:nF)) == 0)

        ############ Objective ############

        #Penalty for the bilinear terms
        @expression(model, bilinear_cost, sum(prob[w] * penalty_t *t[w,p] for w in 1:nS, p in 1:nP))

        #Cost of the pollution - with one s 
        @expression(model, penalty_cost, sum(prob[w] *((Q_down)*k_down[w,p] + (Q_up)*k_up[w,p])  for w in 1:nS, p in 1:nP))

        @expression(model, ppDA_cost, sum(ppDA[p] * Cp[p] for p in 1:nP))
        @expression(model, pfDA_cost, sum(pfDA[f] * Cf[f] for f in 1:nF))
        @expression(model, pgDA_cost, sum(pgDA[r] * Cr[r] for r in 1:nR))

        @expression(model, ppDA_single_cost[p=1:nP], ppDA[p] * Cp[p])

        @expression(model, ppB_single_cost[w=1:nS, p=1:nP], prob[w] * Cp[p] * ppB[w,p])

        @expression(model, pfB_single_cost[w=1:nS, f=1:nF], prob[w] * Cf[f] * (pfB_up[w, f] - pfB_down[w, f]))

        @objective(model, Min,
            pfDA_cost + pgDA_cost + ppDA_cost
            + sum(ppB_single_cost[w, p] for w in 1:nS, p in 1:nP) 
            + sum(pfB_single_cost[w, f] for w in 1:nS, f in 1:nF)
            +  bilinear_cost +  penalty_cost)

        optimize!(model)
        
        ########## Model solved ############
        if termination_status(model) != MOI.OPTIMAL
            #println("Model did not solve to optimality.")
            return nothing
        end

        zk_val_up = value.(z_up) .* JuMP.value.(k_up)
        zk_val_down = value.(z_down) .* JuMP.value.(k_down)

        pct_diff_val_up = 100 .* (zk_val_up .- JuMP.value.(ppB_up)) ./ (JuMP.value.(ppB_up).*value.(s_up) .+ eps())
        pct_diff_val_down = 100 .* (zk_val_down .- JuMP.value.(ppB_down)) ./ (JuMP.value.(ppB_down).*value.(s_down) .+ eps())

        Lambda_t = dual.(t_constraint)
        delta_1 = 0.01

        diff_matrix_up = zk_val_up .- JuMP.value.(ppB_up).*s_up
        diff_matrix_down = zk_val_down .- JuMP.value.(ppB_down).*s_down

        push!(iteration_history, (
            iter = iter,
            ppB = deepcopy(JuMP.value.(ppB)),
            z_up = deepcopy(value.(z_up)),
            k_up = deepcopy(JuMP.value.(k_up)),
            ppB_up = deepcopy(JuMP.value.(ppB_up)),
            diff_up_pct = deepcopy(pct_diff_val_up),
            z_down = deepcopy(value.(z_down)),
            k_down = deepcopy(JuMP.value.(k_down)),
            ppB_down = deepcopy(JuMP.value.(ppB_down)),
            diff_down_pct = deepcopy(pct_diff_val_down),
            ppDA = deepcopy(JuMP.value.(ppDA)),
            pfDA = deepcopy(JuMP.value.(pfDA)),
            pgDA = deepcopy(JuMP.value.(pgDA)),
            pfB_up = deepcopy(JuMP.value.(pfB_up)),
            pfB_down = deepcopy(JuMP.value.(pfB_down)),
            s_up = deepcopy(value.(s_up)),
            s_down = deepcopy(value.(s_down)),
            Imb = deepcopy(value.(Imb)),
            Imb_abs = deepcopy(value.(Imb_abs)),
            y = deepcopy(value.(y)),
            obj = objective_value(model),
            ppDA_cost = deepcopy(value.(ppDA_cost)),
            ppDA_single_cost = deepcopy(value.(ppDA_single_cost)),
            pfDA_cost = deepcopy(JuMP.value.(pfDA_cost)),
            pgDA_cost = deepcopy(value.(pgDA_cost)),
            #ppB_cost = deepcopy(value.(ppB_cost)),
            ppB_single_cost = deepcopy(JuMP.value.(ppB_single_cost)),
            #pf1B_cost = deepcopy(value.(pf1B_cost)),
            pfB_single_cost = deepcopy(JuMP.value.(pfB_single_cost)),
            bilinear_cost = deepcopy(value.(bilinear_cost)),
            penalty_cost = deepcopy(value.(penalty_cost)),
            Lambda_DA = dual(demand_con),
            Lambda_B = dual.(balancing_cons),
            t = deepcopy(value.(t)),
            rf_up = deepcopy(rf_up),
            rf_down = deepcopy(rf_down), 
            diff_up = deepcopy(diff_matrix_up),
            diff_down = deepcopy(diff_matrix_down)
        ))


        ######### Convergence check ############
        

        
        step_norm = sqrt(
            norm(value.(k_up) - k_up_prev)^2 +
            #norm(value.(z_up) - z_up_prev)^2 +
            #norm(value.(ppB_up) - ppB_up_prev)^2 +
            norm(value.(k_down) - k_down_prev)^2 
            #norm(value.(z_down) - z_down_prev)^2 +
            #norm(value.(ppB_down) - ppB_down_prev)^2
        )
        #step_norm < tol &&
        if  all(abs.(diff_matrix_up) .< tol) && all(abs.(diff_matrix_down) .< tol) && all(value.(t) .< tol)
            if any(abs.(diff_matrix_up) .> tol)
                #println("⚠️ Some values in diff_matrix_up exceed tolerance $tol")
            end
            if any(abs.(diff_matrix_up) .> tol)
                #println("⚠️ Some values in diff_matrix_down exceed tolerance $tol")
            end
            #println("TA bin LP Converged after  $iter iterations")
            return iteration_history
        end

        inv_step = 1.0 / max(step_norm, eps())
        rk = min(inv_step, norm(Lambda_t, 1) + delta_1)
        #violation_up = maximum(abs.(diff_matrix_up))
        if (penalty_t < rk)
            #println("Penalty update in iter $iter")
            penalty_t = min(max_penalty, penalty_t * mult_factor)
        end

        if iter % 1000 == 0  # print every 100 iterations
            #println("Iteration $iter summary:")
        end

        ############# Update k, z, ppB and y ############
        k_up_prev    .= (1-alpha).*k_up_prev    .+ alpha.*value.(k_up)
        z_up_prev    .= (1-alpha).*z_up_prev    .+ alpha.*value.(z_up)
        ppB_up_prev  .= (1-alpha).*ppB_up_prev  .+ alpha.*value.(ppB_up)
   
        k_down_prev    .= (1-alpha).*k_down_prev    .+ alpha.*value.(k_down)
        z_down_prev    .= (1-alpha).*z_down_prev    .+ alpha.*value.(z_down)
        ppB_down_prev  .= (1-alpha).*ppB_down_prev  .+ alpha.*value.(ppB_down)
    
    end

    #println("TA Did not converge after $max_iter iterations") 
    return iteration_history
end


####################################################################################

                                    # BENCHMARK

################# Centralized Optimization 
function Benchmark(Pw, nS, prob, rf_up, rf_down, NI_up, NI_down)
    iteration_history = []
    ########### Model ############
    model = Model(HiGHS.Optimizer)
    set_silent(model)

    @variable(model, 0 <= ppDA[p=1:nP] <= Pmaxp[p])
    @variable(model, 0 <= pfDA[f = 1:nF] <= Pmaxf[f])
    @variable(model, 0 <= pgDA[r = 1:nR] <= Pmaxg[r])

    @variable(model, ppB[1:nS, p=1:nP])
    @variable(model, 0 <= ppB_up[1:nS, p=1:nP] <= Pmaxp[p])
    @variable(model, 0 <= ppB_down[1:nS, p=1:nP] <= Pmaxp[p])

    @variable(model, 0 <= pfB_up[1:nS, f=1:nF] <= rf_up[f])
    @variable(model, 0 <= pfB_down[1:nS, f=1:nF] <= rf_down[f])

    ############ Constraints ############
    ##########Polluter ###############
    @constraint(model, [w in 1:nS, p in 1:nP], ppDA[p] + ppB[w, p] - Pw[w, p] == 0)
    
    @constraint(model, [w in 1:nS, p in 1:nP], ppB[w, p] == ppB_down[w, p] - ppB_up[w, p])

    @constraint(model, [w in 1:nS], sum(ppB_up[w, :]) <= NI_up)
    @constraint(model, [w in 1:nS], sum(ppB_down[w, :]) <= NI_down)

    ########## Flex ###############

    @constraint(model,[f=1:nF], rf_down[f] <= pfDA[f])
    @constraint(model,[f=1:nF], pfDA[f] <= Pmaxf[f] - rf_up[f])

    ############### DA market clearing ###############

    demand_con = @constraint(model, 
    - D +(sum(ppDA[p] for p in 1:nP) + sum(pfDA[f] for f in 1:nF) + sum(pgDA[r] for r in 1:nR)) == 0)

    ############### Balancing market clearing ###############

    balancing_cons = @constraint(model, [w in 1:nS], 
        prob[w] * (sum(ppB[w, p] for p in 1:nP)- (sum(pfB_down[w,f] - pfB_up[w,f] for f in 1:nF))) == 0)

    ############ Objective ############

    @objective(model, Min, 
        sum(pfDA[f] * Cf[f] for f in 1:nF) + sum(ppDA[p] * Cp[p] for p in 1:nP) + sum(pgDA[r] * Cr[r] for r in 1:nR)
        + sum(prob[w]*Cf[f] * (pfB_up[w,f] - pfB_down[w,f]) for w in 1:nS, f in 1:nF)
        +sum(prob[w]*Cp[p] *(ppB[w,p]) for w in 1:nS,  p in 1:nP))

    optimize!(model)
    
    ########## Model solved ############
    if termination_status(model) != MOI.OPTIMAL
        #println("Model did not solve to optimality.")
        return nothing
    end
    
    push!(iteration_history, (
        ppB = deepcopy(JuMP.value.(ppB)),
        ppB_up = deepcopy(JuMP.value.(ppB_up)),
        ppB_down = deepcopy(JuMP.value.(ppB_down)),
        ppDA = deepcopy(JuMP.value.(ppDA)),
        pfDA = deepcopy(JuMP.value.(pfDA)),
        pgDA = deepcopy(JuMP.value.(pgDA)),
        pfB_up = deepcopy(JuMP.value.(pfB_up)),
        pfB_down = deepcopy(JuMP.value.(pfB_down)),
        obj = objective_value(model),
        Lambda_DA = dual(demand_con),
        Lambda_B = dual.(balancing_cons), 
        rf_up = deepcopy(rf_up),
        rf_down = deepcopy(rf_down)
            ))
    return iteration_history
end


function polluter_benchmark(Pw, nS, prob, pol, central_sol)
    iteration_history = []
    ########### Model ############
    model = Model(HiGHS.Optimizer)
    set_silent(model)

    @variable(model, 0 <= ppDA <= Pmaxp[pol])

    @variable(model, ppB[1:nS])
    @variable(model, 0 <= ppB_up[1:nS] <= Pmaxp[pol])
    @variable(model, 0 <= ppB_down[1:nS] <= Pmaxp[pol])

    ############ Constraints ############
    ##########Polluter ###############
    @constraint(model, [w in 1:nS],
        ppDA + ppB[w] - Pw[w, pol] == 0)

    ############ Objective ############
    @objective(model, Max, 
        ppDA * (central_sol.Lambda_DA - Cp[pol])
        +sum(prob[w]*(central_sol.Lambda_B[w] -  Cp[pol])*ppB[w] for w in 1:nS))

    optimize!(model)
    
    ########## Model solved ############
    if termination_status(model) != MOI.OPTIMAL
       # println("Model did not solve to optimality.")
        return nothing
    end
    
    push!(iteration_history, (
        Gen = pol,
        ppB = deepcopy(value.(ppB)),
        ppB_up = deepcopy(value.(ppB_up)),
        ppB_down = deepcopy(value.(ppB_down)),
        ppDA = deepcopy(value(ppDA)),
        obj = objective_value(model),
    ))
    return iteration_history
end


################# Flexible's Individual Optimization 
function flex_benchmark(f, central_sol, rf_up, rf_down)
    
    iteration_history = []

    ########### Model ############
    model = Model(HiGHS.Optimizer)
    set_silent(model)
    @variable(model, 0 <= pfDA)

    @variable(model, 0 <= pfB_up[1:nS] <= rf_up[f])
    @variable(model, 0 <= pfB_down[1:nS] <= rf_down[f])

    ############ Constraints ############        

    ########## Flex ###############

    @constraint(model, rf_down[f] <= pfDA)
    @constraint(model, pfDA <= Pmaxf[f] - rf_up[f])

    ############ Objective ############
    #println(central_sol.Lambda_B)

    @objective(model, Max, 
        pfDA * (central_sol.Lambda_DA - Cf[f])
        + sum(prob[w]*(central_sol.Lambda_B[w] -  Cf[f])*(pfB_up[w]-pfB_down[w]) for w in 1:nS))

    optimize!(model)
    
    ########## Model solved ############
    if termination_status(model) != MOI.OPTIMAL
        #println("Model did not solve to optimality.")
        return nothing
    end

    push!(iteration_history, (
        pfDA = deepcopy(value(pfDA)),
        pfB_up = deepcopy(value.(pfB_up)),
        pfB_down = deepcopy(value.(pfB_down)),
        obj = objective_value(model),
    ))

    return iteration_history
end

################# Regular's Individual Optimization 
function reg_opt(r, central_sol)
    
    iteration_history = []

    ########### Model ############
    model = Model(HiGHS.Optimizer)
    set_silent(model)
    @variable(model, 0 <= pgDA <= Pmaxg[r])

    ############ Constraints ############        

    ############ Objective ############

    @objective(model, Max, 
        pgDA * (central_sol.Lambda_DA - Cr[r]))

    optimize!(model)
    
    ########## Model solved ############
    if termination_status(model) != MOI.OPTIMAL
        #println("Model did not solve to optimality.")
        return nothing
    end

    push!(iteration_history, (
        pgDA = deepcopy(value(pgDA)),
        obj = objective_value(model)
    ))

    return iteration_history
end


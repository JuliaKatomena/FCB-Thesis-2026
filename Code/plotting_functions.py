import numpy as np
import pandas as pd
from io import StringIO
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from matplotlib.patches import Patch
from matplotlib.lines import Line2D
from mpl_toolkits.axes_grid1.inset_locator import inset_axes, mark_inset
import matplotlib.colors as mcolors
import re

#### profits plotting functions
def stacked_posneg_matrix(bench, model, polluters, markets, participant_type="Polluters",
                          figsize = (10,8), fontsize = 17, ylim = (-50, 5000), ylabel = "Profit [€]"):
    # u helps for the tunning of the y axis
    nP, nM = bench.shape
    x = np.arange(nP)
    width = 0.35

    fig, ax = plt.subplots(figsize = figsize)

    # colors for markets
    blue = ["#08306b", "#3f8fc4", "#2257a5"]
    orange = ["#EE7518", "#EBA067", "#D05C03"]

    def plot_stack(mat, xpos, colors):
        bot_pos = np.zeros(nP)
        bot_neg = np.zeros(nP)

        for m in range(mat.shape[1]):
            vals = mat[:, m]

            pos = np.where(vals > 0, vals, 0)
            neg = np.where(vals < 0, vals, 0)

            # positive stack
            ax.bar(xpos, pos, width, bottom=bot_pos, color=colors[m])
            bot_pos += pos

            # negative stack
            ax.bar(xpos, neg, width, bottom=bot_neg, color=colors[m])
            bot_neg += neg

        return bot_pos + bot_neg  # net profit

    # plot left (benchmark)
    net_bench = plot_stack(bench, x - width/2, blue)
    # plot right (model)
    net_model = plot_stack(model, x + width/2, orange)

    # net profit dots
    ax.scatter(x - width/2, net_bench, color="black", s=20, label="BM – Net profit")
    ax.scatter(x + width/2, net_model, color="black", s=20, label="FCB – Net profit")
    print(min(np.min(net_bench), np.min(net_model)))
    ax.set_xticks(x)
    ax.set_xticklabels(polluters, rotation=45, ha='right', fontsize = fontsize)
    ax.set_ylim(ylim)
    ax.set_ylabel(ylabel, fontsize = fontsize)
    ax.tick_params(axis="y", labelsize=fontsize, rotation = 45)
    #ax.set_title(f"Profit decomposition per {participant_type}")
    ax.axhline(0, color="black", linewidth=0.8)
    #ax.set_ylim(min(np.min(net_bench), np.min(net_model)) * 1.1, max(np.max(net_bench), np.max(net_model)) * 1.1)

    # --- Create legend manually ---
    legend_handles = []

 # Benchmark colors
    if participant_type == "Polluters":
        for i, market in enumerate(markets[:2]):
            legend_handles.append(Patch(facecolor=blue[i], label=f"BM – {market}"))
    else:
        for i, market in enumerate(markets):
            legend_handles.append(Patch(facecolor=blue[i], label=f"BM – {market}"))

    # Model colors
    for i, market in enumerate(markets):
        legend_handles.append(Patch(facecolor=orange[i], label=f"FCB – {market}"))

    # Add net profit dot
    legend_handles.append(Line2D([0], [0], marker='o', color='black', linestyle='', markersize=6, label='Net profit'))

    ax.legend(handles=legend_handles, loc='upper left', frameon=True, fontsize = fontsize-4)
    plt.grid(True, alpha = 0.5)
    plt.tight_layout()
    
    return fig

def stacked_posneg_matrix(bench, model, polluters, markets, participant_type="Polluters",
                          figsize = (10,8), fontsize = 17, ylim = (-50, 5000), ylabel = "Profit [€]"):
    # u helps for the tunning of the y axis
    nP, nM = bench.shape
    x = np.arange(nP)
    width = 0.35

    fig, ax = plt.subplots(figsize = figsize)

    # colors for markets
    blue = ["#08306b", "#3f8fc4", "#2257a5"]
    orange = ["#EE7518", "#EBA067", "#D05C03"]

    def plot_stack(mat, xpos, colors):
        bot_pos = np.zeros(nP)
        bot_neg = np.zeros(nP)

        for m in range(mat.shape[1]):
            vals = mat[:, m]

            pos = np.where(vals > 0, vals, 0)
            neg = np.where(vals < 0, vals, 0)

            # positive stack
            ax.bar(xpos, pos, width, bottom=bot_pos, color=colors[m])
            bot_pos += pos

            # negative stack
            ax.bar(xpos, neg, width, bottom=bot_neg, color=colors[m])
            bot_neg += neg

        return bot_pos + bot_neg  # net profit

    # plot left (benchmark)
    net_bench = plot_stack(bench, x - width/2, blue)
    # plot right (model)
    net_model = plot_stack(model, x + width/2, orange)

    # net profit dots
    ax.scatter(x - width/2, net_bench, color="black", s=20, label="BM – Net profit")
    ax.scatter(x + width/2, net_model, color="black", s=20, label="FCB – Net profit")
    print(min(np.min(net_bench), np.min(net_model)))
    ax.set_xticks(x)
    ax.set_xticklabels(polluters, rotation=45, ha='right', fontsize = fontsize)
    ax.set_ylim(ylim)
    ax.set_ylabel(ylabel, fontsize = fontsize)
    ax.tick_params(axis="y", labelsize=fontsize, rotation = 45)
    #ax.set_title(f"Profit decomposition per {participant_type}")
    ax.axhline(0, color="black", linewidth=0.8)
    #ax.set_ylim(min(np.min(net_bench), np.min(net_model)) * 1.1, max(np.max(net_bench), np.max(net_model)) * 1.1)

    # --- Create legend manually ---
    legend_handles = []

 # Benchmark colors
    if participant_type == "Polluters":
        for i, market in enumerate(markets[:2]):
            legend_handles.append(Patch(facecolor=blue[i], label=f"BM – {market}"))
    else:
        for i, market in enumerate(markets):
            legend_handles.append(Patch(facecolor=blue[i], label=f"BM – {market}"))

    # Model colors
    for i, market in enumerate(markets):
        legend_handles.append(Patch(facecolor=orange[i], label=f"FCB – {market}"))

    # Add net profit dot
    legend_handles.append(Line2D([0], [0], marker='o', color='black', linestyle='', markersize=6, label='Net profit'))

    ax.legend(handles=legend_handles, loc='lower left', frameon=True, fontsize = fontsize-4)
    plt.grid(True, alpha = 0.5)
    plt.tight_layout()
    
    return fig

#plots power imbalance up/down for multiple models
def plot_power_imbalances_updown_sep_models(
    deviations_up_list,     # list of np.ndarray (nP,), all >=0
    deviations_down_list,   # list of np.ndarray (nP,), all >=0
    NI_up,                  # list of floats per model
    NI_down,                # list of floats per model
    polluters,
    model_names=None,
    base_color="Blues"
):
    """
    Multi-model up/down power imbalances.
    Now takes deviations_up/deviations_down directly (both positive arrays).
    Plots 2 subplots per model (up/down), side by side.
    """
    n_models = len(deviations_up_list)
    nP = len(polluters)
    print(type(polluters))
    print(type(nP))
    print("nP value:", nP)
    # Colors for polluters
    cmap = cm.get_cmap(base_color, nP)
    bar_colors = [cmap(i) for i in range(nP)]

    fig, axes = plt.subplots(
        1, 2 * n_models, figsize=(4 * n_models, 5), squeeze=False
    )

    if model_names is None:
        model_names = [f"Model {i+1}" for i in range(n_models)]

    for m in range(n_models):

        dev_up = np.asarray(deviations_up_list[m]).reshape(-1)
        dev_down = np.asarray(deviations_down_list[m]).reshape(-1)

        sys_imbalance = dev_down.sum() - dev_up.sum()   # positive = excess, negative = deficit

        # --- UP REGULATION subplot ---
        # (previously negative deviations)
        ax = axes[0, 2 * m]
        bottom = 0
        for i in range(nP):
            ax.bar(0, - dev_up[i], width=0.8, bottom=bottom, color=bar_colors[i])
            bottom += - dev_up[i]

        ax.axhline(NI_up, color='green', linestyle='--', linewidth=1)
        if sys_imbalance < 0:
            ax.axhline(-sys_imbalance, color='orange', linestyle='-', linewidth=1.5)

        ax.set_ylim(0, max(NI_up, NI_down)*1.1)
        ax.set_xticks([0])
        ax.set_xticklabels([model_names[m]], rotation=45)
        ax.set_title("Power Deficit")

        if m == 0:
            ax.set_ylabel("Power [PW]")

        # --- DOWN REGULATION subplot ---
        # (previously positive deviations)
        ax = axes[0, 2 * m + 1]
        bottom = 0
        for i in range(nP):
            ax.bar(0, dev_down[i], width=0.8, bottom=bottom, color=bar_colors[i])
            bottom += dev_down[i]

        ax.axhline(NI_down, color='black', linestyle='--', linewidth=1)
        if sys_imbalance > 0:
            ax.axhline(sys_imbalance, color='orange', linestyle='-', linewidth=1.5)

        ax.set_ylim(0, max(NI_up, NI_down)*1.1)
        ax.set_xticks([0])
        ax.set_xticklabels([model_names[m]], rotation=45)
        ax.set_title("Power Excess")

    # --- Legend (last subplot) ---
    polluter_handles = [Patch(facecolor=bar_colors[i], label=polluters[i]) for i in range(nP)]
    NI_handles = [
        Line2D([0], [0], color="black", linestyle="--", label=r"NI$\downarrow$"),
        Line2D([0], [0], color="green", linestyle="--", label=r"NI$\uparrow$")
    ]
    SI_handles = [Line2D([0], [0], color="orange", linewidth=1.5, label="SI")]

    axes[0, -1].legend(
        handles=polluter_handles + NI_handles + SI_handles,
        title="Legend",
        bbox_to_anchor=(1.05, 1),
        loc='upper left'
    )

    plt.tight_layout()
    plt.show()

# Profits for each fold increasingly
def profits_across_folds(means, stds, profits, B_bids_all_folds_dev_diff,
                         figsize = (10, 8), 
                         fontsize = 18):
    def get_total_sum_per_config(jl_matrix, num_rows=12, num_cols=10):
        total_sums = np.zeros(num_cols)
        for c in range(num_cols): 
            config_sum = 0
            for r in range(num_rows): 
                config_sum += np.sum(jl_matrix[r, c])
            total_sums[c] = config_sum
        return total_sums

    mean_on_in, mean_off_in, mean_on_out, mean_off_out = means
    std_on_in, std_off_in, std_on_out, std_off_out = stds
    summed_profits_pol_in, summed_profits_pol_out, summed_profits_flex_in, summed_profits_flex_out, summed_profits_reg_in, summed_profits_reg_out = profits

    total_pol_in = get_total_sum_per_config(summed_profits_pol_in, num_rows=12, num_cols=10)
    total_pol_out = get_total_sum_per_config(summed_profits_pol_out, num_rows=12, num_cols=10)
    total_flex_in = get_total_sum_per_config(summed_profits_flex_in, num_rows=8, num_cols=10)
    total_flex_out = get_total_sum_per_config(summed_profits_flex_out, num_rows=8, num_cols=10)
    total_reg_in = get_total_sum_per_config(summed_profits_reg_in, num_rows=3, num_cols=10)
    total_reg_out = get_total_sum_per_config(summed_profits_reg_out, num_rows=3, num_cols=10)
    total_in = total_pol_in + total_flex_in + total_reg_in
    total_out = total_pol_out + total_flex_out + total_reg_out
    diff_values = total_pol_in - total_pol_out

    # 2. Sort after mean onshore in-sample Pw
    sort_indices = np.argsort(mean_on_in)

    # Sort everything based on the sorted indices
    sorted_diff_on = np.array(mean_on_in)[sort_indices]
    sorted_diff_off = np.array(mean_off_in)[sort_indices]
    sorted_mean_out_on = np.array(mean_on_out)[sort_indices]
    sorted_mean_out_off = np.array(mean_off_out)[sort_indices]
    sorted_std_in = np.array(std_on_in)[sort_indices]
    sorted_std_out = np.array(std_on_out)[sort_indices]
    sorted_total_in = total_in[sort_indices]
    sorted_total_out = total_out[sort_indices]
    sorted_diffs = diff_values[sort_indices]

    # prepare x labels
    sorted_labels = [f"Fold {i+1}" for i,  on, std_in, on_out, std_out in zip(sort_indices, sorted_diff_on, sorted_std_in, sorted_mean_out_on, sorted_std_out)]
    x = np.arange(len(sorted_labels))

    # 3. (Dual Axis)
    fig, ax1 = plt.subplots(figsize=figsize)

    # --- Left axis (Axis 1): Total Profit ---
    ln1 = ax1.plot(x, sorted_total_in, 'o--', color='#08306b', label='Total In-Sample Profit', linewidth=2.5)
    ln2 = ax1.plot(x, sorted_total_out, 's-', color='#EE7518', label='Total Out-Sample Profit', linewidth=2.5)

    ax1.set_ylabel("Total profit [€]", color='black', fontsize = fontsize)
    ax1.tick_params(axis='y', labelcolor='black')
    ax1.set_xticks(x)
    ax1.set_xticklabels(sorted_labels, rotation=45)

    # --- Right axis (Axis 2): Sum of FCB Balancing Bids as bars ---
    ax2 = ax1.twinx()

    # Compute sum of FCB balancing bids for each fold (sorted)
    sum_fcb_bids = np.array([np.sum(B_bids_all_folds_dev_diff[f'Fold_{i+1}']['FCB']) for i in sort_indices])

    bars = ax2.bar(x, sum_fcb_bids, alpha=0.15, color='gray', label='Difference in deviations', width=0.7)
    ax2.set_ylabel(rf"$\Delta$ Expected deviations (In - Out) [MWh]", color='gray', fontsize = fontsize)
    ax2.tick_params(axis='y', labelcolor='gray')

    ax1.grid(True, linestyle=':', alpha=0.5)

    # Combine legends: profit lines + FCB bars
    lines = ln1 + ln2 + [bars]
    labels_legend = [l.get_label() for l in lines]
    ax1.legend(lines, labels_legend, loc='upper left')

    # Add horizontal arrow on ax1
    ax1.annotate(
        '',  # no text, just arrow
        xy=(9, 39650),  # arrow tip (right end)
        xytext=(1, 39650),    # arrow start (left end)
        arrowprops=dict(facecolor='black', arrowstyle='-|>', lw=2)
    )

    # Add text label
    ax1.text(
        len(x)/2,  # roughly center
        39660,  # slightly below the arrow
        rf"Increasing Expected $P_{{g,\omega}}^{{RT}}$",
        ha='center', va='bottom', fontsize=fontsize
    )

    return fig

# Fct used in systems_costs()
def split_costs(system_costs):
    pol  = [np.sum(c["pol"])  for c in system_costs]
    flex = [np.sum(c["flex"]) for c in system_costs]
    reg  = [np.sum(c["reg"])  for c in system_costs]
    nP = len(system_costs[0]['pol'])
    tot = [np.sum(c["reg"]) + np.sum(c["flex"]) + np.sum(c["pol"]) for c in system_costs]
    procurement_cost_up = [np.sum([c['pol'][p][2] for p in range(nP)]) for c in system_costs]
    procurement_cost_down = [np.sum([c['pol'][p][3] for p in range(nP)]) for c in system_costs]

    return np.array(pol), np.array(flex), np.array(reg), np.array(tot), np.array(procurement_cost_up), np.array(procurement_cost_down)

# Plotting the systems cost
def systems_costs(system_costs_BM_in, system_costs_BM_out, system_costs_FCB_in, system_costs_FCB_out, figsize = (10, 8), fontsize = 18):
    BM_in  = split_costs(system_costs_BM_in)
    BM_out = split_costs(system_costs_BM_out)
    FCB_in  = split_costs(system_costs_FCB_in)
    FCB_out = split_costs(system_costs_FCB_out)

    components = ["Deviators", "Flexible", "Base Loads", "Total"]
    n_folds = len(BM_in[0])  # assuming BM_in[i] is length n_folds

    x = np.arange(len(components))

    fig = plt.figure(figsize=figsize)

    for i, comp in enumerate(components):
        # In-sample
        y_in = (np.array(FCB_in[i]) - np.array(BM_in[i]))/np.array(BM_in[i])*100
        plt.scatter(np.full(n_folds, x[i]-0.1), 
                    y_in, color='#08306b', alpha=0.6, label='In-sample' if i==0 else "")
        # Out-of-sample
        y_out = (np.array(FCB_out[i]) - np.array(BM_out[i]))/np.array(BM_out[i])*100
        plt.scatter(np.full(n_folds, x[i]+0.1), 
                    y_out, color='#EE7518', alpha=0.6, label='Out-of-sample' if i==0 else "")
    
    # Optional: mean markers
    for i in range(4):
        plt.plot(x[i]-0.1, np.mean((np.array(FCB_in[i]) - np.array(BM_in[i]))/np.array(BM_in[i]))*100, 'o', color='#08306b', markersize=10)
        plt.plot(x[i]+0.1, np.mean((np.array(FCB_out[i]) - np.array(BM_out[i]))/np.array(BM_out[i]))*100, 'o', color='#EE7518', markersize=10)
    
    print("total diff in : ", np.mean(np.array(FCB_in[0]) - np.array(BM_in[0])))
    print("total diff out: ", np.mean(np.array(FCB_out[0]) - np.array(BM_out[0])))

    penalty_in = np.mean((np.array(FCB_in[0]) - np.array(BM_in[0]))/np.array(BM_in[0]))*100
    penalty_out = np.mean((np.array(FCB_out[0]) - np.array(BM_out[0]))/np.array(BM_out[0]))*100
    penalty_up = penalty_out - penalty_in

    penalty_down_in = np.mean(np.array(FCB_in[5]) - np.array(BM_in[5]))
    penalty_down_out = np.mean(np.array(FCB_out[5]) - np.array(BM_out[5]))

    penalty_up_2 = np.mean(np.array(FCB_out[4]) - np.array(BM_out[4]))
                           
    print("penalty down in insample", penalty_down_in)
    print("penalty down in outsample", penalty_down_out)
    print("penalty up in outsample", penalty_up_2)
    print("total diff penalty: ", penalty_down_out + penalty_up_2)
    print("systems costs difference in outsample in %", np.mean((np.array(FCB_out[3]) - np.array(BM_out[3]))/np.array(BM_out[i]))*100)
    print("TOtal systems cost in FCB out", np.mean(np.array(FCB_out[3])))

    plt.annotate(
    '',  # no text, just arrow
    xy=(0-0.1, penalty_in),  # arrow tip (right end)
    xytext=(0-0.1, 0),    # arrow start (left end)
    arrowprops=dict(color='gray', arrowstyle='-|>', lw=1)
    )

    # Add text label
    plt.text(
        0,  # roughly center
        penalty_in/2,  # slightly below the arrow
        f"Due to: procurement \ncost for down-regulation",
        ha='left', va='center', fontsize=fontsize-2, color = 'gray'
    )


    plt.annotate(
        '',  # no text, just arrow
        xy=(0 + 0.1, penalty_in + penalty_up),  # arrow tip (right end)
        xytext=(0 + 0.1, penalty_in),    # arrow start (left end)
        arrowprops=dict(color='gray', arrowstyle='-|>', lw=1)
    )

    # Add text label
    plt.text(
        0.2,  # roughly center
        penalty_in + (penalty_up)/2,  # slightly below the arrow
        f"Due to: procurement \ncost for up-regulation",
        ha='left', va='center', fontsize=fontsize-2, color = 'gray'
    )

    plt.xticks(x, components, rotation = 45)
    plt.ylabel(rf"$\Delta$ System Costs (FCB − BM) [%]", fontsize = fontsize)
    plt.axhline(0, color='black', linestyle='--', lw=1)
    plt.legend(fontsize = fontsize, loc='center right')
    plt.grid(axis='y', linestyle='--', alpha=0.5)
    plt.margins(0.2, 0.2)
    return fig

#plots power imbalances and difference between two models
def plot_power_deviations_models(
    model_names,
    deviations_list,
    polluters,
    base_color="Blues",
    figsize=(9.0, 7.2),

    # geometry
    bar_width=0.55,
    x_step=0.80,

    # fonts
    title_fontsize=16,
    label_fontsize=14,
    tick_fontsize=12,
    legend_fontsize=12,
    legend_title_fontsize=13,

    # styling
    grid_alpha=0.25,
    separator_color="0.25",
    separator_lw=0.6,
    legend_title="Generators",
    arrow_tol=1e-2,   # minimum difference to show arrow
    show_arrows=True, # NEW: control arrow display
):
    """
    Two-model stacked deviations + stacked difference (FCB − benchmark).
    Arrows show significant contributors if show_arrows=True.
    Y-limits are adaptive and consistent across subplots.
    """

    if len(model_names) != 2:
        raise ValueError("Exactly two models are required.")

    deviations_list = [np.asarray(d, float).ravel() for d in deviations_list]
    nP = len(polluters)

    for d in deviations_list:
        if d.size != nP:
            raise ValueError("Deviation vectors must match number of polluters.")

    diff_deviations = deviations_list[1] - deviations_list[0]
    polluters = list(polluters)

    # ---------------- Colors ----------------
    name_order = np.argsort([p.lower() for p in polluters])
    polluters_sorted = [polluters[i] for i in name_order]

    cmap = plt.get_cmap(base_color)
    xs = np.linspace(0.05, 0.95, nP)
    colors_sorted = [cmap(x) for x in xs]
    color_by_name = {name: colors_sorted[i] for i, name in enumerate(polluters_sorted)}

    # ---------------- Figure ----------------
    fig, (ax, ax_diff) = plt.subplots(
        1, 2,
        figsize=figsize,
        sharey=False,
        gridspec_kw={"width_ratios": [2, 1]}
    )

    for a in (ax, ax_diff):
        a.set_facecolor("white")
        a.set_frame_on(True)
        for spine in a.spines.values():
            spine.set_visible(True)
            spine.set_linewidth(0.8)
            spine.set_color("0.3")
    fig.patch.set_facecolor("white")

    # ---------------- Stacked bar helper ----------------
    def stacked_bar(ax_, xpos, values, half=1):
        pos, neg = 0.0, 0.0
        bottoms = []
        for val, name in zip(values, polluters):
            bottom = pos if val >= 0 else neg
            ax_.bar(
                xpos,
                val,
                bar_width*half,
                bottom=bottom,
                color=color_by_name[name],
                edgecolor=separator_color,
                linewidth=separator_lw,
            )
            bottoms.append(bottom)
            if val >= 0:
                pos += val
            else:
                neg += val
        return bottoms, pos, neg  # store bottoms for arrows

    # ---------------- Left: model deviations ----------------
    x_positions = np.arange(2) * x_step
    bottoms0, pos0, neg0 = stacked_bar(ax, x_positions[0], deviations_list[0])
    bottoms1, pos1, neg1 = stacked_bar(ax, x_positions[1], deviations_list[1])

    ax.set_xticks(x_positions)
    ax.set_xticklabels(model_names, fontsize=tick_fontsize)
    ax.set_ylabel(r"Power Deviations $\Delta p_{{g, \omega}}$ [MW]", fontsize=label_fontsize)

    # ---------------- Determine y-limits ----------------
    ymin_left = min(0.0, neg0, neg1)
    ymax_left = max(pos0, pos1)
    buffer = 0.05 * max(abs(ymax_left), abs(ymin_left), 1)
    ax.set_ylim(ymin_left - buffer, ymax_left + buffer)

    # ---------------- Right: difference ----------------
    bottoms_diff, pos_diff, neg_diff = stacked_bar(ax_diff, 0.0, diff_deviations, half=1.6)
    ax_diff.set_xlim(-bar_width, bar_width)
    ax_diff.set_xticks([0.0])
    ax_diff.set_xticklabels([f"{model_names[1]} − {model_names[0]}"], fontsize=tick_fontsize)

    # Scale right subplot to match left subplot
    ylim_max = max(ax.get_ylim()[1], pos_diff)
    ylim_min = min(ax.get_ylim()[0], neg_diff)
    ax_diff.set_ylim(ylim_min - buffer, ylim_max + buffer)

    # ---------------- Styling ----------------
    for a in (ax, ax_diff):
        a.yaxis.grid(True, alpha=grid_alpha)
        a.set_axisbelow(True)
        a.tick_params(axis="y", labelsize=tick_fontsize)
        a.tick_params(axis="x", labelsize=tick_fontsize)

    # ---------------- Legend ----------------
    sorted_idx = np.argsort([p.lower() for p in polluters])
    handles = [Patch(facecolor=color_by_name[polluters[i]], edgecolor=separator_color, label=polluters[i])
               for i in sorted_idx]

    fig.legend(
        handles=handles,
        title=legend_title,
        loc="center left",
        bbox_to_anchor=(1.02, 0.5),
        frameon=False,
        fontsize=legend_fontsize,
        title_fontsize=legend_title_fontsize,
        handlelength=1.2,
        labelspacing=0.8,
    )

    # ---------------- Arrows for main contributors ----------------
    if show_arrows:
        for i, name in enumerate(polluters):
            diff_val = diff_deviations[i]
            if abs(diff_val) < arrow_tol:
                continue  # skip negligible differences

            # arrow from left model top to right model top
            y_start = deviations_list[0][i]/2 + bottoms0[i]
            y_end = deviations_list[1][i]/2 + bottoms1[i]

            ax.annotate(
                "",
                xy=(x_positions[1], y_end),
                xytext=(x_positions[0], y_start),
                arrowprops=dict(
                    arrowstyle="->",
                    color="black",
                    lw=1.2,
                ),
            )

            # add label near left model bar
            ax.text(
                x_positions[0]-0.05,
                y_start + 2,
                name,
                fontsize=legend_fontsize-2,
                ha="center",
                va="bottom" if diff_val > 0 else "top"
            )

    plt.tight_layout()
    plt.show()
    return fig

#plots the stacked DA bids by technology
def plot_DA_bids_transposed_models_by_type_generators(
    DA_bids,
    types_to_plot=None,
    names_by_type=None,
    tech_colors=None,
    bar_height=0.38,
    spacing=0.08,
    sep_color="white",
    sep_lw=1.0,
    model_edgecolor="black",
    model_linewidth=1.1,
    sort_types_desc=True,
    sort_tech_desc=True,
    sort_units_desc=True,
    label_models_at_end=True,
    label_x_pad=0.01,
    # NEW: sizing + legend placement
    base_fontsize=18,        # overall font size
    tick_fontsize=18,
    label_fontsize=18,
    title_fontsize=18,
    figsize=(10, 8),
    legend_fontsize=17,
    legend_title_fontsize=17,
    legend_loc="upper right",  # inside the plot
):
    # Apply bigger fonts for this plot only
    plt.rcParams.update({
        "font.size": base_fontsize,
        "axes.titlesize": title_fontsize,
        "axes.labelsize": label_fontsize,
        "xtick.labelsize": tick_fontsize,
        "ytick.labelsize": tick_fontsize,
        "legend.fontsize": legend_fontsize,
    })

    def _tech_of(name):
        try:
            return tech_label(name)  # noqa: F821
        except NameError:
            import re
            s = str(name).strip()
            return re.sub(r"([_\-\s]*)\d+$", "", s)

    def _type_total(gtype, model_names):
        unit_names = list(names_by_type.get(gtype, []))
        tot = 0.0
        for model in model_names:
            vals = np.asarray(DA_bids[model].get(gtype, [])).reshape(-1)
            if len(vals) == len(unit_names):
                tot += float(np.sum(vals))
        return tot

    if names_by_type is None:
        raise ValueError("names_by_type is required (dict: type -> list of unit names).")

    model_names = list(DA_bids.keys())
    M = len(model_names)
    if M == 0:
        raise ValueError("DA_bids is empty (no models).")

    if tech_colors is None:
        tech_colors = {}

    all_types = sorted({gtype for model in DA_bids.values() for gtype in model.keys()})
    if types_to_plot is None:
        types_to_plot = all_types
    else:
        types_to_plot = [t for t in types_to_plot if t in all_types]

    if sort_types_desc:
        types_to_plot = sorted(types_to_plot, key=lambda t: _type_total(t, model_names), reverse=True)

    offsets = -np.arange(M) * (bar_height + spacing)
    T = len(types_to_plot)
    y_centers = np.arange(T) * (M * (bar_height + spacing) + 0.35)

    fig, ax = plt.subplots(figsize=figsize)

    xmax = 0.0
    techs_present, seen_techs = [], set()

    for tidx, gtype in enumerate(types_to_plot):
        if gtype not in names_by_type:
            raise ValueError(f"names_by_type missing key: {gtype}")

        unit_names_base = list(names_by_type[gtype])
        name_to_tech = {n: _tech_of(n) for n in unit_names_base}

        for midx, model in enumerate(model_names):
            vals = np.asarray(DA_bids[model].get(gtype, [])).reshape(-1)
            if len(vals) != len(unit_names_base):
                raise ValueError(
                    f"Length mismatch: DA_bids['{model}']['{gtype}'] has {len(vals)} "
                    f"but names_by_type['{gtype}'] has {len(unit_names_base)}"
                )

            items = [(name, float(v), name_to_tech[name]) for name, v in zip(unit_names_base, vals) if float(v) > 0]

            tech_to_items = {}
            for name, v, tech in items:
                tech_to_items.setdefault(tech, []).append((name, v))

            tech_list = list(tech_to_items.keys())
            if sort_tech_desc:
                tech_list.sort(key=lambda t: sum(v for _, v in tech_to_items[t]), reverse=True)
            else:
                tech_list.sort()

            y = y_centers[tidx] + offsets[midx]
            left = 0.0

            for tech in tech_list:
                unit_list = tech_to_items[tech]
                if sort_units_desc:
                    unit_list = sorted(unit_list, key=lambda x: x[1], reverse=True)

                if tech not in seen_techs:
                    techs_present.append(tech)
                    seen_techs.add(tech)

                color = tech_colors.get(tech, "gray")

                for _, v in unit_list:
                    ax.barh(
                        y=y, width=v, left=left, height=bar_height,
                        color=color, edgecolor=sep_color, linewidth=sep_lw
                    )
                    left += v

            ax.barh(
                y=y, width=left, left=0.0, height=bar_height,
                color="none", edgecolor=model_edgecolor, linewidth=model_linewidth
            )

            xmax = max(xmax, left)

            if label_models_at_end and left > 0:
                ax.text(left, y, f" {model}", va="center", ha="left", fontsize=tick_fontsize, color = "#545454")

    if label_models_at_end:
        pad = xmax * float(label_x_pad)
        for txt in ax.texts:
            x, y = txt.get_position()
            txt.set_position((x + pad, y))

    group_mid = y_centers - (M - 1) * (bar_height + spacing) / 2
    ax.set_yticks(group_mid)
    ax.set_yticklabels([t.replace("_", " ").title() for t in types_to_plot])

    ax.set_xlabel("DA Bids [MW]", fontsize = base_fontsize)
    ax.set_xlim(0, xmax * 1.2 if xmax > 0 else 1)
    ax.grid(axis="x", alpha=0.15, linewidth=0.8)

    # Legend INSIDE upper-right
    tech_handles = [Patch(facecolor=tech_colors.get(t, "gray"), label=t) for t in techs_present]
    if tech_handles:
        leg = ax.legend(
            handles=tech_handles,
            title="Technology",
            loc=legend_loc,          # inside
            frameon=False,
            fontsize=legend_title_fontsize,
        )
        leg.get_frame().set_alpha(0.9)

    plt.tight_layout()
    return fig

#Used in the DA bids transposed fct
def tech_label(name: str) -> str: 
    # strip trailing digits: "Onshore3" -> "Onshore", "CCGT2" -> "CCGT" 
    return re.sub(r"\d+$", "", name)

#pltots the difference in DA bids between two models, using a single-hue blue colormap
def plot_DA_bid_difference_one_color_from_nested_dict(
    DA_bids,
    generators_by_type,
    model_names=("Benchmark", "FCB"),
    types_to_plot=("Deviators", "Flexibles", "Base loads"),
    tech_colors=None,        # Excel-ish blue (change if you want darker)
    tol=1e-6,
    sort_by="abs",              # "abs" or "diff"
    figsize=(10.0, 8.0),
    title="Difference in Day-Ahead Bids (FCB − Benchmark)",
    xlabel="",
    gap_between_types=0.6,
    show_type_headers=True,
    fontsize = 20,
    title_fontsize=18,
    grid_color="#D9D9D9",
    grid_lw=1.1,
    zero_line_color="#A6A6A6",
    #zero_line_lw=1.2,
    group_box=True,             # draw faint group boxes on the left like Excel
    asy = False,
):
    """
    Excel-like horizontal bar chart of per-generator DA bid changes across generator types.

    diff = DA_bids[model2][type] - DA_bids[model1][type]
    """
    # ---- validate ----
    if len(model_names) != 2:
        raise ValueError("model_names must contain exactly 2 model names (model1, model2).")
    m1, m2 = model_names
    if m1 not in DA_bids or m2 not in DA_bids:
        raise ValueError(f"Both models {model_names} must exist in DA_bids keys: {list(DA_bids.keys())}")

    # ---- collect rows ----
    rows = []  # (type, generator_name, diff)
    for gtype in types_to_plot:
        if gtype not in generators_by_type:
            raise ValueError(f"generators_by_type missing key: {gtype}")
        if gtype not in DA_bids[m1] or gtype not in DA_bids[m2]:
            continue

        b1 = np.asarray(DA_bids[m1][gtype], dtype=float).ravel()
        b2 = np.asarray(DA_bids[m2][gtype], dtype=float).ravel()
        gens = list(generators_by_type[gtype])
        print(gens)
        if b1.size != len(gens) or b2.size != len(gens):
            raise ValueError(
                f"Size mismatch for type '{gtype}': "
                f"{m1} has {b1.size}, {m2} has {b2.size}, generators list has {len(gens)}"
            )

        diff = b2 - b1
        for name, d in zip(gens, diff):
            d = float(d)
            if abs(d) > tol:
                if asy == True and gtype == 'Deviators':
                    m = re.search(r'(On|Off)\d+.*?(\d+\.\d+)', name)
                    tech = f"{m.group(1)}_{m.group(2)}"
                else:
                    tech = re.sub(r'\d+$', '', name)
                color = tech_colors[tech]
                rows.append((gtype, str(name), d, color))
    if not rows:
        raise ValueError("No non-zero differences found (after tol).")

    # ---- sort within each type ----
    ordered_rows = []
    for gtype in types_to_plot:
        sub = [r for r in rows if r[0] == gtype]
        if not sub:
            continue
        if sort_by == "diff":
            sub = sorted(sub, key=lambda r: r[2], reverse=True)
        else:
            sub = sorted(sub, key=lambda r: abs(r[2]), reverse=True)
        ordered_rows.extend(sub)

    # ---- y positions with gaps between types ----
    y_positions, y_labels, diffs, bar_color = [], [], [], []
    # store start/end y (in data coords) for each group for separators/labels
    group_spans = []  # (gtype, y_start, y_end)

    y = 0.0
    last_type = None
    y_start = None
    
    for gtype, name, d, color in ordered_rows:
        if last_type is None or gtype != last_type:
            # close previous span
            if last_type is not None:
                group_spans.append((last_type, y_start, y - 1.0))
                y += gap_between_types
            y_start = y
            last_type = gtype
        
        bar_color.append(color)
        y_positions.append(y)
        y_labels.append(name)
        diffs.append(d)
        y += 1.0

    # close last span
    group_spans.append((last_type, y_start, y - 1.0))

    diffs = np.asarray(diffs, dtype=float)
    vmax = float(np.max(np.abs(diffs)))
    vmax_pos = float(np.max(diffs))
    vmax_neg = float(np.min(diffs))

    # ---- plot ----
    fig, ax = plt.subplots(figsize=figsize)
    fig.patch.set_facecolor("white")
    ax.set_facecolor("white")

    # bars
    ax.barh(y_positions, diffs, color=bar_color, edgecolor="none", height=0.62)

    # Add values on top of bars
    for y, d in zip(y_positions, diffs):
        ax.text(d + d*0.04, y, f"{d:+.1f}",
                va='center', ha='right' if d < 0 else 'left', color='black', fontsize=fontsize)

    ax.axvline(0, color=zero_line_color, zorder=0)

    # ticks/labels
    ax.set_yticks(y_positions)
    ax.set_yticklabels(y_labels, fontsize=fontsize)
    ax.invert_yaxis()

    if not xlabel:
        xlabel = rf"$\Delta$DA bids [MW]"
    ax.set_xlabel(xlabel, fontsize=fontsize)

    ax.set_title(title, fontsize=title_fontsize, pad=18, color="black")

    maxval = max(abs(vmax_neg), vmax_pos)
    # x limits and Excel-like vertical grid
    ax.set_xlim(-maxval * 1.25, maxval * 1.3)
    ax.xaxis.grid(True, color = grid_color, linewidth=grid_lw)
    ax.yaxis.grid(False)
    ax.set_axisbelow(True)
    ax.set_frame_on(True)                     # ensure frame is on
    ax.patch.set_edgecolor('black')           # frame color

    ax.tick_params(axis="x", labelsize=fontsize, colors="black")
    ax.tick_params(axis="y", colors="black", length=0)

    # add group labels in the left gutter using axes coords (robust)
    if show_type_headers:
        # widen left margin for the group label column
        fig.subplots_adjust(left=0.22, right=0.86, top=0.88, bottom=0.15)

        for (gtype, y0, y1) in group_spans:
            y_mid = (y0 + y1) / 2.0

            # convert y_mid (data) -> axes fraction for stable placement
            y_axes = ax.transAxes.inverted().transform(ax.transData.transform((0, y_mid)))[1]

            # group label
            # ax.text(
            #     -0.22, y_axes+0.05, gtype,
            #     transform=ax.transAxes,
            #     rotation=90,
            #     va="center", ha="center",
            #     fontsize=tick_fontsize,
            #     color="black",
            # )

            # faint group separator line (between groups)
            ax.hlines(y1 + 0.5, *ax.get_xlim(), colors="#E6E6E6", linewidth=1)

    else:
        fig.subplots_adjust(right=0.95, top=0.88, bottom=0.15)

    if asy :
        good_y = [y for y, lbl in zip(y_positions, y_labels)
                if ("0.05" in lbl and "Deviators" in [r[0] for r in rows])]
        bad_y  = [y for y, lbl in zip(y_positions, y_labels)
                if ("0.25" in lbl and "Deviators" in [r[0] for r in rows])]

        y_good_start, y_good_end = min(good_y), max(good_y)
        y_bad_start,  y_bad_end  = min(bad_y),  max(bad_y)

        x_brace =  max(diffs)*0.165
        pad = 0.4
        add_vertical_brace(ax, y_good_start - pad, y_good_end + pad , x_brace, "Good \nforecasts")
        add_vertical_brace(ax, y_bad_start - pad,  y_bad_end + pad,  x_brace, "Bad \nforecasts")

    #{g: {name: d for (gt, name, d) in rows if gt == g} for g in types_to_plot}
    return fig

# Used in DA bid differnce 
def add_vertical_brace(ax, y0, y1, x, text, width=0.02, text_offset=8):
    ax.annotate(
        text,
        xy=(x*4.5 + width, (y0 + y1) / 2),
        xytext=(text_offset, 0),
        textcoords="offset points",
        va="center", ha="left",
        fontsize=15,
    )
    ax.plot(
        [x, x + width, x + width, x],
        [y0, y0, y1, y1],
        transform=ax.get_yaxis_transform(),
        clip_on=False,
        linewidth=0.6,
        color="black",
    )

def make_palette_from_anchors(hex_anchors, n):
    """
    Build an n-color palette by interpolating between anchor colors.
    hex_anchors: list like ["#EBA067", "#EE7518", "#D05C03"]
    """
    rgb = np.array([mcolors.to_rgb(h) for h in hex_anchors], dtype=float)

    # positions of anchors in [0,1]
    xp = np.linspace(0, 1, len(hex_anchors))
    x = np.linspace(0, 1, n)

    # interpolate each channel
    out = np.column_stack([
        np.interp(x, xp, rgb[:, 0]),
        np.interp(x, xp, rgb[:, 1]),
        np.interp(x, xp, rgb[:, 2]),
    ])
    return [mcolors.to_hex(c) for c in out]


#plot to show the penalty recovered 
def plot_socialized_vs_individualized_FCB_only(
    pBIN_fcb,                 
    Q_up, Q_down,             
    polluters,                
    tech_colors=None,         # if provided and asy=True
    asy=False,                
    social_color="gray",
    prev_social_color="#f7f7f7",
    show_prev_all_socialized=True,
    figsize=(10, 5.2),
    anchors=("#EBA067", "#EE7518", "#D05C03"),
    color_order="as_is",      
    fontsize = 18
):
    pBIN = np.asarray(pBIN_fcb, dtype=float)
    nP = pBIN.shape[0]
    if len(polluters) != nP:
        raise ValueError("polluters length must match pBIN rows (nP).")
    if pBIN.shape[1] < 4:
        raise ValueError("pBIN_fcb must have at least 4 columns (needs indices 2 and 3).")

    polluters = list(map(str, polluters))

    pen_up = np.abs(pBIN[:, 2])
    pen_down = np.abs(pBIN[:, 3])

    # --- define bar colors ---
    if asy and tech_colors is not None:
        colors_for_polluters = []
        for name in polluters:
            m = re.search(r'(On|Off)\d+.*?(\d+\.\d+)', name)
            if m:
                tech = f"{m.group(1)}_{m.group(2)}"  # e.g., "On_0.05"
            else:
                tech = name  # fallback
            colors_for_polluters.append(tech_colors.get(tech, "#999999"))

        legend_names = polluters
        legend_colors = colors_for_polluters
        
    else:
        bar_colors = make_palette_from_anchors(list(anchors), nP)
        if color_order == "sorted":
            poll_sorted = sorted(polluters, key=lambda s: s.lower())
            color_by_name = {name: bar_colors[i] for i, name in enumerate(poll_sorted)}
            colors_for_polluters = [color_by_name[name] for name in polluters]
            legend_names = poll_sorted
            legend_colors = [color_by_name[name] for name in poll_sorted]
        else:
            colors_for_polluters = bar_colors
            legend_names = polluters
            legend_colors = bar_colors

    fig, axes = plt.subplots(1, 2, figsize=figsize)
    ymax = max(Q_up, Q_down) * 1.20 if max(Q_up, Q_down) > 0 else 1.0

    def _stack_one(ax, values, Q, title_txt, q_label, fontsize):
        if show_prev_all_socialized:
            ax.bar(0, Q, width=1.0, color=prev_social_color, edgecolor="none", zorder=0)

        bottom = 0.0
        for i in range(nP):
            v = float(values[i])
            if v <= 0:
                continue
            ax.bar(
                0, v, width=0.55, bottom=bottom,
                color=colors_for_polluters[i],
                edgecolor="white", linewidth=0.6, zorder=2
            )
            bottom += v

        social = max(Q - bottom, 0.0)
        if social > 0:
            ax.bar(
                0, social, width=0.55, bottom=bottom,
                color=social_color, edgecolor="white", linewidth=0.6, zorder=2
            )
            bottom += social

        # --- Add values on top of bars ---
        ax.axhline(Q, color="black", linestyle="--", linewidth=1.2)
        ax.text(0, Q * 1.02, f"{q_label} = {Q:.2f} €", color="black", fontsize=fontsize, ha = 'center')

        # --- Add total % paid by polluters ---
        polluter_total = sum(values)
        pct_paid = polluter_total / Q * 100 if Q > 0 else 0
        ax.text(
            0, polluter_total,  # slightly above the top of stacked bar
            f"{pct_paid:.1f}% recovered",
            ha="center", va="bottom",
            fontsize=fontsize,
            color="white",
        )

        ax.set_title(title_txt, fontsize = fontsize)
        ax.set_xticks([])
        ax.tick_params(axis="y", labelsize=fontsize)
        #ax.set_xticklabels([""])
        ax.set_ylim(0, ymax)
        ax.grid(axis="y", alpha=0.2)

        return social

    axes[0].set_ylabel("Procurement cost [€]", fontsize = fontsize+2)
    _stack_one(axes[0], pen_up, Q_up, "Up-regulation", "Q↑", fontsize=fontsize)
    _stack_one(axes[1], pen_down, Q_down, "Down-regulation", "Q↓", fontsize=fontsize)

    # --- Legend ---
    polluter_handles = [Patch(facecolor=legend_colors[i], label=legend_names[i]) for i in range(nP)]
    social_handle = [Patch(facecolor=social_color, label="Socialized")]
    prev_handle = [Patch(facecolor=prev_social_color, label="Previously:\n100% socialized")]
    line_handles = [Line2D([0], [0], color="black", linestyle="--", label="Q↑ / Q↓")]

    axes[1].legend(
        handles=polluter_handles + social_handle + (prev_handle if show_prev_all_socialized else []) + line_handles,
        bbox_to_anchor=(1.02, 1.01),
        loc="upper left",
        frameon=True,
        fontsize=fontsize,
    )

    #plt.tight_layout()
    return fig

### plot that shows the profit vs. difference of Pw 
def plot_profit_vs_diff_on_off(pBIN_list, diff_on, diff_off, title="Summed profit vs (diff_on, diff_off)"):
    """
    pBIN_list: list of arrays, each array shape (nP, >=1) with profit in column 0
    diff_on, diff_off: arrays/list length = len(pBIN_list)
    """
    diff_on = np.asarray(diff_on, dtype=float).ravel()
    diff_off = np.asarray(diff_off, dtype=float).ravel()
    K = len(pBIN_list)

    if len(diff_on) != K or len(diff_off) != K:
        raise ValueError(f"diff_on/diff_off must have length {K}")

    profits = np.array([np.asarray(p)[:, 0].sum() for p in pBIN_list], dtype=float)

    fig, ax = plt.subplots(figsize=(7, 5))
    sc = ax.scatter(diff_on, diff_off, c=profits)  # color = profit

    cb = plt.colorbar(sc, ax=ax)
    cb.set_label("Summed profit (polluters)")

    ax.set_xlabel("diff_on")
    ax.set_ylabel("diff_off")
    ax.set_title(title)
    ax.grid(alpha=0.25)

    # optional: label points with fold index
    for i, (x, y) in enumerate(zip(diff_on, diff_off), start=1):
        ax.text(x, y, str(i), fontsize=9, ha="left", va="bottom")

    plt.tight_layout()
    plt.show()

    return profits


#Showing count of SI up/down
def plot_SI_in_out_BM_vs_BIN_one_axes(
    system_imbalance_BM,
    system_imbalance_BIN,
    model_labels=("BM", "FCB"),
    component_names=("Up", "Down", "None"),
    tol=0.0,                 # 0.0 = exactly zero; use 1e-6 to ignore tiny noise
    figsize=(10, 8),
    legend_anchor=(-0.1, -0.01), # as you requested
    legend_ncol=3,
    fontsize = 18
):
    """
    Two figures (OUT and IN). Each figure has one axes with up/down/none for BM and BIN.
    A specific series is NOT plotted if its whole array is (near-)zero.
    Legend is placed below using bbox_to_anchor=(0, -0.1).
    """

    def _to_mats(sysimb, name):
        in_mat = np.array([t[0] for t in sysimb], dtype=float)   # (K,3)
        out_mat = np.array([t[1] for t in sysimb], dtype=float)  # (K,3)
        if in_mat.ndim != 2 or out_mat.ndim != 2 or in_mat.shape[1] != 3 or out_mat.shape[1] != 3:
            raise ValueError(f"{name} must be list of (in_vec,out_vec) with vec length 3.")
        return in_mat, out_mat

    BM_in, BM_out = _to_mats(system_imbalance_BM, "system_imbalance_BM")
    BIN_in, BIN_out = _to_mats(system_imbalance_BIN, "system_imbalance_BIN")

    if BM_in.shape[0] != BIN_in.shape[0]:
        raise ValueError("BM and BIN must have same number of folds.")

    folds = np.arange(1, BM_in.shape[0] + 1)

    # styles: line style by model; marker by component
    model_style = {model_labels[0]: "-", model_labels[1]: "--"}
    comp_marker = {0: "o", 1: "s", 2: "^"}  # Up, Down, None
    cmap = {0: "#EE7518", 1: "#08306b", 2: "#2E8B57"}  # red, blue, gray
    def _series_is_all_zero(y):
        y = np.asarray(y, dtype=float)
        return np.all(np.abs(y) <= tol)

    def _plot_one(A, B, title):
        fig, ax = plt.subplots(figsize=figsize)

        for j in range(3):
            y_BM = A[:, j]
            y_BIN = B[:, j]

            # Plot BM component only if BM series not all-zero
            ax.plot(
                    folds, y_BM,
                    linestyle=model_style[model_labels[0]],
                    marker=comp_marker[j],
                    label=f"{model_labels[0]} – {component_names[j]}",
                    color =cmap[j],
                )

            # Plot BIN component only if BIN series not all-zero
            ax.plot(
                    folds, y_BIN,
                    linestyle=model_style[model_labels[1]],
                    marker=comp_marker[j],
                    label=f"{model_labels[1]} – {component_names[j]}",
                    color =cmap[j],
                )

        ax.axhline(0, color="0.35", linewidth=2)
        ax.set_xticks(folds)
        ax.set_ylim(bottom = 0)
        ax.set_xlabel("Fold", fontsize = fontsize)
        ax.set_ylabel("Count of scenarios by regulation type", fontsize = fontsize)
        ax.set_title(title)
        ax.grid(axis="y", alpha=0.5)
        ax.tick_params(labelsize=fontsize-1)

        # Legend below plot
        handles, labels = ax.get_legend_handles_labels()
        if handles:
            ax.legend(
                loc="upper left",
                bbox_to_anchor=legend_anchor,   # (0, -0.1)
                ncol=legend_ncol,
                frameon=True,
                fontsize=fontsize -1,
            )

        plt.tight_layout()
        return fig

    # _plot_one(BM_out, BIN_out, "System imbalance (Out-of-sample): BM vs FCB")
    # _plot_one(BM_in,  BIN_in,  "System imbalance (In-sample): BM vs FCB")
    fig_out = _plot_one(BM_out, BIN_out, "")
    fig_in = _plot_one(BM_in,  BIN_in,  "")
    return fig_out, fig_in


def plot_profit_in_out_vs_diff_on_off(
    diff_on, diff_off,
    profit_in, profit_out,
    title="Summed profit vs (diff_on, diff_off)",
    labels=("In-sample profit", "Out-of-sample profit"),
    figsize=(11, 4.8),
):
    diff_on = np.asarray(diff_on, dtype=float).ravel()
    diff_off = np.asarray(diff_off, dtype=float).ravel()
    profit_in = np.asarray(profit_in, dtype=float).ravel()
    profit_out = np.asarray(profit_out, dtype=float).ravel()

    if not (len(diff_on) == len(diff_off) == len(profit_in) == len(profit_out)):
        raise ValueError("diff_on, diff_off, profit_in, profit_out must all have the same length (n_folds).")

    n = len(diff_on)

    # shared limits so both panels are comparable
    xpad = 0.05 * (diff_on.max() - diff_on.min() if diff_on.max() != diff_on.min() else 1.0)
    ypad = 0.05 * (diff_off.max() - diff_off.min() if diff_off.max() != diff_off.min() else 1.0)
    xlim = (diff_on.min() - xpad, diff_on.max() + xpad)
    ylim = (diff_off.min() - ypad, diff_off.max() + ypad)

    fig, axes = plt.subplots(1, 2, figsize=figsize, sharex=True, sharey=True)

    # --- IN ---
    ax = axes[0]
    sc1 = ax.scatter(diff_on, diff_off, c=profit_in)
    cb1 = plt.colorbar(sc1, ax=ax)
    cb1.set_label(labels[0])
    ax.set_title("In-sample")
    ax.set_xlabel("diff_on")
    ax.set_ylabel("diff_off")
    ax.set_xlim(*xlim)
    ax.set_ylim(*ylim)
    ax.grid(alpha=0.25)

    # annotate fold numbers
    for i, (x, y) in enumerate(zip(diff_on, diff_off), start=1):
        ax.text(x, y, str(i), fontsize=9, ha="left", va="bottom")

    # --- OUT ---
    ax = axes[1]
    sc2 = ax.scatter(diff_on, diff_off, c=profit_out)
    cb2 = plt.colorbar(sc2, ax=ax)
    cb2.set_label(labels[1])
    ax.set_title("Out-of-sample")
    ax.set_xlabel("diff_on")
    ax.set_xlim(*xlim)
    ax.set_ylim(*ylim)
    ax.grid(alpha=0.25)

    for i, (x, y) in enumerate(zip(diff_on, diff_off), start=1):
        ax.text(x, y, str(i), fontsize=9, ha="left", va="bottom")

    fig.suptitle(title)
    plt.tight_layout()
    plt.show()

#penalty for deviators across folds, so far not used but maybe useful for basecase?
def plot_penalty_violins(in_data, out_data, Deviators, fontsize = 18, figsize = (10, 8), ylim = None, title_suffix="Regulation"):
    """
    in_data: 12x10 matrix (12 Deviators, 10 Polluters/Folds)
    out_data: 12x10 matrix (12 Deviators, 10 Polluters/Folds)
    """
    labels = ["In-Sample", "Out-Sample"]
    colors = ['#08306b', '#EE7518'] # Your specific Dark Blue and Crimson
    
    nD = 6*2 # 12 deviators divided into 2 for up and down regulation
    base = np.arange(0, nD)
    step = 0.3
    offsets = [-step/2, step/2]

    fig, ax = plt.subplots(figsize=figsize)
    legend_handles = []

    for idx, (matrix, label, color, offset) in enumerate(zip([in_data, out_data], labels, colors, offsets)):
        pos = base + offset
        # Each Deviator (row) has 10 values (columns) to form the violin shape
        formatted_data = [np.array(matrix)[i, :] for i in range(nD)]

        vp = ax.violinplot(
            formatted_data,
            positions=pos,
            widths=step * 0.9,
            showmeans=True,
            showextrema=True
        )

        # Apply Colors
        for b in vp['bodies']:
            b.set_alpha(0.5)
            b.set_facecolor(color)
            b.set_edgecolor(color)

        for part in ('cmeans', 'cmins', 'cmaxes', 'cbars'):
            if part in vp:
                vp[part].set_edgecolor(color)
                vp[part].set_linewidth(1.5)

        legend_handles.append(Line2D([0], [0], color=color, lw=6, label=label))

    # Formatting
    ax.set_xticks(base)
    ax.set_xticklabels(Deviators, rotation=45, ha='right', fontsize=fontsize)
    ax.set_xlabel("Deviators", fontsize=fontsize)
    ax.set_ylabel("Penalty Value", fontsize=fontsize)
    ax.set_title(f"Comparison of {title_suffix} Penalties", fontsize=fontsize)
    ax.legend(handles=legend_handles, loc="upper right")
    ax.grid(True, alpha=0.6)
    if ylim is not None: 
        ax.set_ylim(ylim)
    plt.tight_layout()
    plt.show()

#Plot penalty violins for asymmetry
def plot_penalty_violins_combined(
    up_in, up_out,
    down_in, down_out,
    Deviators,
    fontsize=18,
    figsize=(10, 6),
    ylim=None,
    gap=1.5
):
    labels = ["In-Sample", "Out-Sample"]
    colors = ['#08306b', '#EE7518']

    nD = len(Deviators)
    step = 0.3
    offsets = [-step/2, step/2]

    # Base positions
    base_up = np.arange(nD)
    base_down = base_up + nD + gap

    fig, ax = plt.subplots(figsize=figsize)
    legend_handles = []

    def plot_block(data_in, data_out, base_positions):
        for matrix, label, color, offset in zip(
            [data_in, data_out], labels, colors, offsets
        ):
            pos = base_positions + offset
            formatted_data = [np.array(matrix)[i, :] for i in range(nD)]

            vp = ax.violinplot(
                formatted_data,
                positions=pos,
                widths=step * 0.9,
                showmeans=True,
                showextrema=True
            )

            for b in vp['bodies']:
                b.set_alpha(0.5)
                b.set_facecolor(color)
                b.set_edgecolor(color)

            for part in ('cmeans', 'cmins', 'cmaxes', 'cbars'):
                if part in vp:
                    vp[part].set_edgecolor(color)
                    vp[part].set_linewidth(1.5)

            if not any(h.get_label() == label for h in legend_handles):
                legend_handles.append(Line2D([0], [0], color=color, lw=6, label=label))

    # Plot both blocks
    plot_block(up_in, up_out, base_up)
    plot_block(down_in, down_out, base_down)

    # Vertical separator
    sep_x = nD - 0.5 + gap / 2
    ax.axvline(sep_x, color='k', linestyle='--', alpha=0.6)

    # X ticks and labels
    xticks = np.concatenate([base_up, base_down])
    xticklabels = Deviators + Deviators
    ax.set_xticks(xticks)
    ax.set_xticklabels(xticklabels, rotation=45, ha='right', fontsize=fontsize-3)

    # Group labels
    ax.text(base_up.mean(), ylim[1] if ylim else ax.get_ylim()[1],
            "Up-Regulation", ha='center', va='bottom', fontsize=fontsize)
    ax.text(base_down.mean(), ylim[1] if ylim else ax.get_ylim()[1],
            "Down-Regulation", ha='center', va='bottom', fontsize=fontsize)

    # Formatting
    ax.set_ylabel("Procurement Cost [€]", fontsize=fontsize)
    ax.legend(handles=legend_handles, loc="lower left", fontsize = fontsize)
    ax.grid(True, alpha=0.6)

    if ylim is not None:
        ax.set_ylim(ylim)

    plt.tight_layout()
    return fig

#profits difference plotting function
def stacked_posneg_diff_matrix(bench, model, polluters, markets, participant_type="Polluters",
                               figsize=(10,5), 
                               labelsize=12, 
                               titlesize=14,
                               fontsize=12, 
                               loc ="lower left"):
    nP, nM = bench.shape
    x = np.arange(nP)
    width = 0.35

    fig, ax = plt.subplots(figsize=figsize)

    # colors for markets
    blue = [ "#aec7e8","#16669f", "#428ac6"]
 
    def plot_stack(mat, xpos, colors):
        bot_pos = np.zeros(nP)
        bot_neg = np.zeros(nP)

        for m in range(mat.shape[1]):
            vals = mat[:, m]

            pos = np.where(vals > 0, vals, 0)
            neg = np.where(vals < 0, vals, 0)

            # positive stack
            ax.bar(xpos, pos, width, bottom=bot_pos, color=colors[m])
            bot_pos += pos

            # negative stack
            ax.bar(xpos, neg, width, bottom=bot_neg, color=colors[m])
            bot_neg += neg

        return bot_pos + bot_neg  # net profit

    diff = model - bench
    # plot left (benchmark)
    net_diff = plot_stack(diff, x, blue)
 
    ax.scatter(x, net_diff, color="black", s=20, label="Net profit difference")

    ax.set_xticks(x)
    ax.set_xticklabels(polluters, rotation=45, ha='right', fontsize=fontsize)
    ax.set_ylim(-300, 300)
    ax.tick_params(axis='y', labelsize=fontsize)
    ax.set_ylabel(rf"$\Delta$ Profit [€]", fontsize=labelsize)
    #ax.set_title(f"Profit decomposition per {participant_type}", fontsize=titlesize)
    ax.axhline(0, color="black", linewidth=0.8)
    #ax.set_ylim(min(np.min(net_bench), np.min(net_model)) * 1.1, max(np.max(net_bench), np.max(net_model)) * 1.1)

    # --- Create legend manually ---
    legend_handles = []


    # Model colors
    for i, market in enumerate(markets):
        legend_handles.append(Patch(facecolor=blue[i], label=f"{market}"))

    legend_handles.append(Line2D([0], [0], marker='o', color='w', label='Net profit difference',
                              markerfacecolor='black', markersize=8))

    ax.legend(handles=legend_handles, loc=loc, frameon=True, fontsize=fontsize)
    plt.tight_layout()
    return fig


#power deviations for different models (Q sensitivity)
def plot_power_deviation_differences_multi_models(
    model_names,
    deviations_list,
    generators,
    benchmark_name=None,
    benchmark_idx=0,
    diff_mode="model_minus_benchmark",  # "model_minus_benchmark" or "benchmark_minus_model"
    base_cmap="Blues",
    legend_mode="nonzero",   # "positive", "nonzero", "all"
    tol=1e-6,
    sort_by="abs",           # "abs", "diff", or None (sorting per bar uses same generator order)
    figsize=(10, 5),
    title="Difference in Mean Deviator Bids (Out-of-sample)",
    ylabel="Δ mean balancing bids [PW]",
):
    # --- checks ---
    if len(model_names) < 2:
        raise ValueError("Provide at least 2 models (benchmark + others).")
    if len(model_names) != len(deviations_list):
        raise ValueError("model_names and deviations_list must have the same length.")

    nG = len(generators)
    devs = [np.asarray(d, dtype=float).ravel() for d in deviations_list]
    if any(d.size != nG for d in devs):
        raise ValueError("Each deviation vector must match number of generators.")

    # --- benchmark selection ---
    if benchmark_name is not None:
        if benchmark_name not in model_names:
            raise ValueError("benchmark_name not found in model_names.")
        b_idx = model_names.index(benchmark_name)
    else:
        b_idx = benchmark_idx

    if not (0 <= b_idx < len(model_names)):
        raise ValueError("benchmark index out of range.")

    bench = devs[b_idx]
    bench_name = model_names[b_idx]

    # models to plot
    plot_ids = [i for i in range(len(model_names)) if i != b_idx]
    if len(plot_ids) == 0:
        raise ValueError("Nothing to plot: only benchmark provided.")

    # --- compute diffs matrix: (n_models_to_plot, nG) ---
    diffs = []
    plot_names = []
    for i in plot_ids:
        if diff_mode == "model_minus_benchmark":
            diffs.append(devs[i] - bench)
            plot_names.append(f"{model_names[i]} − {bench_name}")
        elif diff_mode == "benchmark_minus_model":
            diffs.append(bench - devs[i])
            plot_names.append(f"{bench_name} − {model_names[i]}")
        else:
            raise ValueError("diff_mode must be 'model_minus_benchmark' or 'benchmark_minus_model'.")

    diffs = np.vstack(diffs)  # shape: (nBars, nG)

    # --- pick generator order once (global), so stacks are comparable across bars ---
    # Here: sort by aggregated absolute diff across bars (or aggregated diff)
    agg = np.sum(np.abs(diffs), axis=0) if sort_by == "abs" else np.sum(diffs, axis=0)
    if sort_by in ("abs", "diff"):
        order = np.argsort(agg)[::-1]
    else:
        order = np.arange(nG)

    diffs = diffs[:, order]
    gens = [generators[i] for i in order]

    # --- filter generators that are ~0 across ALL bars (optional but usually helpful) ---
    keep = [j for j in range(nG) if np.max(np.abs(diffs[:, j])) > tol]
    if len(keep) == 0:
        raise ValueError("All differences are ~0 after tol (across all models).")

    diffs = diffs[:, keep]
    gens = [gens[j] for j in keep]

    nBars, nF = diffs.shape

    # --- colors for generators (consistent across bars) ---
    cmap = plt.colormaps.get_cmap(base_cmap)
    if nF == 1:
        colors = [cmap(0.75)]
    else:
        colors = [cmap(0.45 + 0.4*i/(nF-1)) for i in range(nF)]

    # --- plot ---
    fig, ax = plt.subplots(figsize=figsize)

    width = 0.6
    x_positions = np.arange(nBars)

    # track global ylim
    global_pos = np.zeros(nBars)
    global_neg = np.zeros(nBars)

    legend_handles = []

    for j, (gname, col) in enumerate(zip(gens, colors)):
        vals = diffs[:, j]

        for b in range(nBars):
            v = vals[b]
            if abs(v) <= tol:
                continue

            if v >= 0:
                ax.bar(
                    x_positions[b], v, width,
                    bottom=global_pos[b],
                    color=col, edgecolor="#f2f2f2", linewidth=0.6
                )
                global_pos[b] += v
            else:
                ax.bar(
                    x_positions[b], v, width,
                    bottom=global_neg[b],
                    color=col, edgecolor="#f2f2f2", linewidth=0.6
                )
                global_neg[b] += v

        include = (
            (legend_mode == "all") or
            (legend_mode == "positive" and np.any(vals > tol)) or
            (legend_mode == "nonzero" and np.any(np.abs(vals) > tol))
        )
        if include:
            legend_handles.append(Patch(facecolor=col, label=gname))

    ax.axhline(0, linewidth=1, color="0.25")

    ymax = max(np.max(global_pos), 0) * 1.2 if np.max(np.abs(global_pos)) > tol else 1
    ymin = min(np.min(global_neg), 0) * 1.2 if np.max(np.abs(global_neg)) > tol else -1
    ax.set_ylim(ymin, ymax)

    ax.set_xticks(x_positions)
    ax.set_xticklabels(plot_names, rotation=45)
    ax.set_ylabel(ylabel)
    ax.set_title(title, fontsize=14)

    ax.legend(
        handles=legend_handles,
        title="Generators",
        fontsize=12,
        title_fontsize=13,
        handlelength=1.2,
        labelspacing=0.4,
        bbox_to_anchor=(1.05, 0.5),
        loc="center left",
        frameon=False,
    )

    plt.tight_layout()
    plt.show()

    return diffs, gens, plot_names


#Plotting the da and balancing bids for each fold for on and offshore
def plot_all_folds_two_deviators_FCB_DAplusB(
    pol_out_bids,
    deviator_names,                 # full list (len 12)
    onshore_idx=0,                  # index of Onshore1 in deviator_names
    offshore_idx=6,                 # index of Offshore1 in deviator_names
    model="FCB",
    folds=None,                     # None -> all folds in sorted order
    figsize=(14, 6),
    colors = ["#EE7518", "#1F77B4"],                # FCB color, another color
    market_alpha=(0.55, 1.0),       # DA alpha, B alpha
    ylim=None,
    title="FCB bids: DA + Balancing (stacked) across folds",
    DA = True,
    ):
    all_folds = pol_out_bids["DA_bids"].keys()
    if folds is None:
        folds = all_folds
    else:
        for f in folds:
            if f not in all_folds:
                raise ValueError(f"Unknown fold '{f}'. Available: {all_folds}")

    F = len(folds)
    x = np.arange(F)

    # collect series
    DA_on, B_on, DA_off, B_off = [], [], [], []
    for fold in folds:
        DA_vec = np.asarray(pol_out_bids["DA_bids"][fold][model], dtype=float).reshape(-1)
        B_vec  = np.asarray(pol_out_bids["B_bids"][fold][model],  dtype=float).reshape(-1)

        DA_on.append(DA_vec[onshore_idx])
        B_on.append(B_vec[onshore_idx])
        DA_off.append(DA_vec[offshore_idx])
        B_off.append(B_vec[offshore_idx])

    DA_on = np.array(DA_on); B_on = np.array(B_on)
    DA_off = np.array(DA_off); B_off = np.array(B_off)

    D_on_avg = np.mean(DA_on + B_on)
    D_off_avg = np.mean(DA_off + B_off)

    D_on_DA_avg = np.mean(DA_on)
    D_off_DA_avg = np.mean(DA_off)
    D_on_B_avg = np.mean(B_on)
    D_off_B_avg = np.mean(B_off)

    D_on_DA_std = np.std(DA_on)
    D_off_DA_std = np.std(DA_off)
    D_on_B_std = np.std(B_on)
    D_off_B_std = np.std(B_off)

    print(f"Onshore1 DA avg: {D_on_DA_avg:.2f} MW (std: {D_on_DA_std:.2f})")
    print(f"Onshore1 B  avg: {D_on_B_avg:.2f} MW (std: {D_on_B_std:.2f})")
    print(f"Offshore1 DA avg: {D_off_DA_avg:.2f} MW (std: {D_off_DA_std:.2f})")
    print(f"Offshore1 B  avg: {D_off_B_avg:.2f} MW (std: {D_off_B_std:.2f})")

    #Calculate diifference for all generators 
    

    fig, ax = plt.subplots(figsize=figsize)

    width = 0.35
    x_on  = x - width/2
    x_off = x + width/2

    # Onshore1 stacked bar
    if DA:
        ax.bar(x_on, DA_on, width=width, color=colors[0], alpha=market_alpha[0],
            edgecolor="white", linewidth=0.6, label="DA")
    ax.bar(x_on, B_on, width=width, bottom=DA_on, color=colors[0], alpha=market_alpha[1],
           edgecolor="white", linewidth=0.6, label="Balancing")
    
    #add horizontal lines for averages
    if DA:
        ax.axhline(D_on_avg, color=colors[0], linestyle='--', linewidth=1.5, label=f'Onshore1 Avg: {D_on_avg:.2f} MW')
        ax.axhline(D_off_avg, color=colors[1], linestyle='--', linewidth=1.5, label=f'Offshore1 Avg: {D_off_avg:.2f} MW')

    # Offshore1 stacked bar
    if DA:
        ax.bar(x_off, DA_off, width=width, color=colors[1], alpha=market_alpha[0],
            edgecolor="white", linewidth=0.6)
    ax.bar(x_off, B_off, width=width, bottom=DA_off, color=colors[1], alpha=market_alpha[1],
           edgecolor="white", linewidth=0.6)

    # x labels
    ax.set_xticks(x)
    ax.set_xticklabels(folds, rotation=0)
    ax.set_ylabel("Bid [MW]")
    ax.set_title(title)
    ax.grid(axis="y", alpha=0.2)

    if ylim is not None:
        ax.set_ylim(*ylim)

    # Legends: deviators + markets (alpha)
    dev_handles = []

    if DA:
        dev_handles.append(
            Patch(facecolor=colors[0], alpha=market_alpha[0],
                label="On - DA segment")
        )

    dev_handles.append(
        Patch(facecolor=colors[0], alpha=market_alpha[1],
            label="On - Balancing segment")
    )

    if DA:
        dev_handles.append(
            Patch(facecolor=colors[1], alpha=market_alpha[0],
                label="Off - DA segment")
        )

    dev_handles.append(
        Patch(facecolor=colors[1], alpha=market_alpha[1],
            label="Off - Balancing segment")
    )

    ax.legend(handles=dev_handles, title="Deviator", loc="lower right", frameon=True)

    plt.tight_layout()
    plt.show()

# Plot for merit order Balancing
def merit_order_price_balancing(
    DA_bids, Pw, nS, rf_up, rf_down, Cf, nF, nP,
    nb_w=0, tol=1e-6, figsize = (10, 6), fontsize = 18
):
    """
    Compute balancing prices and plot activated bids with SCGT/CCGT colors and style
    similar to reserve market plot.
    """
    prices = np.zeros(nS)
    activated_bids_up = np.zeros((nS, nF))
    activated_bids_down = np.zeros((nS, nF))
    ppB = np.zeros((nS, nP))
    ppB_up_contr = np.zeros((nS, nP))
    ppB_down_contr = np.zeros((nS, nP))

    # Flexible units
    flex_names = [f"SCGT{i}" for i in range(1, 5)] + [ f"CCGT{i}" for i in range(1,5)] # Example names for 5 flexibles
    name2idx = {name: i for i, name in enumerate(flex_names)}

    # --- Compute deviations for each scenario ---
    for w in range(nS):
        for p in range(nP):
            ppB[w, p] = Pw[w, p] - DA_bids[0][p]
        total_imbalance = np.sum(ppB[w, :])
        for p in range(nP):
            imb = ppB[w, p]
            if total_imbalance > tol and imb > 0:
                ppB_down_contr[w, p] = imb
            elif total_imbalance < -tol and imb < 0:
                ppB_up_contr[w, p] = imb

    # --- Plot scenario nb_w ---
    w = nb_w
    total_req = np.sum(ppB[w, :])

    bids_reserve = [(rf_up[i], rf_down[i], Cf[i], flex_names[i]) for i in range(nF)]

    if total_req > 0:
        selected_bids = [(x[1], x[2], x[3]) for x in bids_reserve]
        regulation_type = "down"
    elif total_req < 0:
        selected_bids = [(x[0], x[2], x[3]) for x in bids_reserve]
        regulation_type = "up"
    else:
        selected_bids = []
        regulation_type = "none"

    sorted_bids = sorted(selected_bids, key=lambda x: x[1])

    # --- Determine activated bids and balancing price ---
    cumulative = 0.0
    balancing_price = 0.0
    last_cumulative = 0.0
    for bid, cost, name in sorted_bids:
        last_cumulative = cumulative
        cumulative += bid
        f = name2idx[name]

        if round(cumulative, 3) < round(abs(total_req), 3):
            if total_req > 0:
                activated_bids_down[w, f] = bid
            else:
                activated_bids_up[w, f] = bid
        else:
            balancing_price = cost
            activated = abs(total_req) - last_cumulative
            if total_req > 0:
                activated_bids_down[w, f] = activated
            else:
                activated_bids_up[w, f] = activated
            break

    prices[w] = balancing_price

    # --- Prepare activated plotting list ---
    activated_array = activated_bids_down[w, :] if total_req > 0 else activated_bids_up[w, :]
    activated_plotting_bids = [
        (activated_array[name2idx[name]], cost, name)
        for bid, cost, name in sorted_bids
        if activated_array[name2idx[name]] > 0
    ]

    # --- Plot ---
    fig, ax = plt.subplots(figsize=figsize)
    ax.set_xlim(0, 140)
    ax.set_ylim(0, 100)
    ax.set_xlabel("Aggregated activated balancing (MW)", fontsize=fontsize)
    ax.set_ylabel("Cost (€/MWh)", fontsize=fontsize)
    ax.tick_params(labelsize=fontsize)

    cumulative = 0.0
    last_cost = None
    for bid, cost, name in activated_plotting_bids:
        tech = "SCGT" if "SCGT" in name else "CCGT"
        color = "#3f8fc4" if tech == "SCGT" else "#08306b"

        # vertical step line if cost repeats
        if last_cost is not None and cost == last_cost:
            ax.vlines(cumulative, cost-cost*0.05, cost+ cost*0.05, color=color, linewidth=2)

        # vertical line for step
        if last_cost is not None:
            ax.vlines(cumulative, last_cost, cost, color=color, linewidth=2)

        # horizontal block
        ax.hlines(cost, cumulative, cumulative+bid, color=color, linewidth=3)

        # index-only label
        ax.text(cumulative + bid/2, cost + cost*0.08, name,
                ha='center', color=color, fontsize=fontsize)
        
        cumulative += bid
        last_cost = cost

    # --- Legend ---
    #ax.plot([], [], color="#3f8fc4", linewidth=4, label="SCGT")
    #ax.plot([], [], color="#08306b", linewidth=4, label="CCGT")

    # Total activated and balancing price lines
    ax.vlines(abs(total_req), 0, balancing_price*1.5, color="#2E6F40", linestyle="--", label="Total Deviations", linewidth = 2)
    ax.axhline(balancing_price, color="black", linestyle="--", label="Balancing price", linewidth = 2)
    # Annotations
    price_label = rf"$\lambda_{{\omega = {nb_w +1}}}^{{\up*}}$" if total_req < 0 else rf"$\lambda_{{\omega = {nb_w+ 1}}}^{{\downarrow*}}$"
    ax.text(cumulative-cumulative*0.1, balancing_price*1.3, f"{int(cumulative)} MW", color="#388E3C", fontsize=fontsize)
    ax.text(-14, balancing_price, price_label, color="black", fontsize=fontsize)

    ax.legend(loc="lower right", fontsize=fontsize)
    ax.grid(True, alpha=0.5)
    plt.tight_layout()

    return fig

# Used in merit order DA plot
def wind_location(name):
        return "onshore" if name.startswith("On") else "offshore"

#Plot for merit order DA and balancing
def merit_order_price_DA(DA_bids, Cp, Cf, Cr, Deviators, Flexibles, Base_loads, title = None, fontsize = 18):
    """
    Plot merit order curve for day-ahead market and return market price.

    Parameters
    ----------
    DA_bids : list of lists
        [polluter_bids, flexible_bids, regular_bids], each a list of MW values
    Cp, Cf, Cr : lists
        Corresponding cost for each generator type
    """
    
    # --- Collect bids and costs ---
    bids = []
    # Polluters
    for i, b in enumerate(DA_bids[0]):
        bids.append((b, Cp[i], 'polluter', Deviators[i]))
    # Flexibles
    for i, b in enumerate(DA_bids[1]):
        bids.append((b, Cf[i], 'flexible', Flexibles[i]))
    # Regular
    for i, b in enumerate(DA_bids[2]):
        bids.append((b, Cr[i], 'regular', Base_loads[i]))

    # --- Sort by cost ascending ---
    sorted_bids = sorted(bids, key=lambda x: x[1])

    # --- Total demand ---
    total_demand = sum(DA_bids[0]) + sum(DA_bids[1]) + sum(DA_bids[2])

    # --- Colors ---
    colors = {'onshore': "#EBA067", 'offshore': "#D05C03", 'flexible': "#2257a5", 'regular': "#388E3C"}

    # --- Initialize plot ---
    fig, ax = plt.subplots(figsize=(10,8))
    ax.set_xlabel("Aggregated generation [MW]", fontsize=fontsize)
    ax.set_ylabel("Cost [€/MWh]", fontsize=fontsize)
    ax.set_title(title, fontsize=fontsize)
    plt.xticks(fontsize=fontsize, rotation=45)
    plt.yticks(fontsize=fontsize)
    
    # --- Draw each generator block ---
    cumulative = 0.0
    last_cost = None
    counts_by_cost = {}
    x_vals, y_vals = [], []

    for bid, cost, typ, name in sorted_bids:
        if typ == "polluter":
            loc = wind_location(name)
            color = colors[loc]
        else:
            color = colors[typ]
            
        if last_cost is not None and cost == last_cost and cumulative > 0:
            ax.plot([cumulative, cumulative], [cost-1, cost+ 1],
                    color=color, lw=2, alpha=1)

        # vertical line (cost step) if not first generator
        if last_cost is not None:
            ax.plot([cumulative, cumulative], [last_cost, cost],
                    color='gray', lw=2)

        # horizontal block
        ax.plot([cumulative, cumulative + bid], [cost, cost],
                color=color, lw=3)

        # add label at midpoint of block
        counts_by_cost[cost] = counts_by_cost.get(cost, 0) + 1
        if bid > 0 and typ != 'polluter':
            ax.text(cumulative + bid/2, cost + 0.6, name, fontsize=fontsize, color=color,
                    ha='center')
       
       # store for inset
        x_vals.append(cumulative + bid)
        y_vals.append(cost)

        cumulative += bid
        last_cost = cost

    # --- Find market clearing price ---
    market_price = 0.0
    cumulative = 0.0
    for bid, cost, _, _ in sorted_bids:
        cumulative += bid
        if round(cumulative, 3) >= round(total_demand, 3):
            market_price = cost
            break

    # --- Draw market lines ---
    ax.axvline(total_demand, color='#2E6F40', linestyle='--', label='Total demand')
    ax.axhline(market_price, color='black', linestyle='--', label='Market price')
    ax.text(total_demand, market_price + 5, f'{round(total_demand)} MW', fontsize=fontsize, color="#2E6F40", ha='center')
    ax.text(-160, market_price, rf'$\lambda^{{DA *}}$',
            fontsize=fontsize, color='black')
    # --- Create inset axes for zoomed region ---
    axins = inset_axes(ax, width="45%", height="45%", loc='lower right', borderpad=4.5)
    cumulative_inset = 0.0
    last_cost = None
    counts_by_cost_inset = {}

    for bid, cost, typ, name in sorted_bids:
        if typ == "polluter":
            loc = wind_location(name)
            color = colors[loc]
        else:
            color = colors[typ]
        if last_cost is not None and cost == last_cost and cumulative_inset > 0:
            axins.plot([cumulative_inset, cumulative_inset], [cost-0.1, cost+ 0.1],
                       color=color, lw=2, alpha=1)
        if last_cost is not None:
            axins.plot([cumulative_inset, cumulative_inset], [last_cost, cost],
                       color='gray', lw=2)
        axins.plot([cumulative_inset, cumulative_inset + bid], [cost, cost],
                   color=color, lw=3)
        counts_by_cost_inset[cost] = counts_by_cost_inset.get(cost, 0) + 1
        if bid > 0:
            gen_id = ''.join(filter(str.isdigit, name))
            axins.text(cumulative_inset + bid/2, cost + 0.3, gen_id, fontsize=fontsize, color=color,
                        ha='center')
        cumulative_inset += bid
        last_cost = cost
    
    legend_elements = [
    Line2D([0], [0], color=colors["onshore"], lw=4, label="Onshore wind"),
    Line2D([0], [0], color=colors["offshore"], lw=4, label="Offshore wind"),
    ]
    
    # zoom limits: first ~20 MW (adjust as needed)
    axins.set_xlim(-10, 540)
    axins.set_ylim(0, 7)
    axins.set_xticks([0, 100, 200, 300, 400, 500])
    axins.set_yticks([0, 2, 4, 6])

    axins.tick_params(axis='x', labelsize=fontsize, rotation = 45)
    axins.tick_params(axis='y', labelsize=fontsize)
    axins.legend(handles=legend_elements, bbox_to_anchor=(0.35, 0.4), fontsize=fontsize-2)
    # Enable axes (ticks, labels)
#    axins.tick_params(axis='both', labelsize=10, direction='in', which='both')
    axins.grid(True, alpha = 0.3)

    # optional: connect inset with main plot
    mark_inset(ax, axins, loc1=3, loc2=2, fc="none", ec="0.7")

    ax.grid(True, alpha = 0.3)
    ax.set_xlim(-10, cumulative + 130)
    ax.set_ylim(0, 85)
    
    plt.tight_layout()

    return market_price, fig

# Function for power deviations grouped


def plot_power_deviations_grouped(
    model_names,
    deviations_list,
    polluters,
    base_color="Blues",
    figsize=(9.0, 6.0),
    # Geometry
    bar_width=0.55,
    # Fonts
    title_fontsize=16,
    label_fontsize=14,
    tick_fontsize=12,
    legend_fontsize=12,
    legend_title_fontsize=13,
    # Styling
    grid_alpha=0.25,
    separator_color="0.25",
    separator_lw=0.6,
    legend_title="Generators"
    ):
    """
    Grouped stacked bar chart for power deviations. 
    Models with identical results are merged into one bar.
    Professional formatting matched to the 2-model comparison style.
    """
    
    # 1. Grouping identical models
    groups = []  # (deviation_vector, list_of_names)
    seen_values = {}
    
    for name, dev in zip(model_names, deviations_list):
        dev_arr = np.asarray(dev, float).ravel()
        v_key = tuple(np.round(dev_arr, 4))
        
        if v_key in seen_values:
            groups[seen_values[v_key]][1].append(name)
        else:
            seen_values[v_key] = len(groups)
            groups.append((dev_arr, [name]))

    nP = len(polluters)
    M_unique = len(groups)
    x_positions = np.arange(M_unique)

    # 2. Colors (matching the professional sorted style)
    polluters = list(polluters)
    name_order = np.argsort([p.lower() for p in polluters])
    polluters_sorted = [polluters[i] for i in name_order]

    cmap = plt.get_cmap(base_color)
    xs = np.linspace(0.15, 0.90, nP) # Avoiding extreme white/dark
    colors_dict = {name: cmap(x) for name, x in zip(polluters_sorted, xs)}

    # 3. Figure & Styling
    fig, ax = plt.subplots(figsize=figsize)
    ax.set_facecolor("white")
    for spine in ax.spines.values():
        spine.set_visible(True)
        spine.set_linewidth(0.8)
        spine.set_color("0.3")
    fig.patch.set_facecolor("white")

    # 4. Drawing Stacked Bars
    all_pos_totals = []
    all_neg_totals = []

    for idx, (values, names) in enumerate(groups):
        pos_bottom, neg_bottom = 0.0, 0.0
        for val, p_name in zip(values, polluters):
            bottom = pos_bottom if val >= 0 else neg_bottom
            
            ax.bar(
                x_positions[idx],
                val,
                bar_width,
                bottom=bottom,
                color=colors_dict[p_name],
                edgecolor=separator_color,
                linewidth=separator_lw
            )
            
            if val >= 0: pos_bottom += val
            else: neg_bottom += val
        
        all_pos_totals.append(pos_bottom)
        all_neg_totals.append(neg_bottom)

    # 5. Formatting Axes
    ax.set_xticks(x_positions)
    # Combine model names for labels (e.g., "Q1, Q2")
    combined_labels = [", ".join(g[1]) for g in groups]
    ax.set_xticklabels(combined_labels, fontsize=tick_fontsize)
    
    ax.set_ylabel(r"Power Deviations $\Delta p_{g, \omega}$ [MW]", fontsize=label_fontsize)
    ax.set_title("Power Deviations Comparison Across Models", fontsize=title_fontsize, pad=15)
    
    # Y-limits with buffer
    ymax = max(all_pos_totals + [0])
    ymin = min(all_neg_totals + [0])
    buffer = 0.1 * max(abs(ymax), abs(ymin), 1)
    ax.set_ylim(ymin - buffer, ymax + buffer)

    ax.yaxis.grid(True, alpha=grid_alpha)
    ax.set_axisbelow(True)
    ax.tick_params(axis="both", labelsize=tick_fontsize)
    ax.axhline(0, color="0.3", linewidth=0.8) # Reference line at 0

    # 6. Legend (External like in your reference)
    handles = [Patch(facecolor=colors_dict[p], edgecolor=separator_color, label=p)
               for p in polluters_sorted]

    ax.legend(
        handles=handles,
        title=legend_title,
        loc="center left",
        bbox_to_anchor=(1.02, 0.5),
        frameon=False,
        fontsize=legend_fontsize,
        title_fontsize=legend_title_fontsize,
        handlelength=1.2,
        labelspacing=0.8,
    )

    plt.tight_layout()
    plt.show()
    return fig

# DA bids identical for diff models
def plot_DA_bids_grouped_identical(
    DA_bids,
    names_by_type=None,
    types_to_plot=None,
    tech_colors=None,
    bar_height=0.4,
    spacing=0.1,
    group_gap=0.6,          # Kontroluje przerwę między sekcjami
    sep_color="white",
    sep_lw=1.0,
    model_edgecolor="black",
    model_linewidth=1.1,
    sort_types_desc=True,
    sort_tech_desc=True,
    sort_units_desc=True,
    label_x_pad=0.01,
    base_fontsize=16,
    tick_fontsize=14,
    label_fontsize=16,
    title_fontsize=18,
    legend_fontsize=14,
    legend_title_fontsize=14,
    legend_loc="upper right",
    figsize=(10, 6)
    ):
    # --- 1. Konfiguracja stylów ---
    plt.rcParams.update({
        "font.size": base_fontsize,
        "axes.titlesize": title_fontsize,
        "axes.labelsize": label_fontsize,
        "xtick.labelsize": tick_fontsize,
        "ytick.labelsize": tick_fontsize,
        "legend.fontsize": legend_fontsize,
    })

    model_names = list(DA_bids.keys())
    if not model_names: raise ValueError("DA_bids jest pusty.")
    if names_by_type is None: raise ValueError("names_by_type jest wymagany.")

    all_types = sorted({gtype for model in DA_bids.values() for gtype in model.keys()})
    types_to_plot = [t for t in (types_to_plot or all_types) if t in all_types]

    def _type_total(gtype):
        return sum(np.sum(DA_bids[m].get(gtype, 0)) for m in model_names)

    if sort_types_desc:
        types_to_plot = sorted(types_to_plot, key=_type_total, reverse=True)

    # --- 2. Inicjalizacja wykresu ---
    fig, ax = plt.subplots(figsize=figsize)
    
    current_y = 0.0      # Startowa pozycja pionowa
    type_labels_y = []   # Pozycje dla etykiet "Base", "Deviators" itd.
    xmax = 0.0
    techs_present, seen_techs = [], set()

    # Rysujemy od góry do dołu, więc przechodzimy przez typy
    for gtype in types_to_plot:
        unit_names_base = list(names_by_type[gtype])
        
        # --- KROK A: GRUPOWANIE MODELI ---
        groups = []
        val_to_group_idx = {}
        for model in model_names:
            v_arr = np.asarray(DA_bids[model].get(gtype, [])).reshape(-1)
            v_key = tuple(np.round(v_arr, 4)) # Odporność na błędy numeryczne
            
            if v_key in val_to_group_idx:
                groups[val_to_group_idx[v_key]][1].append(model)
            else:
                val_to_group_idx[v_key] = len(groups)
                groups.append((v_arr, [model]))

        # Zapamiętujemy górną krawędź grupy
        group_top_y = current_y

        # --- KROK B: RYSOWANIE PASKÓW W GRUPIE ---
        for i, (vals, models_in_group) in enumerate(groups):
            # Obliczamy y dla konkretnego paska
            y = current_y - i * (bar_height + spacing)
            
            items = []
            for name, v in zip(unit_names_base, vals):
                if v > 1e-3: # Pomijamy zerowe oferty
                    items.append((name, float(v), tech_key(name)))

            tech_to_items = {}
            for name, v, tech in items:
                tech_to_items.setdefault(tech, []).append((name, v))

            sorted_techs = list(tech_to_items.keys())
            if sort_tech_desc:
                sorted_techs.sort(key=lambda t: sum(v for _, v in tech_to_items[t]), reverse=True)

            left = 0.0
            for tech in sorted_techs:
                if tech not in seen_techs:
                    techs_present.append(tech)
                    seen_techs.add(tech)
                
                unit_list = tech_to_items[tech]
                if sort_units_desc:
                    unit_list = sorted(unit_list, key=lambda x: x[1], reverse=True)

                for _, v in unit_list:
                    ax.barh(y, v, left=left, height=bar_height,
                            color=tech_colors.get(tech, "gray"), 
                            edgecolor=sep_color, linewidth=sep_lw)
                    left += v

            # Ramka modelu i etykieta tekstowa
            if left > 0:
                ax.barh(y, left, left=0.0, height=bar_height,
                        color="none", edgecolor=model_edgecolor, linewidth=model_linewidth)
                
                combined_names = ", ".join(models_in_group)
                # Dynamiczny pad na podstawie aktualnego xmax
                actual_pad = (xmax if xmax > 0 else left) * label_x_pad
                ax.text(left + actual_pad, y, f" {combined_names}", 
                        va="center", ha="left", fontsize=tick_fontsize, color="#948f8f")
                
                xmax = max(xmax, left)

        # --- KROK C: AKTUALIZACJA POZYCJI ---
        group_bottom_y = current_y - (len(groups) - 1) * (bar_height + spacing)
        # Środek grupy dla etykiety osi Y
        type_labels_y.append((group_top_y + group_bottom_y) / 2)
        
        # Przesuwamy current_y pod ostatni narysowany pasek plus przerwa między sekcjami
        current_y = group_bottom_y - (bar_height + spacing + group_gap)

    # --- 3. Finalizacja osi ---
    ax.set_yticks(type_labels_y)
    ax.set_yticklabels([t.replace("_", " ").title() for t in types_to_plot])

    ax.set_xlabel("DA Bids [MW]")
    ax.set_xlim(0, xmax * 1.3) # Miejsce na nazwy modeli
    ax.grid(axis="x", alpha=0.15)
    
    # Odwracamy osie, aby pierwszy typ był na górze
    ax.invert_yaxis()

    # Legenda technologii
    tech_handles = [
        Patch(facecolor=tech_colors.get(t, "gray"), label=pretty_tech_label(t))
        for t in techs_present
    ]
    if tech_handles:
        ax.legend(handles=tech_handles, title="Technology", loc=legend_loc, frameon=False)

    plt.tight_layout()
    plt.show()
    return fig


def tech_key(label: str) -> str:
    # LaTeX subscript: PREFIX$_{sub}$
    m = re.match(r"^(?P<prefix>[^$]+)\$_\{(?P<sub>[^}]*)\}\$$", str(label).strip())
    if m:
        prefix = m.group("prefix")
        sub = m.group("sub")

        # Deviators: good1/mid2/bad3
        m_sub = re.match(r"^(good|mid|bad)\d+$", sub)
        if m_sub:
            return f"{prefix}_{m_sub.group(1)}"

        # Indexed techs: SCGT$_{1}$
        if sub.isdigit():
            return prefix

    s = str(label).strip()

    # Plain deviators with index: On_good2 -> On_good
    m = re.match(r"^(?P<prefix>.+?)_(?P<qual>good|mid|bad)\d+$", s)
    if m:
        return f"{m.group('prefix')}_{m.group('qual')}"

    # Plain indexed tech: SCGT1 -> SCGT
    m = re.match(r"^(?P<tech>[A-Za-z_]+)\d+$", s)
    if m:
        return m.group("tech")

    return s

def pretty_tech_label(tech: str) -> str:
    """
    Convert tech keys to matplotlib-mathtext labels for the legend.
    Examples:
      On_good -> On$_{good}$
      Off_mid -> Off$_{mid}$
      CCGT -> CCGT
    """
    m = re.match(r"^(On|Off)_(good|mid|bad)$", tech)
    if m:
        prefix, qual = m.group(1), m.group(2)
        return f"{prefix}$_{{{qual}}}$"
    return tech


#plot to show the costs for many models:
def plot_socialized_vs_individualized_all_models_no_lines(
    pBIN_pen_by_model,         # dict: model -> (nP x 2) [up, down] per polluter (can be negative; we use abs)
    Q_up_by_model,             # dict: model -> scalar
    Q_down_by_model,           # dict: model -> scalar
    polluters,                 # list length nP
    model_order=None,          # e.g. [1,2,3,4,5,6,7]
    figsize=(16, 5.6),
    percent_fontsize=10,
    xticks_fontsize=12,

    # styling aligned with your FCB-only function
    anchors=("#EBA067", "#EE7518", "#D05C03"),
    color_order="as_is",       # "as_is" or "sorted"
    social_color="gray",
    prev_social_color="#928f8f",
    show_prev_all_socialized=True,

    bar_width=0.60,
    edgecolor="white",
    lw=0.6,
    grid_alpha=0.2,

    legend_title="Legend",
    legend_fontsize=12,
    legend_title_fontsize=13,
    ):
    if model_order is None:
        model_order = list(pBIN_pen_by_model.keys())

    polluters = list(map(str, polluters))

    # infer nP
    first = np.asarray(pBIN_pen_by_model[model_order[0]], dtype=float)
    if first.ndim != 2 or first.shape[1] != 2:
        raise ValueError("Each pBIN_pen_by_model[model] must be shape (nP, 2) [up, down].")
    nP = first.shape[0]
    if len(polluters) != nP:
        raise ValueError("polluters length must match nP.")

    # palette (same logic as your FCB-only plot)
    bar_colors = make_palette_from_anchors(list(anchors), nP)

    if color_order == "sorted":
        poll_sorted = sorted(polluters, key=lambda s: s.lower())
        color_by_name = {name: bar_colors[i] for i, name in enumerate(poll_sorted)}
        colors_for_polluters = [color_by_name[name] for name in polluters]
        legend_names = poll_sorted
        legend_colors = [color_by_name[name] for name in poll_sorted]
    else:
        colors_for_polluters = bar_colors
        legend_names = polluters
        legend_colors = bar_colors

    def _model_label(k):
        if isinstance(k, int) and k == 3:
            return r"$Q_{FCB}$"
        
        if isinstance(k, int) and k == 5:
            return r"$Q_{7}$"
        
        if isinstance(k, int):
            return rf"$Q_{{{k}}}$"
        
        return str(k)
    models = list(model_order)
    K = len(models)
    x = np.arange(K)

    # compute ymax across all models and both directions
    ymax = 1.0
    for m in models:
        Qup = float(Q_up_by_model[m-1])
        Qdn = float(Q_down_by_model[m-1])
        ymax = max(ymax, Qup, Qdn)
    ymax *= 1.20

    fig, axes = plt.subplots(1, 2, figsize=figsize, sharey=True)
    axes[0].set_ylabel("Penalty cost [€]")

    def _stack_models(ax, col_idx, title_txt, Q_dict):
        # optional background "100% socialized" (full Q) behind each model bar
        if show_prev_all_socialized:
            ax.bar(x, [float(Q_dict[m-1]) for m in models],
                   width=0.85, color=prev_social_color, edgecolor="none", zorder=0)

        bottoms = np.zeros(K)

        # stack polluters
        for p in range(nP):
            vals_p = []
            for m in models:
                mat = np.asarray(pBIN_pen_by_model[m], dtype=float)
                if mat.shape != (nP, 2):
                    raise ValueError(f"Model {m} has shape {mat.shape}, expected {(nP,2)}")
                vals_p.append(abs(float(mat[p, col_idx])))

            ax.bar(
                x, vals_p, width=bar_width, bottom=bottoms,
                color=colors_for_polluters[p],
                edgecolor=edgecolor, linewidth=lw, zorder=2
            )
            bottoms += np.array(vals_p)

        # add "socialized remainder" per model up to Q
        socials = np.array([max(float(Q_dict[m-1]) - bottoms[j], 0.0) for j, m in enumerate(models)], dtype=float)
        ax.bar(
            x, socials, width=bar_width, bottom=bottoms,
            color=social_color, edgecolor=edgecolor, linewidth=lw, zorder=2
        )

        # annotate socialized %
        for j, m in enumerate(models):
            Q = float(Q_dict[m-1])
            if Q > 0:
                pct = socials[j] / Q * 100
                ax.text(x[j], (bottoms[j] + socials[j]) * 1.01, f" Socialized \n {pct:.1f}%",
                        ha="center", va="bottom", fontsize=percent_fontsize)

        ax.set_title(title_txt)
        ax.set_xticks(x)
        ax.tick_params(axis='x', labelsize=xticks_fontsize)
        ax.tick_params(axis='y', labelsize=xticks_fontsize)
        ax.set_xticklabels([_model_label(m) for m in models])
        ax.set_ylim(0, ymax)
        ax.grid(axis="y", alpha=grid_alpha)

    _stack_models(axes[0], 0, "Up-regulation procurement cost", Q_up_by_model)
    _stack_models(axes[1], 1, "Down-regulation procurement cost", Q_down_by_model)

    # legend (no dashed-line handle)
    polluter_handles = [Patch(facecolor=legend_colors[i], label=legend_names[i]) for i in range(nP)]
    social_handle = [Patch(facecolor=social_color, label="Socialized")]
    prev_handle = [Patch(facecolor=prev_social_color, label="Previously:\n100% socialized")] if show_prev_all_socialized else []

    axes[1].legend(
        handles=polluter_handles + social_handle + prev_handle,
        title=legend_title,
        bbox_to_anchor=(1.02, 1.01),
        loc="upper left",
        frameon=True,
        fontsize=legend_fontsize,
        title_fontsize=legend_title_fontsize,
    )

    plt.tight_layout()
    plt.show()
    return fig 


import PATHSolver
PATHSolver.c_api_License_SetString("1259252040&Courtesy&&&USR&GEN2035&5_1_2026&1000&PATH&GEN&31_12_2035&0_0_0&6000&0_0")
using Pkg
Pkg.activate(".")
Pkg.instantiate()

using MPSGE
using DataFrames
using PlotlyJS
import JuMP.Containers: DenseAxisArray, @container
import JuMP
using Ipopt

include("Ramsey-dyn-4.jl")

#
#  DYN4 uses balanced growth

d4_data  = IndexedData();
d4_model = dynamic_dyn4_model(d4_data);

solve!(d4_model, cumulative_iteration_limit = 0)


goods        = d4_data.goods           # [:A, :B]
time_periods = d4_data.time_periods    # 2000:2100
QREF         = d4_data.QREF            # (1.02)^(t-2000) for all t


#  In GAMS DYN4 the benchmark is just QREF(T) because the balanced-growth
#  economy sits exactly on QREF.  We store the solved values so the
#  % deviation calculation uses the same basis as the GAMS output:
#    MACRO(T,"INVEST") = 100*(I.L(T)/QREF(T) - 1)

println("\n=== Step 2: Store benchmark values ===")
benchmark_values = DataFrame(
    t  = collect(time_periods),
    Y  = [value(d4_model[:Y][t])        for t in time_periods],
    C  = [value(d4_model[:C][t])        for t in time_periods],
    I  = [value(d4_model[:I][t])        for t in time_periods],
    K  = [value(d4_model[:K][t])        for t in time_periods],
    X1 = [value(d4_model[:X][:A, t])   for t in time_periods],
    X2 = [value(d4_model[:X][:B, t])   for t in time_periods],
) |> x -> stack(x, Not(:t), variable_name = :variable, value_name = :benchmark)

#  Policy shock: 10% tax on investment from 2010 onwards

println("\n=== Step 3: Apply 10% investment tax from 2010 and solve ===")
set_value!.(d4_model[:TAX][2010:end], 0.1)
solve!(d4_model)

shock_values = DataFrame(
    t  = collect(time_periods),
    Y  = [value(d4_model[:Y][t])        for t in time_periods],
    C  = [value(d4_model[:C][t])        for t in time_periods],
    I  = [value(d4_model[:I][t])        for t in time_periods],
    K  = [value(d4_model[:K][t])        for t in time_periods],
    X1 = [value(d4_model[:X][:A, t])   for t in time_periods],
    X2 = [value(d4_model[:X][:B, t])   for t in time_periods],
) |> x -> stack(x, Not(:t), variable_name = :variable, value_name = :shock)


#  Compute % deviations from benchmark

df = leftjoin(benchmark_values, shock_values, on = [:t, :variable]) |>
    x -> transform(x,
        [:benchmark, :shock] => ByRow((b, s) -> 100 * (s / b - 1)) => :value
    ) |>
    x -> select(x, :t, :variable, :value) |>
    x -> transform(x, :variable => ByRow(string) => :variable)

# ── Quick console summary (every 10 years) ───────────────────────────────
println("\n========== % Deviations from Benchmark (every 10 years) ==========")
println(rpad("Year", 6),
        rpad("Output",10), rpad("Invest",10),
        rpad("Cons",  10), rpad("Capital",10),
        rpad("X_A",   10), rpad("X_B",   10))
for t in 2000:10:2100
    row = Dict(
        r.variable => r.value
        for r in eachrow(subset(df, :t => ByRow(==(t))))
    )
    println(rpad(t, 6),
        rpad(round(get(row, "Y",  NaN), digits=3), 10),
        rpad(round(get(row, "I",  NaN), digits=3), 10),
        rpad(round(get(row, "C",  NaN), digits=3), 10),
        rpad(round(get(row, "K",  NaN), digits=3), 10),
        rpad(round(get(row, "X1", NaN), digits=3), 10),
        rpad(round(get(row, "X2", NaN), digits=3), 10),
    )
end


# Interactive Plotly chart
mode_order = ["Y", "C", "I", "K", "X1", "X2"]
modes = Dict(
    "Y"  => "Aggregate Output",
    "C"  => "Consumption",
    "I"  => "Investment",
    "K"  => "Capital",
    "X1" => "Sector A Output",
    "X2" => "Sector B Output",
)

layout = Layout(
    updatemenus = [
        attr(
            active  = 0,
            x       = 0.1,
            y       = 1.15,
            buttons = [
                attr(
                    label  = modes[mode],
                    method = "update",
                    args   = [
                        attr(visible = [m == mode for m in mode_order]),
                        attr(title   = modes[mode]),
                    ],
                ) for mode in mode_order
            ],
        )
    ],
    title        = modes[mode_order[1]],
    yaxis_title  = "% deviation from benchmark",
    xaxis_title  = "Year",
    legend_title = "Variable",
    template     = "plotly_white",
)

fig = plot(
    vec([
        scatter(
            subset(df, :variable => ByRow(==(mode))),
            x       = :t,
            y       = :value,
            mode    = "lines+markers",
            name    = modes[mode],
            visible = mode == mode_order[1],
        ) for mode in mode_order
    ]),
    layout,
)

display(fig)

# all variables on one plot
fig2 = plot(
    vec([
        scatter(
            subset(df, :variable => ByRow(==(mode))),
            x    = :t,
            y    = :value,
            mode = "lines",
            name = modes[mode],
        ) for mode in mode_order
    ]),
    Layout(
        title       = "% Deviation from Benchmark — All Variables",
        yaxis_title = "% deviation from benchmark",
        xaxis_title = "Year",
        template    = "plotly_white",
    )
)

display(fig2)
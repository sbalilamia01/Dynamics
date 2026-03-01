
struct IndexedData
    Y0::Float64
    C0::Float64
    LS0::Float64
    KS0::Float64
    I0::Float64

    goods::Vector{Symbol}
    X::DenseAxisArray
    LD::DenseAxisArray
    KD::DenseAxisArray

    G::Float64
    R::Float64
    D::Float64
    K0::Float64
    rk0::Float64
    pk0::Float64

    time_periods::UnitRange{Int}
    time_periods_horizon::UnitRange{Int}
    QREF::DenseAxisArray{Float64, 1, Tuple{UnitRange{Int}}}
    PREF::DenseAxisArray{Float64, 1, Tuple{UnitRange{Int}}}
end


function IndexedData(;
        Y0         = 200,
        C0         = 180,
        LS0        = 150,
        KS0        = 50,
        I0         = 20,
        goods      = [:A, :B],
        X          = DenseAxisArray([120.0,  80.0], goods),
        LD         = DenseAxisArray([100.0,  50.0], goods),
        KD         = DenseAxisArray([ 20.0,  30.0], goods),
        G          = 0.02,
        R          = 0.05,
        D          = 0.02,
        start_year = 2000,
        end_year   = 2100,
    )

    time_periods         = start_year:end_year
    time_periods_horizon = start_year:(end_year + 1)

    # Steady-state objects 
    K0   = I0 / (G + D)     
    rk0  = R + D              # in GAMS (R-R*D+D)/(1-R) 
    pk0  = 1 + R              # in GAMS 1/(1-R)

    QREF = DenseAxisArray((1 + G).^(eachindex(time_periods) .- 1),  time_periods)
    PREF = DenseAxisArray((1 / (1 + R)).^(eachindex(time_periods_horizon) .- 1), time_periods_horizon)

    adjusted_KS0 = K0 * rk0   # calibrated total capital earnings

    # Weighted-least-squares capital redistribution to match the balanced growth reference path
    M = JuMP.Model(Ipopt.Optimizer)
    JuMP.set_silent(M)
    JuMP.@variables(M, begin
        VK[g = goods] >= 0
        VL[g = goods] >= 0
    end)
    JuMP.@objective(M, Min,
        sum((1 / KD[g]) * (VK[g] - KD[g])^2 for g in goods))
    JuMP.@constraints(M, begin
        VABAL[g = goods], VK[g] + VL[g] == KD[g] + LD[g]
        VKBAL,            sum(VK[g] for g in goods) == adjusted_KS0
    end)
    JuMP.optimize!(M)

    adjusted_LD = JuMP.value.(VL)
    adjusted_KD = JuMP.value.(VK)

    @info "Raw factor demands"      LD = LD  KD = KD
    @info "Calibrated factor demands" LD = adjusted_LD  KD = adjusted_KD

    return IndexedData(
        Y0, C0, LS0, adjusted_KS0, I0,
        goods, X, adjusted_LD, adjusted_KD,
        G, R, D, K0, rk0, pk0,
        time_periods, time_periods_horizon,
        QREF, PREF,
    )
end


function dynamic_dyn4_model(data::IndexedData)

    Y0   = data.Y0
    C0   = data.C0
    I0   = data.I0
    KS0  = data.KS0         

    time_periods         = data.time_periods
    time_periods_horizon = data.time_periods_horizon

    D    = data.D
    K0   = data.K0
    pk0  = data.pk0
    r    = data.R

    QREF = data.QREF
    PREF = data.PREF

    goods = data.goods
    X0    = data.X
    LD    = data.LD
    KD    = data.KD


    DYN4 = MPSGEModel()


    @parameters(DYN4, begin
        TAX[t = time_periods], 0      
    end)

    @sectors(DYN4, begin
        Y[t = time_periods],          (start = QREF[t], description = "Macro output / transitory utility")
        X[g = goods, t = time_periods],(start = QREF[t], description = "production")
        I[t = time_periods],          (start = QREF[t], description = "Investment")
        K[t = time_periods],          (start = QREF[t], description = "Capital stock")
        C[t = time_periods],          (start = QREF[t], description = "Consumption index")
    end)

    @commodities(DYN4, begin
        PY[t = time_periods],          (start = PREF[t],      description = "Price of macro output")
        PX[g = goods, t = time_periods],(start = PREF[t],      description = "Sectoral output price")
        PC[t = time_periods],          (start = PREF[t],      description = "Consumption price index")
        RK[t = time_periods],          (start = PREF[t],      description = "Present-value return to capital")
        PL[t = time_periods],          (start = PREF[t],      description = "Present-value wage")
        PK[t = time_periods_horizon],  (start = PREF[t]*pk0,  description = "Price of capital stock")

    end)

    @consumer(DYN4, RA)

    @auxiliary(DYN4, TCAP,
        start = I0 * QREF[time_periods[end]] + K0 * (1 - D) * QREF[time_periods[end]])


    @production(DYN4, Y[t = time_periods], [t = 0, s = 1], begin
        @output(PY[t], Y0, t)
        @input(PX[g = goods, t], X0[g], s)
    end)


    @production(DYN4, X[g = goods, t = time_periods], [t = 0, s = 1], begin
        @output(PX[g, t], X0[g], t)
        @input(RK[t],    KD[g], s)
        @input(PL[t],    LD[g], s)
    end)

    @production(DYN4, I[t = time_periods], [t = 0, s = 0], begin
        @output(PK[t + 1], I0, t)
        @input(PY[t], I0, s, taxes = [Tax(RA, TAX[t])])
    end)

    @production(DYN4, K[t = time_periods], [t = 0, s = 0], begin
        @output(PK[t + 1], K0 * (1 - D), t)
        @output(RK[t],     KS0,           t)
        @input( PK[t],     K0,            s)
    end)

    @production(DYN4, C[t = time_periods], [t = 0, s = 0], begin
        @output(PC[t], C0, t)
        @input( PY[t], C0, s)
    end)


    @demand(DYN4, RA, begin
        @final_demand(PC[t = time_periods], C0 * QREF[t], reference_price = PREF[t])
        @endowment(PK[time_periods[begin]], K0)
        @endowment(PL[t = time_periods], sum(LD[g] for g in goods) * QREF[t])
        @endowment(PK[time_periods_horizon[end]], -TCAP)
    end)

    #  TRANSVERSALITY CONDITION
    @aux_constraint(DYN4, TCAP,
        C[time_periods[end - 1]] * I[time_periods[end]] -
        I[time_periods[end - 1]] * C[time_periods[end]]
    )

    return DYN4
end


function dyn4_report(model::MPSGEModel, data::IndexedData; baseline = nothing)
    time_periods = data.time_periods
    QREF         = data.QREF

    raw = DataFrame(
        t  = collect(time_periods),
        Y  = [value(model[:Y][t])        for t in time_periods],
        C  = [value(model[:C][t])        for t in time_periods],
        I  = [value(model[:I][t])        for t in time_periods],
        K  = [value(model[:K][t])        for t in time_periods],
        X1 = [value(model[:X][:A, t])   for t in time_periods],
        X2 = [value(model[:X][:B, t])   for t in time_periods],
    )

    return raw
end
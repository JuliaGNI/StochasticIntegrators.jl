"""
`SolutionPSDE`: Solution of a partitioned stochastic differential equation

Contains all fields necessary to store the solution of a PSDE or SPSDE

### Fields
* `conv`: type of the solution: :strong or :weak
* `nd`: dimension of the dynamical variable ``q``
* `nm`: dimension of the Wiener process
* `nt`: number of time steps to store
* `ns`: number of sample paths
* `ni`: number of initial conditions
* `t`:  time steps
* `q`:  solution `q[nd, nt+1, ns, ni]` with `q[:,0,:,:]` the initial conditions
* `p`:  solution `p[nd, nt+1, ns, ni]` with `p[:,0,:,:]` the initial conditions
* `W`:  Wiener process driving the stochastic processes q and p
* `K`:  integer parameter defining the truncation of the increments of the Wiener process (for strong solutions),
*       A = √(2 K Δt |log Δt|) due to Milstein & Tretyakov; if K=0 no truncation
* `ntime`: number of time steps to compute
* `nsave`: save every nsave'th time step

"""
abstract type SolutionPSDE{dType, tType, wType, NQ, NW, CONV} <: StochasticSolution{dType, tType, wType, NQ, NW} end

# Create SolutionPSDEs with serial and parallel data structures.
for (TSolution, TDataSeries, Tdocstring) in
    ((:SSolutionPSDE, :SDataSeries, "Serial Solution of a partitioned stochastic differential equation."),
     (:PSolutionPSDE, :PDataSeries, "Parallel Solution of a partitioned stochastic differential equation."))
    @eval begin
        $Tdocstring
        mutable struct $TSolution{dType, tType, wType, NQ, NW, CONV} <: SolutionPSDE{dType, tType, wType, NQ, NW, CONV}
            nd::Int
            nm::Int
            nt::Int
            ns::Int
            t::TimeSeries{tType}
            q::$TDataSeries{dType,NQ}
            p::$TDataSeries{dType,NQ}
            W::WienerProcess{wType,tType,NW,CONV}
            K::Int
            ntime::Int
            nsave::Int
            nwrite::Int
            counter::Vector{Int}
            woffset::Int
            ioffset::Int
            periodicity::dType
            h5::HDF5.File

            function $TSolution(t::TimeSeries{TT}, q::$TDataSeries{DT,NQ}, p::$TDataSeries{DT,NQ}, W::WienerProcess{WT,TT,NW,CONV}; K::Int=0) where {DT,TT,WT,NQ,NW,CONV}
                # extract parameters
                nd = q.nd
                ns = q.ni
                nt = q.nt
                nm = W.nd
                ntime = W.nt
                nsave = t.step
                nwrite = ntime

                @assert CONV==:strong || (CONV==:weak && K==0) || CONV==:null
                @assert ntime==nt*nsave
                @assert q.ni == p.ni
                @assert q.nd == p.nd
                @assert q.nt == p.nt == t.n
                @assert W.ns == q.ni

                new{DT,TT,WT,NQ,NW,CONV}(nd, nm, nt, ns, t, q, p, W, K, ntime, nsave, nwrite, zeros(Int, ns), 0, 0)
            end

            function $TSolution(::Type{dType}, nd::Int, nm::Int, nt::Int, ns::Int, ni::Int, Δt::tType,
                        W::WienerProcess{wType,tType,NW,CONV}, K::Int, ntime::Int, nsave::Int, nwrite::Int, periodicity=zeros(dType, nd)) where {dType <:  Union{Number,AbstractArray}, tType <: Real, wType <: Number, NW, CONV}

                @assert CONV==:strong || (CONV==:weak && K==0) || CONV==:null

                @assert nsave > 0
                @assert ntime == 0 || ntime ≥ nsave
                @assert nwrite == 0 || nwrite ≥ nsave
                @assert mod(ntime, nsave) == 0

                if nwrite > 0
                    @assert mod(nwrite, nsave) == 0
                    @assert mod(ntime, nwrite) == 0
                end

                @assert nd > 0
                @assert ns > 0
                @assert ni > 0
                @assert ni == 1 || ns == 1
                @assert nt ≥ 0

                @assert NW ∈ (2,3)

                t = TimeSeries{tType}(nt, Δt, nsave)
                q = $TDataSeries(dType, nd, nt, max(ns,ni))
                p = $TDataSeries(dType, nd, nt, max(ns,ni))
                NQ = ns==ni==1 ? 1 : 2

                new{Vector{dType}, tType, wType, NQ, NW, CONV}(nd, nm, nt, max(ns,ni), t, q, p, W, K, ntime, nsave, nwrite, zeros(Int, max(ns,ni)), 0, 0, periodicity)
            end
        end


        function $TSolution(equation::Union{PSDE{DT,TT},SPSDE{DT,TT}}, Δt::TT, ntime::Int; nsave::Int=DEFAULT_NSAVE, nwrite::Int=DEFAULT_NWRITE, K::Int=0, conv=DEFAULT_SCONV, filename=nothing) where {DT,TT}
            nd = equation.d
            nm = equation.m
            ns = equation.ns
            ni = nsamples(equation)
            nt = div(ntime, nsave)
            nt = (nwrite == 0 ? nt : div(nwrite, nsave))
            nw = (nwrite == 0 ? ntime : nwrite)

            # Holds the Wiener process data for ALL computed time steps
            # Wiener process increments are automatically generated here
            W = WienerProcess(DT, nm, nw, max(ni,ns), Δt, conv)

            s = $TSolution(DT, nd, nm, nt, ns, ni, Δt, W, K, ntime, nsave, nw, periodicity(equation))
            set_initial_conditions!(s, equation)

            if !isnothing(filename)
                isfile(filename) ? @warn("Overwriting existing HDF5 file.") : nothing
                create_hdf5!(s, filename)
            end

            return s
        end


        function $TSolution(equation::Union{PSDE{DT,TT},SPSDE{DT,TT}}, Δt::TT, dW::Array{DT, NW}, dZ::Array{DT, NW}, ntime::Int; nsave::Int=DEFAULT_NSAVE, nwrite::Int=DEFAULT_NWRITE, K::Int=0, conv=DEFAULT_SCONV, filename=nothing) where {DT,TT,NW}
            nd = equation.d
            nm = equation.m
            ns = equation.ns
            ni = nsamples(equation)
            nt = div(ntime, nsave)
            nt = (nwrite == 0 ? nt : div(nwrite, nsave))
            nw = (nwrite == 0 ? ntime : nwrite)

            @assert size(dW) == size(dZ)
            @assert NW ∈ (2,3)

            @assert nm == size(dW,1)
            @assert ntime == size(dW,2)
            @assert max(ni,ns) == size(dW,3)

            # Holds the Wiener process data for ALL computed time steps
            # Wiener process increments are prescribed by the arrays dW and dZ
            W = WienerProcess(Δt, dW, dZ, conv)

            s = $TSolution(DT, nd, nm, nt, ns, ni, Δt, W, K, ntime, nsave, nw, periodicity(equation))
            set_initial_conditions!(s, equation)

            if !isnothing(filename)
                isfile(filename) ? @warn("Overwriting existing HDF5 file.") : nothing
                create_hdf5!(s, filename)
            end

            return s
        end


        # If the Wiener process W data are not available, creates a one-element zero array instead
        # For instance used when reading a file with no Wiener process data saved
        function $TSolution(t::TimeSeries{TT}, q::$TDataSeries{DT,NQ}, p::$TDataSeries{DT,NQ}; K::Int=0, conv=DEFAULT_SCONV) where {DT,TT,NQ}
            # extract parameters
            nd = q.nd
            ns = q.ni
            nt = t.n
            nsave = t.step

            ΔW = (ns == 1 ? zeros(DT,0,0) : zeros(DT,0,0,0))
            ΔZ = (ns == 1 ? zeros(DT,0,0) : zeros(DT,0,0,0))

            W = WienerProcess{DT,TT,ns == 1 ? 2 : 3,:null}(nd, nt ,ns, t.Δt, ΔW, ΔZ)

            # create solution
            $TSolution(t, q, p, W, K=K)
        end


        function $TSolution(file::String)
            # open HDF5 file
            get_config(:verbosity) > 1 ? @info("Reading HDF5 file ", file) : nothing
            h5 = h5open(file, "r")

            # read attributes
            ntime = read(attributes(h5)["ntime"])
            nsave = read(attributes(h5)["nsave"])

            # reading data arrays
            t = TimeSeries(read(h5["t"]), nsave)

            if haskey(attributes(h5),"conv")
                conv = Symbol(read(attributes(h5)["conv"]))
            else
                conv = DEFAULT_SCONV
            end

            W_exists = haskey(h5, "ΔW") && haskey(h5, "ΔZ")

            if W_exists == true
                W = WienerProcess(t.Δt, read(h5["ΔW"]), read(h5["ΔZ"]), conv)
            end

            if haskey(attributes(h5),"K")
                K = read(attributes(h5)["K"])
            else
                K=0
            end

            q_array = read(h5["q"])
            p_array = read(h5["p"])

            close(h5)

            q = $TDataSeries(q_array)
            p = $TDataSeries(p_array)

            # create solution
            if W_exists == true
                $TSolution(t, q, p, W; K=K)
            else
                $TSolution(t, q, p; K=K, conv=conv)
            end
        end
    end
end


Base.:(==)(sol1::SolutionPSDE{DT1,TT1,NQ1,NW1,C1}, sol2::SolutionPSDE{DT2,TT2,NQ2,NW2,C2}) where {DT1,TT1,NQ1,NW1,C1,DT2,TT2,NQ2,NW2,C2} = (
                                DT1 == DT2
                             && TT1 == TT2
                             && NQ1 == NQ2
                             && NW1 == NW2
                             && C1  == C2
                             && sol1.nd == sol2.nd
                             && sol1.nm == sol2.nm
                             && sol1.nt == sol2.nt
                             && sol1.ns == sol2.ns
                             && sol1.t  == sol2.t
                             && sol1.q  == sol2.q
                             && sol1.p  == sol2.p
                             && sol1.W  == sol2.W
                             && sol1.K  == sol2.K
                             && sol1.ntime == sol2.ntime
                             && sol1.nsave == sol2.nsave
                             && sol1.nwrite == sol2.nwrite
                             && sol1.counter == sol2.counter
                             && sol1.woffset == sol2.woffset
                             && sol1.periodicity == sol2.periodicity)

@inline Solutions.hdf5(sol::SolutionPSDE) = sol.h5
@inline Solutions.timesteps(sol::SolutionPSDE)  = sol.t
@inline Solutions.nsave(sol::SolutionPSDE) = sol.nsave
@inline Solutions.counter(sol::SolutionPSDE) = sol.counter
@inline Solutions.offset(sol::SolutionPSDE) = sol.woffset
@inline ioffset(sol::SolutionPSDE) = sol.ioffset
@inline Solutions.lastentry(sol::SolutionPSDE) = sol.ni == 1 ? sol.counter[1] - 1 : sol.counter .- 1
@inline conv(sol::SolutionPSDE{DT,TT,NQ,NW,CONV}) where {DT,TT,NQ,NW,CONV} = CONV
@inline Common.ntime(sol::SolutionPSDE) = sol.ntime
@inline Common.periodicity(sol::SolutionPSDE) = sol.periodicity


function Solutions.set_initial_conditions!(sol::SolutionPSDE, equ::Union{PSDE,SPSDE})
    set_initial_conditions!(sol, equ.t₀, equ.q₀, equ.p₀)
end

function Solutions.set_initial_conditions!(sol::SolutionPSDE{AT,TT,WT,1}, t₀::TT, q₀::AT, p₀::AT) where {AT,TT,WT}
    set_data!(sol.q, q₀, 0)
    set_data!(sol.p, p₀, 0)
    compute_timeseries!(sol.t, t₀)
    sol.counter .= 1
end

function Solutions.set_initial_conditions!(sol::SolutionPSDE{AT,TT,WT,2}, t₀::TT, q₀::AT, p₀::AT) where {AT,TT,WT}
    for k in 1:sol.ns
        set_data!(sol.q, q₀, 0, k)
        set_data!(sol.p, p₀, 0, k)
    end
    compute_timeseries!(sol.t, t₀)
    sol.counter .= 1
end

function Solutions.set_initial_conditions!(sol::SolutionPSDE{AT,TT,WT}, t₀::TT, q₀::AbstractVector{AT}, p₀::AbstractVector{AT}) where {AT,TT,WT}
    # Sets the initial conditions sol.q[0] with the data from q₀
    # Here, q₀ is 1D (nd elements) representing a single deterministic or
    # multiple random initial conditions.
    # Similar for sol.p[0].
    if length(eachindex(q₀)) == length(eachindex(p₀)) == 1
        for k in 1:sol.ns
            set_data!(sol.q, q₀[begin], 0, k)
            set_data!(sol.p, p₀[begin], 0, k)
        end
    else
        @assert length(eachindex(q₀)) == length(axes(sol.q,2))
        @assert length(eachindex(p₀)) == length(axes(sol.p,2))
        for k in eachindex(q₀,p₀)
            set_data!(sol.q, q₀[k], 0, k)
            set_data!(sol.p, p₀[k], 0, k)
        end
    end
    compute_timeseries!(sol.t, t₀)
    sol.counter .= 1
end


function Solutions.get_initial_conditions!(sol::SolutionPSDE{AT,TT}, asol::AtomicSolutionPSDE{DT,TT,AT}, k, n=1) where {DT, TT, AT <: AbstractArray{DT}}
    get_solution!(sol, asol.q, asol.p, n-1, k)
    asol.t  = sol.t[n-1]
    asol.q̃ .= 0
    asol.p̃ .= 0
end

function Solutions.get_initial_conditions!(sol::SolutionPSDE{AT}, q::AT, p::AT, k, n=1) where {AT}
    get_solution!(sol, q, p, n-1, k)
end

function Solutions.get_initial_conditions(sol::SolutionPSDE, k, n=1)
    get_solution(sol, n-1, k)
end


function Solutions.set_solution!(sol::SolutionPSDE, t, q, p, n, k=1)
    set_solution!(sol, q, p, n, k)
end

function Solutions.set_solution!(sol::SolutionPSDE{AT,TT}, asol::AtomicSolutionPSDE{DT,TT,AT}, n, k=1) where {DT, TT, AT <: AbstractArray{DT}}
    set_solution!(sol, asol.t, asol.q, asol.p, n, k)
end

function Solutions.set_solution!(sol::SolutionPSDE{AT}, q::AT, p::AT, n, k=1) where {AT}
    @assert n <= sol.ntime
    @assert k <= sol.ns
    if mod(n, sol.nsave) == 0
        if sol.counter[k] > sol.nt
            @error("Solution overflow. Call write_to_hdf5() and reset!() before continuing the simulation.")
        end
        set_data!(sol.q, q, sol.counter[k], k)
        set_data!(sol.p, p, sol.counter[k], k)
        sol.counter[k] += 1
    end
end


function Solutions.get_solution!(sol::SolutionPSDE{AT}, q::AT, p::AT, n, k=1) where {AT}
    get_data!(sol.q, q, n, k)
    get_data!(sol.p, p, n, k)
end

function Solutions.get_solution(sol::SolutionPSDE{AT,TT,WT,1}, n, k=1) where {AT,TT,WT}
    @assert k == 1
    (sol.t[n], sol.q[n], sol.p[n])
end

function Solutions.get_solution(sol::SolutionPSDE{AT,TT,WT,2}, n, k=1) where {AT,TT,WT}
    (sol.t[n], sol.q[n,k], sol.p[n,k])
end


# copy increments of the Brownian Process for multidimensional Brownian motion, 1 sample path
function get_increment(sol::SolutionPSDE{AT,TT,WT,NQ,2}, n, k=1) where {AT,TT,WT,NQ}
    @assert k==1
    return (sol.W.ΔW[:,n], sol.W.ΔZ[:,n])
end

# copy increments of the Brownian Process for multidimensional Brownian motion, r-th sample path
function get_increment(sol::SolutionPSDE{AT,TT,WT,NQ,3}, n, k) where {AT,TT,WT,NQ}
    return (sol.W.ΔW[:,n,k], sol.W.ΔZ[:,n,k])
end

# copy increments of the Brownian Process for multidimensional Brownian motion, 1 sample path
function get_increments!(sol::SolutionPSDE{AT,TT,WT,NQ,2}, asol::AtomicSolutionPSDE{DT,TT,AT}, n, k=1) where {DT,TT,AT,WT,NQ}
    @assert k==1
    for l = 1:sol.nm
        asol.ΔW[l] = sol.W.ΔW[l,n]
        asol.ΔZ[l] = sol.W.ΔZ[l,n]
    end
end

# copy increments of the Brownian Process for multidimensional Brownian motion, r-th sample path
function get_increments!(sol::SolutionPSDE{AT,TT,WT,NQ,3}, asol::AtomicSolutionPSDE{DT,TT,AT}, n, k) where {DT,TT,AT,WT,NQ}
    for l = 1:sol.nm
        asol.ΔW[l] = sol.W.ΔW[l,n,k]
        asol.ΔZ[l] = sol.W.ΔZ[l,n,k]
    end
end


function Common.reset!(sol::SolutionPSDE)
    reset!(sol.q)
    reset!(sol.p)
    compute_timeseries!(sol.t, sol.t[end])
    generate_wienerprocess!(sol.W)
    sol.counter .= 1
    sol.woffset += sol.nt
    sol.ioffset += sol.nwrite
end

function Base.close(solution::SolutionPSDE)
    close(solution.h5)
end


"save_attributes: Saves attributes of Stochastic Solutions to HDF5 file."
function Solutions.save_attributes(h5::HDF5.File, sol::StochasticSolution)
    attributes(h5)["ntime"] = ntime(sol)
    attributes(h5)["nsave"] = nsave(sol)
    attributes(h5)["conv"]  = string(conv(sol))
    attributes(h5)["nd"] = sol.nd
    attributes(h5)["nm"] = sol.nm
    attributes(h5)["ns"] = sol.ns
    attributes(h5)["K"]  = sol.K
end


function Solutions.init_solution(h5::HDF5.File, solution::SolutionSDE{DT,TT,2}) where {DT,TT}
    q = create_dataset(h5, "q", DT, ((solution.nd, solution.nt+1), (solution.nd, -1)), chunk=(solution.nd,1))
    q[:,1] = solution.q[:,0]
end

function Solutions.init_solution(h5::HDF5.File, solution::SolutionSDE{DT,TT,3}) where {DT,TT}
    q = create_dataset(h5, "q", DT, ((solution.nd, solution.nt+1, solution.ns),(solution.nd, -1, solution.ns)), chunk=(solution.nd,1,1))
    q[:,1,:] = solution.q[:,0,:]
end

function Solutions.init_solution(h5::HDF5.File, solution::SolutionPSDE{DT,TT,2}) where {DT,TT}
    q = create_dataset(h5, "q", DT, ((solution.nd, solution.nt+1), (solution.nd, -1)), chunk=(solution.nd,1))
    p = create_dataset(h5, "p", DT, ((solution.nd, solution.nt+1), (solution.nd, -1)), chunk=(solution.nd,1))
    q[:, 1] = solution.q[:, 0]
    p[:, 1] = solution.p[:, 0]
end

function Solutions.init_solution(h5::HDF5.File, solution::SolutionPSDE{DT,TT,3}) where {DT,TT}
    q = create_dataset(h5, "q", DT, ((solution.nd, solution.nt+1, solution.ns),(solution.nd, -1, solution.ns)), chunk=(solution.nd,1,1))
    p = create_dataset(h5, "p", DT, ((solution.nd, solution.nt+1, solution.ns),(solution.nd, -1, solution.ns)), chunk=(solution.nd,1,1))
    q[:, 1, :] = solution.q[:, 0, :]
    p[:, 1, :] = solution.p[:, 0, :]
end


function init_increments(h5::HDF5.File, solution::StochasticSolution{DT,TT,NQ,2}) where {DT,TT,NQ}
    dW = create_dataset(h5, "ΔW", DT, ((solution.nm, solution.ntime),(solution.nm, -1)), chunk=(solution.nm,1))
    dZ = create_dataset(h5, "ΔZ", DT, ((solution.nm, solution.ntime),(solution.nm, -1)), chunk=(solution.nm,1))
end

function init_increments(h5::HDF5.File, solution::StochasticSolution{DT,TT,NQ,3}) where {DT,TT,NQ}
    dW = create_dataset(h5, "ΔW", DT, ((solution.nm, solution.ntime, solution.ns),(solution.nm, -1, solution.ns)), chunk=(solution.nm,1,1))
    dZ = create_dataset(h5, "ΔZ", DT, ((solution.nm, solution.ntime, solution.ns),(solution.nm, -1, solution.ns)), chunk=(solution.nm,1,1))
end


"Creates HDF5 file and initialises datasets for stochastic solution object."
function create_hdf5(solution::StochasticSolution, file::AbstractString; save_W=true)
    # create HDF5 file and save attributes and common parameters
    h5 = createHDF5(solution, file)

    # save attributes
    save_attributes(h5, solution)

    # create dataset
    init_timeteps(h5, solution)
    init_solution(h5, solution)

    if save_W
        init_increments(h5, solution)
    end

    return h5
end


function Solutions.save_solution(h5::HDF5.File, solution::SolutionSDE{DT,TT,1}, j1, j2, n1, n2) where {DT <: Number, TT}
    if size(h5["q"])[end] < j2
        HDF5.set_extent_dims(h5["q"], (j2,))
    end
    h5["q"][j1:j2] = solution.q[n1:n2]
end

function Solutions.save_solution(h5::HDF5.File, solution::SolutionSDE{AT,TT,1}, j1, j2, n1, n2) where {DT, AT <: Array{DT}, TT}
    elaxes = axes(solution.q[begin])
    if size(h5["q"])[end] < j2
        HDF5.set_extent_dims(h5["q"], (size(h5["q"])[1:end-1]..., j2))
    end
    for i in eachindex(j1:j2, n1:n2)
        j = (j1:j2)[i]
        n = (n1:n2)[i]
        h5["q"][elaxes..., j] = solution.q[n]
    end
end

function Solutions.save_solution(h5::HDF5.File, solution::SolutionSDE{DT,TT,2}, j1, j2, n1, n2) where {DT <: Number, TT}
    if size(h5["q"],1) < j2
        HDF5.set_extent_dims(h5["q"], (j2, size(h5["q"],2)))
    end
    h5["q"][j1:j2, :] = solution.q[n1:n2, :]
end

function Solutions.save_solution(h5::HDF5.File, solution::SolutionSDE{AT,TT,2}, j1, j2, n1, n2) where {DT, AT <: Array{DT}, TT}
    elaxes = axes(solution.q[begin,begin])
    if size(h5["q"])[end-1] < j2
        HDF5.set_extent_dims(h5["q"], (size(h5["q"])[begin:end-2]..., j2, size(h5["q"])[end]))
    end
    for k in 1:nsamples(solution.q)
        for i in eachindex(j1:j2, n1:n2)
            j = (j1:j2)[i]
            n = (n1:n2)[i]
            h5["q"][elaxes..., j, k] = solution.q[n, k]
        end
    end
end


function Solutions.save_solution(h5::HDF5.File, solution::SolutionPSDE{DT,TT,1}, j1, j2, n1, n2) where {DT <: Number, TT}
    if size(h5["q"])[end] < j2
        HDF5.set_extent_dims(h5["q"], (j2,))
    end
    if size(h5["p"])[end] < j2
        HDF5.set_extent_dims(h5["p"], (j2,))
    end
    h5["q"][j1:j2] = solution.q[n1:n2]
    h5["p"][j1:j2] = solution.p[n1:n2]
end

function Solutions.save_solution(h5::HDF5.File, solution::SolutionPSDE{AT,TT,1}, j1, j2, n1, n2) where {DT, AT <: Array{DT}, TT}
    elaxes = axes(solution.q[begin])
    if size(h5["q"])[end] < j2
        HDF5.set_extent_dims(h5["q"], (size(h5["q"])[1:end-1]..., j2))
    end
    for i in eachindex(j1:j2, n1:n2)
        j = (j1:j2)[i]
        n = (n1:n2)[i]
        h5["q"][elaxes..., j] = solution.q[n]
    end

    elaxes = axes(solution.p[begin])
    if size(h5["p"])[end] < j2
        HDF5.set_extent_dims(h5["p"], (size(h5["p"])[1:end-1]..., j2))
    end
    for i in eachindex(j1:j2, n1:n2)
        j = (j1:j2)[i]
        n = (n1:n2)[i]
        h5["p"][elaxes..., j] = solution.p[n]
    end
end

function Solutions.save_solution(h5::HDF5.File, solution::SolutionPSDE{DT,TT,2}, j1, j2, n1, n2) where {DT <: Number, TT}
    if size(h5["p"],1) < j2
        HDF5.set_extent_dims(h5["p"], (j2, size(h5["p"],2)))
    end
    if size(h5["p"],1) < j2
        HDF5.set_extent_dims(h5["p"], (j2, size(h5["p"],2)))
    end
    h5["q"][j1:j2, :] = solution.q[n1:n2, :]
    h5["p"][j1:j2, :] = solution.p[n1:n2, :]
end

function Solutions.save_solution(h5::HDF5.File, solution::SolutionPSDE{AT,TT,2}, j1, j2, n1, n2) where {DT, AT <: Array{DT}, TT}
    elaxes = axes(solution.q[begin,begin])
    if size(h5["q"])[end-1] < j2
        HDF5.set_extent_dims(h5["q"], (size(h5["q"])[begin:end-2]..., j2, size(h5["q"])[end]))
    end
    for k in 1:nsamples(solution.q)
        for i in eachindex(j1:j2, n1:n2)
            j = (j1:j2)[i]
            n = (n1:n2)[i]
            h5["q"][elaxes..., j, k] = solution.q[n, k]
        end
    end

    elaxes = axes(solution.p[begin,begin])
    if size(h5["p"])[end-1] < j2
        HDF5.set_extent_dims(h5["p"], (size(h5["p"])[begin:end-2]..., j2, size(h5["p"])[end]))
    end
    for k in 1:nsamples(solution.p)
        for i in eachindex(j1:j2, n1:n2)
            j = (j1:j2)[i]
            n = (n1:n2)[i]
            h5["p"][elaxes..., j, k] = solution.p[n, k]
        end
    end
end


function save_increments(h5::HDF5.File, solution::StochasticSolution{DT,TT,NQ,2}, j1, j2, n1, n2) where {DT,TT,NQ}
    if size(h5["ΔW"],2) < j2
        HDF5.set_extent_dims(h5["ΔW"], (size(h5["ΔW"],1), j2))
    end
    if size(h5["ΔZ"],2) < j2
        HDF5.set_extent_dims(h5["ΔZ"], (size(h5["ΔZ"],1), j2))
    end
    h5["ΔW"][:,j1:j2] = solution.W.ΔW[:,n1:n2]
    h5["ΔZ"][:,j1:j2] = solution.W.ΔZ[:,n1:n2]
end

function save_increments(h5::HDF5.File, solution::StochasticSolution{DT,TT,NQ,3}, j1, j2, n1, n2) where {DT,TT,NQ}
    if size(h5["ΔW"],2) < j2
        HDF5.set_extent_dims(h5["ΔW"], (size(h5["ΔW"],1), j2, size(h5["ΔW"],3)))
    end
    if size(h5["ΔZ"],2) < j2
        HDF5.set_extent_dims(h5["ΔZ"], (size(h5["ΔZ"],1), j2, size(h5["ΔZ"],3)))
    end
    h5["ΔW"][:,j1:j2,:] = solution.W.ΔW[:,n1:n2,:]
    h5["ΔZ"][:,j1:j2,:] = solution.W.ΔZ[:,n1:n2,:]
end


"""
Append solution to HDF5 file.
  soffset - start writing the solution q at the position soffset+2
  woffset - start writing the increments ΔW, ΔZ at the position woffset+1
"""
function write_to_hdf5(solution::StochasticSolution, h5::HDF5.File=hdf5(solution), soffset=offset(solution), woffset=ioffset(solution))
    # set convenience variables and compute ranges
    js1 = soffset+2
    js2 = soffset+1+solution.nt
    jw1 = woffset+1
    jw2 = woffset+solution.nwrite

    # copy data from solution to HDF5 dataset
    save_timeteps(h5, solution, js1, js2, 1, solution.nt)
    save_solution(h5, solution, js1, js2, 1, solution.nt)

    if exists(h5, "ΔW") && exists(h5, "ΔZ")
        # copy the Wiener process increments from solution to HDF5 dataset
        save_increments(h5, solution, jw1, jw2, 1, solution.nwrite)
    end

    return nothing
end


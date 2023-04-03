using PythonPlot
using .plotting

function upsample(v::AbstractVector, fac::Int,
                  do_plot::Bool = false)
    dv = v[2]-v[1]
    all(abs.(diff(v) - dv)/abs(dv) .<= 1e-6) || error("Not equal spacing between all elements")
    ndv = dv/fac
    nv = v[1]:ndv:v[end]+0.9ndv

    if do_plot
        cfigure("upsample") do
            csubplot(311) do
                plot(v, v, "s-")
                plot(nv, nv, ".")
                axis([0.9v[1], 1.1v[10], 0.9v[1], 1.1v[10]])
            end
            csubplot(312) do
                plot(v, v, "s-")
                plot(nv, nv, ".")
                axis([0.9v[end-10], 1.1v[end], 0.9v[end-10], 1.1v[end]])
            end
            csubplot(313) do
                semilogy(abs(v-nv[1:fac:length(nv)]))
            end
        end
    end

    nsel = filter(eachindex(nv)) do i
        !(i in 1:fac:length(nv))
    end
    nv, nv[nsel]
end

function log_upsample(v::AbstractVector, fac::Int)
    nv,nvv = upsample(log.(v), fac)
    exp.(nv),exp.(nvv)
end

export upsample, log_upsample

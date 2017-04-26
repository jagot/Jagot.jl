using PyPlot

ind{T<:Number,U<:Number}(v::AbstractVector{T}, a::U) = indmin(abs(v-a))

function upsample(v::AbstractVector, fac::Int,
                  do_plot::Bool = false)
    dv = v[2]-v[1]
    all(abs(diff(v) - dv)/abs(dv) .<= 1e-6) || error("Not equal spacing between all elements")
    ndv = dv/fac
    nv = v[1]:ndv:v[end]+0.9ndv

    if do_plot
        figure("upsample")
        clf()
        subplot(311)
        plot(v, v, "s-")
        plot(nv, nv, ".")
        axis([0.9v[1], 1.1v[10], 0.9v[1], 1.1v[10]])
        subplot(312)
        plot(v, v, "s-")
        plot(nv, nv, ".")
        axis([0.9v[end-10], 1.1v[end], 0.9v[end-10], 1.1v[end]])
        subplot(313)
        semilogy(abs(v-nv[1:fac:length(nv)]))
    end

    nsel = filter(eachindex(nv)) do i
        !(i in 1:fac:length(nv))
    end
    nv, nv[nsel]
end

function log_upsample(v::AbstractVector, fac::Int)
    nv,nvv = upsample(log(v), fac)
    exp(nv),exp(nvv)
end

meshgrid(x::AbstractVector,y::AbstractVector) =
    repeat(x', outer=(length(y),1)),repeat(y, outer=(1,length(x)))

export ind, upsample, log_upsample, meshgrid

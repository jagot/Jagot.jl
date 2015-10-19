module plotting
using PyPlot

include("colormaps.jl")

function plot_map(args...; kwargs...)
    keys = map(first,kwargs)
    function push_default!(key, val)
        if !(key in keys)
            push!(kwargs, (key, val))
        end
    end
    push_default!(:cmap, colormaps.viridis)
    push_default!(:rasterized, true)
    p = pcolormesh(args...; kwargs...)
    margins(0,0)
    p
end

export colormaps, plot_map

end

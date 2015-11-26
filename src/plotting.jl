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
    push_default!(:cmap, plotting.get_cmap("viridis"))
    push_default!(:rasterized, true)
    p = pcolormesh(args...; kwargs...)
    margins(0,0)
    p
end

function set_times_new_roman()
    rc("font";
       family = "serif",
       serif = "Times New Roman",
       weight = "regular")
    rc("mathtext";
       rm = "serif",
       it = "serif:italic",
       bf = "serif:bold",
       fontset = "stix")
end

export colormaps, plot_map, set_times_new_roman

end

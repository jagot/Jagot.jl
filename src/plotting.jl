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

function set_font(; kwargs...)
    keys = map(first, kwargs)

    if :serif in keys
        push!(kwargs, (:family, "serif"))
    elseif :sans_serif in keys
        push!(kwargs, (:family, "sans-serif"))
    end

    for i in eachindex(kwargs)
        kwargs[i] = (symbol(replace(string(kwargs[i][1]), "_", "-")), kwargs[i][2])
    end

    rc("font"; kwargs...)
end

function set_times_new_roman(; kwargs...)
    set_font(serif = "Times New Roman", kwargs...)
    rc("mathtext";
       rm = "serif",
       it = "serif:italic",
       bf = "serif:bold",
       fontset = "stix")
end

export colormaps, plot_map, set_times_new_roman
export colormaps, plot_map, set_font, set_times_new_roman

end

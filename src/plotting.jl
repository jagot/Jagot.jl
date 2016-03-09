module plotting
using PyPlot
using PyCall
@pyimport matplotlib.colors as COL

include("colormaps.jl")

function plot_map(args...; kwargs...)
    keys = map(first,kwargs)
    function push_default!(key, val)
        if !(key in keys)
            push!(kwargs, (key, val))
        end
    end
    push_default!(:cmap, get_cmap("viridis"))
    push_default!(:rasterized, true)
    if (i = findfirst(keys, :norm)) != 0 && kwargs[i][2] == :log
        kwargs[i] = (:norm, COL.LogNorm())
    end
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

function latex_base10(v)
    if v==0
        return "0"
    end
    ex = floor(log10(v))
    v /= 10^ex
    @sprintf("%0.3f\\cdot10^{%i}",v,ex)
end

function axis_add_ticks(ticks, labels, ax = :x)
    a = gca()
    a2 = a[ax == :x ? :twiny : :twinx]()
    a2[ax == :x ? :set_xlim : :set_ylim](a[ax == :x ? :get_xlim : :get_ylim]())
    a2[ax == :x ? :set_xticks : :set_yticks](collect(ticks))
    a2[ax == :x ? :set_xticklabels : :set_yticklabels](collect(labels))
    sca(a)
end

export colormaps, plot_map, set_font, set_times_new_roman, latex_base10, axis_add_ticks

end

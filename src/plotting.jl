# * Setup
module plotting
using PyPlot
using PyCall
@pyimport matplotlib.colors as COL
@pyimport numpy.ma as masked_array
@pyimport matplotlib.backends.backend_pgf as pgf_b
using SuperSub

using GSL

# * Colormaps

include("colormaps.jl")

function colorbar_hack(v::AbstractVector,cmap)
    Z = [0 0
         0 0]
    CS3 = contourf(Z, v, cmap=cmap)
    clf()
    CS3
end

# * Map plots

function filter_kwargs(kwargs, sym)
    filter!(f -> f[1] != sym, kwargs)
end

function plot_map(args...; kwargs...)
    keys = map(first,kwargs)
    function push_default!(key, val)
        if !(key in keys)
            push!(kwargs, (key, val))
        end
    end
    push_default!(:cmap, get_cmap("viridis"))
    if (i = findfirst(keys, :cmap)) != 0 && typeof(kwargs[i][2]) <: AbstractString
        kwargs[i] = (:cmap, get_cmap(kwargs[i][2]))
    end
    push_default!(:rasterized, true)
    if (i = findfirst(keys, :norm)) != 0 && kwargs[i][2] == :log
        kwargs[i] = (:norm, COL.LogNorm())
    end
    if (i = findfirst(keys, :align_ticks)) != 0 && (aw=kwargs[i][2]) != false && length(args) >= 3
        kwargs = filter_kwargs(kwargs, :align_ticks)
        x,y,z = args
        if aw != :y
            x -= (x[2]-x[1])/2
        end
        if aw != :x
            y -= (y[2]-y[1])/2
        end
        args = (x,y,z)
    end
    p = pcolormesh(args...; kwargs...)
    margins(0,0)
    p
end

function plot_polar_map(r::AbstractVector, v::AbstractVector, nφ,
                        args...; kwargs...)
    φ = linspace(0,2π,nφ)
    R = repmat(r,1,nφ)
    Φ = repmat(φ',length(r),1)
    X = R.*cos(Φ)
    Y = R.*sin(Φ)
    V = repmat(v,1,nφ)
    plot_map(X,Y,V, args...; kwargs...)
end

function spherical_harmonic_plot(fun::Function,
                                 r::AbstractVector,
                                 J::AbstractVector{Int},
                                 v::AbstractMatrix,
                                 nθ::Int,
                                 args...;
                                 wfn = false,
                                 kwargs...)
    nr = length(r)
    (nr,length(J)) == size(v) || error("Dimension mismatch!")

    θ = linspace(0,2π,nθ)
    cosθ = cos(θ)

    R = repmat(r,1,nθ)
    Θ = repmat(θ',nr,1)

    X = R.*cos(Θ)
    Y = R.*sin(Θ)

    V = zeros(eltype(v), nr, nθ)

    leg_fun = (wfn ? sf_legendre_sphPlm : sf_legendre_Plm)

    for j in eachindex(J)
        PJ = map(x -> leg_fun(J[j], 0, x), cosθ)'
        V += broadcast(*, v[:,j], PJ)
    end

    plot_map(X, Y, fun(V), args...; kwargs...)
end
spherical_harmonic_plot(r::AbstractVector,
                        J::AbstractVector{Int},
                        v::AbstractMatrix,
                        nθ::Int,
                        args...; kwargs...) = spherical_harmonic_plot(identity,
                                                                      r::AbstractVector,
                                                                      J::AbstractVector{Int},
                                                                      v::AbstractMatrix,
                                                                      nθ,
                                                                      args...; kwargs...)

function plot_matrix(a, args...; kwargs...)
    aa = pycall(masked_array.masked_equal, Any, full(a), 0)
    plot_map(aa, args...; kwargs...)
    gca()[:invert_yaxis]()
    square_axs()
end

# * LaTeX/font setup

function set_pgf_to_pdf(preamble=[]; texsystem = "xelatex")
    matplotlib[:rcdefaults]()
    ion()
    matplotlib[:backend_bases][:register_backend]("pdf", pgf_b.FigureCanvasPgf)

    set_latex_serif()
    rc("pgf";
       texsystem = texsystem,
       preamble = preamble)
end

function set_font(; kwargs...)
    keys = map(first, kwargs)

    if :serif in keys
        push!(kwargs, (:family, "serif"))
    elseif :sans_serif in keys
        push!(kwargs, (:family, "sans-serif"))
    end

    for i in eachindex(kwargs)
        kwargs[i] = (Symbol(replace(string(kwargs[i][1]), "_", "-")), kwargs[i][2])
    end

    rc("font"; kwargs...)
end

function set_times_new_roman(; kwargs...)
    set_font(serif = "Times New Roman"; kwargs...)
    rc("mathtext";
       rm = "serif",
       it = "serif:italic",
       bf = "serif:bold",
       fontset = "stix")
end

function set_latex_serif(; kwargs...)
    rc("text";
       usetex = true)
    rc("font";
       family = "serif")
    rc("pgf";
       rcfonts = false)
    rc("mathtext";
       rm = "serif",
       it = "serif:italic",
       bf = "serif:bold",
       fontset = "cm")
    set_font(; kwargs...)
end

function latex(o)
    io = IOBuffer()
    writemime(io, MIME"text/latex"(), o)
    takebuf_string(io)
end

# * Scientific notation

function latex_base10(v)
    if v==0
        return "0"
    end
    ex = floor(log10(abs(v)))
    v /= 10^ex
    if ex == 0
        @sprintf("%0.7f",v)
    elseif ex == 1
        @sprintf("%0.7f\\cdot10",v)
    else
        @sprintf("%0.7f\\cdot10^{%i}",v,ex)
    end
end
latex_base10(v::Vector) = map(latex_base10, v)

function base10(v)
    if v==0
        return "0"
    end
    ex = floor(log10(abs(v)))
    v /= 10^ex
    ex = floor(Int, ex)
    if v==1
        @sprintf("10%s",superscript(ex))
    else
        @sprintf("%0.3f×10%s",v,superscript(ex))
    end
end
base10(v::Vector) = map(base10, v)

# * Axes/tick(label)s

function axis_add_ticks(ticks, labels, ax = :x; kwargs...)
    a = gca()
    a2 = a[ax == :x ? :twiny : :twinx]()
    a2[ax == :x ? :set_xlim : :set_ylim](a[ax == :x ? :get_xlim : :get_ylim]())
    a2[ax == :x ? :set_xscale : :set_yscale](a[ax == :x ? :get_xscale : :get_yscale]())
    a2[ax == :x ? :set_xticks : :set_yticks](vec(collect(ticks)))
    a2[ax == :x ? :set_xticklabels : :set_yticklabels](vec(collect(labels)); kwargs...)
    sca(a)
end

function set_ticklabel_props(ax=:x; kwargs...)
    labels = gca()[ax == :x ? :get_xticklabels : :get_yticklabels]()
    setp(labels; kwargs...)
end

function π_labels(ax = :x, max_count = 10; pi_sym = "\\pi")
    lims = gca()[ax == :x ? :get_xlim : :get_ylim]()
    f = 2
    mi,ma = map(l -> trunc(Int, l/(π/2)), lims)
    if ma-mi > max_count
        f = 1
        mi,ma = map(l -> trunc(Int, l/π), lims)
    end
    d = round(Int, max((ma-mi)/max_count,1))
    r = (mi:d:ma)//f
    gca()[ax == :x ? :set_xticks : :set_yticks](collect(r)*π)

    den_str = v -> den(v) != 1 ? "/$(den(v))" : ""
    tick_labels = map(r) do i
        if i == 0
            L"$0$"
        elseif num(i) == 1
            latexstring("\$$(pi_sym)$(den_str(i))\$")
        elseif num(i) == -1
            latexstring("\$-$(pi_sym)$(den_str(i))\$")
        else
            latexstring("\$$(num(i))$(pi_sym)$(den_str(i))\$")
        end
    end
    gca()[ax == :x ? :set_xticklabels : :set_yticklabels](tick_labels)
end

function frac_ticks(ts::AbstractVector{Rational{Int}}, axis = :x)
    if axis == :x
        xticks(ts)
    else
        yticks(ts)
    end
    tls = map(ts) do t
        if t == 0
            L"0"
        elseif den(t)==1
            latexstring("$(num(t))")
        else
            latexstring("$(num(t))/$(den(t))")
        end
    end
    gca()[axis == :x ? :set_xticklabels : :set_yticklabels](tls)
end

function square_axs()
    gca()[:set_aspect]("equal")
    gca()[:autoscale](tight=true)
end

function axes_labels_opposite(axis = :y, ax = gca())
    if axis == :y
        ax[:yaxis][:tick_right]()
        ax[:yaxis][:set_ticks_position]("both")
        ax[:yaxis][:set_label_position]("right")
    else
        ax[:xaxis][:tick_top]()
        ax[:xaxis][:set_ticks_position]("both")
        ax[:xaxis][:set_label_position]("top")
    end
end

# * Misc

pyslice(args...) = pycall(pybuiltin("slice"), PyObject, args...)

function savefig_f(filename, args...; kwargs...)
    savefig(filename, args...; kwargs...)
    filename
end

# * Exports

export colormaps, colorbar_hack,
plot_map, plot_polar_map, spherical_harmonic_plot, plot_matrix,
set_pgf_to_pdf, set_font, set_times_new_roman, set_latex_serif,
latex, latex_base10, base10,
axis_add_ticks, set_ticklabel_props, π_labels, frac_ticks,
square_axs, axes_labels_opposite,
pyslice, savefig_f

end

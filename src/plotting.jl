# * Setup
module plotting
using PythonPlot
using PythonCall
using Compat

const COL = PythonCall.pynew()
const masked_array = PythonCall.pynew()
const pgf_backend = PythonCall.pynew()
const pyplot_collections = PythonCall.pynew()
const RegularPolyCollection = PythonCall.pynew()
rcParams = PythonCall.pynew()

function __init__()
    PythonCall.pycopy!(COL, pyimport("matplotlib.colors"))
    PythonCall.pycopy!(masked_array, pyimport("numpy.ma"))
    PythonCall.pycopy!(pgf_backend, pyimport("matplotlib.backends.backend_pgf"))
    PythonCall.pycopy!(pyplot_collections, pyimport("matplotlib.collections"))
    PythonCall.pycopy!(RegularPolyCollection, pyplot_collections.RegularPolyCollection)
    global rcParams
    PythonCall.pycopy!(rcParams, pyplot.matplotlib.rcParams)
end

using UnicodeFun

import Jagot: meshgrid

# using GSL
using Printf

plot_style(style::String) = matplotlib.style.use(style)

# * Figure wrappers

function ax_common(;nox=false, noy=false)
    if nox
        no_tick_labels(:x)
        xlabel("")
    end
    if noy
        no_tick_labels(:y)
        ylabel("")
    end
end

function toggle_toolbar(visible)
    rcParams["toolbar"] = visible ? "toolbar2" : "None"
end

toggle_statusbar(visible, fig=gcf()) = hasproperty(fig.canvas, :window) && fig.canvas.window().statusBar().setVisible(visible)

function cfigure(fun::Function, figname;
                 clear=true,
                 nox=false, noy=false,
                 hide_bars=true, tight=true,
                 kwargs...)
    toggle_toolbar(!hide_bars)
    fig = figure(figname; kwargs...)
    clear && clf()
    fun()
    ax_common(nox=nox, noy=noy)
    toggle_statusbar(!hide_bars, fig)
    tight && tight_layout()
    fig
end

function csubplot(fun::Function, ax::Py, args...; kwargs...)
    sca(ax)
    fun()
    ax_common(;kwargs...)
    ax
end

csubplot(fun::Function, args...; projection=nothing, kwargs...) =
    csubplot(fun, subplot(args...; projection=projection), args...; kwargs...)

function csubplot(fun::Function, m::Integer, n::Integer, (i,j)::Tuple{Integer,Integer}, args...; kwargs...)
    LI = LinearIndices((n,m))'
    csubplot(fun, m, n, LI[i,j], args...; kwargs...)
end

function subplot_ratio(N, r; verbosity=0, max_iter=100)
    verbosity > 0 && @info "Generating at least $(N) subplots with approximate width/height ratio of $(r)"
    a = √(N/r)
    b = √(r*N)
    verbosity > 1 && @info "Initial guess" a b a*b b/a
    a = round(Int, a)
    b = round(Int, b)
    verbosity > 1 && @info "Initial guess, rounded" a b a*b b/a
    a*b ≥ N && return a,b
    i = 0
    while a*b < N && i < max_iter
        a, b = if abs((b+1)/a - r) < abs(b/(a+1) - r)
            a, b+1
        else
            a+1, b
        end
        i += 1
        verbosity > 2 && @info "Iteration" a b a*b b/a i
    end
    a, b
end

# ** GridSpec

grid_spec(fig::Figure, args...; kwargs...) =
    fig.add_gridspec(args...; kwargs...)

grid_spec(gr::PythonCall.Py, args...; kwargs...) =
    gr.subgridspec(args...; kwargs...)

gr_sub_plots(gr::PythonCall.Py) = pyconvert(Array, gr.subplots())

clear_figure!(fig::Figure) = fig.clf()
clear_figure!(gr::PythonCall.Py) =
    clear_figure!(gr.get_gridspec().figure)
clear_figure!(::Any) = nothing

tight_layout!(fig::Figure) = fig.tight_layout()
tight_layout!(gr::PythonCall.Py) =
    tight_layout!(gr.get_gridspec().figure)
tight_layout!(::Any) = nothing

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

function filter_kwargs!(kwargs, sym, default=nothing)
    v = get(kwargs, sym, default)
    filter!(kv -> kv[1] != sym, kwargs)
    v
end

# Slightly ugly hack since `size(o::Py)` returns the
# _length_ of `o`.

function shape(pyobj::Py)
    pyhasattr(pyobj, "shape") ||
        throw(ArgumentError("`:shape` not present among `Py`s `attrs`."))
    pyconvert(Tuple, pyobj.shape)
end
shape(pyobj::Py, i) = shape(pyobj)[i]
shape(a, args...) = size(a, args...)

function plot_map(args...; kwargs...)
    kwargs = Dict{Symbol,Any}(kwargs)
    set_default!(key, val) = (kwargs[key] = get(kwargs, key, val))
    kwargs[:cmap] = get_cmap(get(kwargs, :cmap, "viridis"))

    set_default!(:shading, "auto")

    plot_fun = if filter_kwargs!(kwargs, :contour, false)
        contourf
    else
        set_default!(:rasterized, true)
        pcolormesh
    end

    vmin = filter_kwargs!(kwargs, :vmin, nothing)
    vmax = filter_kwargs!(kwargs, :vmax, nothing)
    kwargs[:norm] = (get(kwargs, :norm, :lin) == :lin ?
                     COL.Normalize(vmin, vmax) :
                     COL.LogNorm(vmin, vmax))
    aw = filter_kwargs!(kwargs, :align_ticks, false)
    xtl = filter_kwargs!(kwargs, :xticklabels, nothing)
    ytl = filter_kwargs!(kwargs, :yticklabels, nothing)
    tickrotation = filter_kwargs!(kwargs, :tickrotation, 0)
    xtickrotation = filter_kwargs!(kwargs, :xtickrotation, tickrotation)
    ytickrotation = filter_kwargs!(kwargs, :ytickrotation, tickrotation)
    xticker,yticker = if aw != false
        x,y = if length(args) == 3
            x,y,z = args
            x,y
        elseif length(args) == 1
            z = args[1]
            m,n = shape(z)
            x = 1:n
            y = 1:m
            x,y
        else
            @warn "Don't know how to shift x/y axes in case of $(length(args)) `args`"
            nothing,nothing
        end
        xt = !isnothing(x) && aw != :y ? () -> xticks(x .- (x[2]-x[1])/2,
                                                     isnothing(xtl) ? x : xtl,
                                                     rotation=xtickrotation) : () -> ()
        yt = !isnothing(y) && aw != :x ? () -> yticks(y .- (y[2]-y[1])/2,
                                                     isnothing(ytl) ? y : ytl,
                                                     rotation=ytickrotation) : () -> ()
        xt,yt
    else
        () -> (), () -> ()
    end
    p = plot_fun(args...; kwargs...)
    margins(0,0)
    xticker()
    yticker()
    p
end

function plot_polar_map(r::AbstractVector, ϕ::AbstractVector, V::AbstractMatrix,
                        args...; kwargs...)
    nϕ = length(ϕ)
    R = repeat(r,1,nϕ)
    Φ = repeat(ϕ',length(r),1)
    X = R.*cos.(Φ)
    Y = R.*sin.(Φ)
    plot_map(X,Y,V, args...; kwargs...)
end

plot_polar_map(r::AbstractVector, V::AbstractMatrix, nϕ::Integer, args...; kwargs...) =
    plot_polar_map(r, range(0,stop=2π,length=nϕ), V, args...; kwargs...)

plot_polar_map(r::AbstractVector, v::AbstractVector, nϕ::Integer, args...; kwargs...) =
    plot_polar_map(r, repeat(v,1,nϕ), nϕ, args...; kwargs...)

# function spherical_harmonic_plot(fun::Function,
#                                  r::AbstractVector,
#                                  J::AbstractVector{Int},
#                                  v::AbstractMatrix,
#                                  nθ::Int,
#                                  args...;
#                                  wfn = false,
#                                  J_max = Inf,
#                                  m_min = 0, m_max = Inf,
#                                  kwargs...)
#     nr = length(r)
#     (nr,length(J)) == size(v) || (wfn && (nr,length(J)^2) == size(v)) || error("Dimension mismatch!")

#     θ = range(0,stop=2π,length=nθ)
#     cosθ = cos.(θ)

#     R = repmat(r,1,nθ)
#     Θ = repmat(θ',nr,1)

#     X = R.*cos.(Θ)
#     Y = R.*sin.(Θ)

#     V = zeros(eltype(v), nr, nθ)

#     leg_fun = (wfn ? sf_legendre_sphPlm : sf_legendre_Plm)

#     if size(v,2) == length(J)
#         for j in eachindex(J)
#             J[j] > J_max && break
#             PJ = map(x -> leg_fun(J[j], 0, x), cosθ)'
#             V += broadcast(*, v[:,j], PJ)
#         end
#     elseif !wfn
#         error("3d plot only implemented for wavefunctions")
#     else
#         j = 0
#         for ell in J
#             ell > J_max && break
#             for m in -ell:ell
#                 j += 1
#                 abs(m) > m_max && continue
#                 abs(m) < m_min && continue
#                 Pell_m = (m<0 ? (-1)^abs(m) : 1
#                           )*map(x -> leg_fun(ell, abs(m), x), cosθ)'
#                 V += broadcast(*, v[:,j], Pell_m)
#             end
#         end
#     end

#     plot_map(X, Y, fun.(V), args...; kwargs...)
# end
# spherical_harmonic_plot(r::AbstractVector,
#                         J::AbstractVector{Int},
#                         v::AbstractMatrix,
#                         nθ::Int,
#                         args...; kwargs...) = spherical_harmonic_plot(identity,
#                                                                       r::AbstractVector,
#                                                                       J::AbstractVector{Int},
#                                                                       v::AbstractMatrix,
#                                                                       nθ,
#                                                                       args...; kwargs...)
# * Matrix plots

function plot_matrix(A::AbstractMatrix{T}, args...; align_ticks=true, bias=0.5, kwargs...) where T
    vals = filter(!iszero, unique(A))
    if T <: Complex
        @warn "Plotting real part of complex-valued matrix"
        vals = real(vals)
    end
    vmin,vmax = isempty(vals) ? (1,1) : extrema(vals)
    if vmin == vmax
        δ = 0.1abs(vmin)
        vmin = vmin - 2bias*δ
        vmax = vmax + 2*(1-bias)*δ
    end
    aa = pycall(masked_array.masked_equal, Matrix(A), 0)
    plot_map(aa, args...; align_ticks=align_ticks, vmin=vmin, vmax=vmax, kwargs...)
    gca().invert_yaxis()
    square_axs()
end

# # Inspired by/stolen from
# # https://github.com/tonysyu/mpltools/blob/master/mpltools/special/hinton.py

# @pydef mutable struct SquareCollection <: pyplot_collections.RegularPolyCollection
#     __init__(self; kwargs...) = begin
#         pyplot_collections.RegularPolyCollection.__init__(
#             self, 4, rotation=pi/4; kwargs...
#         )
#     end
#     get_transform(self) = begin
#         """Return transform scaling circle areas to data space."""
#         ax = self.axes
#         pts2pixels = 72.0 / ax.figure.dpi
#         scale_x = pts2pixels * ax.bbox.width / ax.viewLim.width
#         scale_y = pts2pixels * ax.bbox.height / ax.viewLim.height
#         matplotlib.transforms.Affine2D().scale(scale_x, scale_y)
#     end
# end

# function hinton_plot_matrix(a::AbstractMatrix;
#                             max_weight=nothing,
#                             pos_color=nothing,
#                             neg_color=nothing,
#                             bg_color="none",
#                             traditional=false)
#     if traditional
#         pos_color = "white"
#         neg_color = "black"
#         bg_color = "gray"
#     else
#         pos_color == nothing && (pos_color="C0")
#         neg_color == nothing && (neg_color="C1")
#     end

#     ax = gca()
#     ax.set_facecolor(bg_color)

#     height, width = size(a)
#     if max_weight == nothing
#         max_weight = 2^ceil(log2(maximum(abs, a)))
#     end

#     vals = clamp.(a/max_weight, -1, 1)
#     neg = vals .< 0
#     pos = vals .> 0

#     cols,rows = meshgrid(1:width,1:height)

#     for (sel,color) in zip([neg,pos], [neg_color,pos_color])
#         if any(sel)
#             xy = collect(zip(cols[sel], rows[sel]))
#             circle_areas = π/2 * abs.(vals[sel])
#             squares = SquareCollection(sizes=circle_areas,
#                                        offsets=xy, transOffset=ax.transData,
#                                        facecolor=color, edgecolor=color)
#             ax.add_collection(squares, autolim=true)
#         end
#     end
#     square_axs()
#     grid(false)
#     ax.set_xlim(0.5, width+0.5)
#     ax.set_ylim(0.5, height+0.5)
#     gca().invert_yaxis()
#     ax
# end
# hinton_plot_matrix(a;kwargs...) = hinton_plot_matrix(full(a);kwargs...)

# * Patches

function draw_patch(kind, args...; ax=gca(), kwargs...)
    p = getproperty(matplotlib.patches, kind)
    ax.add_patch(p(args...; kwargs...))
end

draw_ellipse((x,y), w, h; kwargs...) = draw_patch(:Ellipse, (x,y), w, h; kwargs...)
draw_circle((x,y), r; kwargs...) = draw_ellipse((x,y), r, r; kwargs...)
draw_arc((x,y), w, h, ϕ, θ₁, θ₂; kwargs...) =
    draw_patch(:Arc, (x,y), w, h, ϕ, θ₁, θ₂; kwargs...)
draw_arc((x,y), r, ϕ, θ₁, θ₂; kwargs...) =
    draw_arc((x,y), r, r, ϕ, θ₁, θ₂; kwargs...)

# * LaTeX/font setup

function set_pgf_to_pdf(preamble=[]; texsystem = "xelatex")
    matplotlib.rcdefaults()
    ion()
    matplotlib.backend_bases.register_backend("pdf", pgf_backend.FigureCanvasPgf)

    set_latex_serif()
    rc("pgf";
       texsystem = texsystem,
       preamble = preamble)
end

function set_font(; kwargs...)
    if :serif in keys(kwargs)
        kwargs[:family] = "serif"
    elseif :sans_serif in keys(kwargs)
        kwargs[:family] = "sans-serif"
    end

    kwargs = Dict(Symbol(replace(string(k), "_" => "-")) => v
                  for (k,v) in kwargs)

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
    show(io, MIME"text/latex"(), o)
    String(take!(io))
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
        @sprintf("10%s",to_superscript(ex))
    else
        @sprintf("%0.3f×10%s",v,to_superscript(ex))
    end
end
base10(v::Vector) = map(base10, v)

# * Axes/tick(label)s

function axis_add_ticks(ticks, labels, ax = :x; kwargs...)
    a = gca()
    a2 = getproperty(a, ax == :x ? :twiny : :twinx)()
    getproperty(a2, ax == :x ? :set_xlim : :set_ylim)(getproperty(a, ax == :x ? :get_xlim : :get_ylim)())
    getproperty(a2, ax == :x ? :set_xscale : :set_yscale)(getproperty(a, ax == :x ? :get_xscale : :get_yscale)())
    getproperty(a2, ax == :x ? :set_xticks : :set_yticks)(vec(collect(ticks)))
    getproperty(a2, ax == :x ? :set_xticklabels : :set_yticklabels)(vec(collect(labels)); kwargs...)
    sca(a)
end

function set_ticklabel_props(ax=:x; kwargs...)
    labels = getproperty(gca(), ax == :x ? :get_xticklabels : :get_yticklabels)()
    setp(labels; kwargs...)
end

function π_frac_string(i; pi_sym = "\\pi")
    den_str = v -> denominator(v) != 1 ? "/$(denominator(v))" : ""
    if i == 0
        L"$0$"
    elseif numerator(i) == 1
        latexstring("\$$(pi_sym)$(den_str(i))\$")
    elseif numerator(i) == -1
        latexstring("\$-$(pi_sym)$(den_str(i))\$")
    else
        latexstring("\$$(numerator(i))$(pi_sym)$(den_str(i))\$")
    end
end

function π_labels(ax = :x; max_count = 10, divisor=4,
                  pi_sym = "\\pi", ca=gca())
    lims = pyconvert(Tuple, getproperty(ca, ax == :x ? :get_xlim : :get_ylim)())
    f = divisor
    mi,ma = map(l -> trunc(Int, l/(π/divisor)), lims)
    if ma-mi > max_count
        f = 1
        mi,ma = map(l -> trunc(Int, l/π), lims)
    end
    d = round(Int, max((ma-mi)/max_count,1))
    r = (mi:d:ma)//f
    getproperty(ca, ax == :x ? :set_xticks : :set_yticks)(collect(r)*π)

    tick_labels = map(r) do i
        π_frac_string(i, pi_sym = pi_sym)
    end
    getproperty(ca, ax == :x ? :set_xticklabels : :set_yticklabels)(tick_labels)
end

function frac_ticks(ts::AbstractVector{Rational{Int}}, axis = :x; sfrac = false)
    fts = float(ts)
    if axis == :x
        xticks(fts)
    else
        yticks(fts)
    end
    tls = map(ts) do t
        if t == 0
            L"0"
        elseif denominator(t)==1
            latexstring("$(numerator(t))")
        else
            if sfrac
                latexstring("\\sfrac{$(numerator(t))}{$(denominator(t))}")
            else
                latexstring("$(numerator(t))/$(denominator(t))")
            end
        end
    end
    getproperty(gca(), axis == :x ? :set_xticklabels : :set_yticklabels)(tls)
end

function sci_ticks(ax, min=0, max=0)
    ticklabel_format(style="sci", axis="$ax", scilimits=(min,max))
end

function colorbar_sci_ticks(cb, min=0, max=0)
    cb.formatter.set_scientific(true)
    cb.formatter.set_powerlimits((min, max))
    cb.update_ticks()
end

function square_axs()
    gca().set_aspect("equal")
    gca().autoscale(tight=true)
end

function axes_labels_opposite(axis = :y, ax = gca())
    if axis == :y
        ax.yaxis.tick_right()
        ax.yaxis.set_label_position("right")
    else
        ax.xaxis.tick_top()
        ax.xaxis.set_label_position("top")
    end
end

function no_tick_labels(axis = :x, ax = gca(), ticks = false)
    a,b,la,lb = if axis == :x
        :bottom,:top,:labelbottom,:labeltop
    else
        :left,:right,:labelleft,:labelright
    end
    ax.tick_params(; which="both",
                   Dict(a => ticks, b => ticks,
                        la => false, lb => false)...)
end

# ** Sqrt axes

# mticker = matplotlib.ticker

# SqrtTransform = pytype("SqrtTransform", (matplotlib.transforms.Transform,), [
#     "__module__" => "__main__",
#     pyfunc(name="__init__",
#            function (self)
#                self.input_dims = 1
#                self.output_dims = 1
#                self.is_separable = true
#                self.has_inverse = true
#                matplotlib.transforms.Transform.__init__(self)
#            end),
#     pyfunc(name="transform_non_affine",
#            function(self, a)
#                b = pyconvert(Array, a)
#                v = similar(b)
#                sel = b .>= 0
#                v[sel] .= .√(b[sel])
#                # v[.!sel] .= NaN
#                v
#            end),
#     pyfunc(name="inverted",
#            (self) -> SquareTransform())
# ])

# SquareTransform = pytype("SquareTransform", (matplotlib.transforms.Transform,), [
#     "__module__" => "__main__",
#     pyfunc(name="__init__",
#            function (self)
#                self.input_dims = 1
#                self.output_dims = 1
#                self.is_separable = true
#                self.has_inverse = true
#                matplotlib.transforms.Transform.__init__(self)
#            end),
#     pyfunc(name="transform_non_affine",
#            function(self, a)
#                pyconvert(Array, a) .^ 2
#            end),
#     pyfunc(name="inverted",
#            (self) -> SqrtTransform())
# ])

# SqrtScale = pytype("SqrtScale", (matplotlib.scale.ScaleBase,), [
#     "__module__" => "__main__",
#     "name" => "sqrt",
#     pyfunc(name="__init__",
#            function(self, axis, args...; kwargs...)
#                matplotlib.scale.ScaleBase.__init__(self, axis)
#            end),
#     pyfunc(name="get_transform",
#            (self) -> SqrtTransform()),
#     pyfunc(name="set_default_locators_and_formatters",
#            (self, axis) -> nothing),
#     pyfunc(name="limit_range_for_scale",
#            function(self, vmin, vmax, minpos)
#                max(vmin, 0), max(vmax, 0)
#            end)])

# # https://github.com/cjdoris/PythonCall.jl/issues/289
# matplotlib.scale.register_scale(SqrtScale)

# * Misc

# pyslice(args...) = pycall(pybuiltin("slice"), Py, args...)

function savefig_f(filename, args...; kwargs...)
    savefig(filename, args...; kwargs...)
    filename
end

reltext(x,y,string,args...;kwargs...) = text(x,y,string,transform=gca().transAxes,args...;kwargs...)

function disp()
    display(gcf())
    println()
end

# GridSpec = matplotlib.gridspec.GridSpec

# https://stackoverflow.com/a/51207905/1079038
function next_color(ax=gca())
    pc = ax._get_lines.prop_cycler
    first(pc)["color"]
end

neutral_color() = matplotlib.rcParams["axes.edgecolor"]

for (f,ff) in ((:hline,:axhline), (:vline,:axvline))
    @eval begin
        $f(args...; kwargs...) =
            $ff(args...; color=neutral_color(), linewidth=1.0, linestyle="--", kwargs...)
    end
end

# * ICC support
include("save_pgf_with_icc.jl")

# * PythonPlot recipes
include("python_plot_recipes.jl")

# * Exports

export plot_style,
    toggle_toolbar, toggle_statusbar,
    cfigure, csubplot, subplot_ratio,
    grid_spec, gr_sub_plots, clear_figure!, tight_layout!,
    colormaps, colorbar_hack,
    plot_map, plot_polar_map, spherical_harmonic_plot, plot_matrix, hinton_plot_matrix,
    draw_patch, draw_ellipse, draw_circle, draw_arc,
    set_pgf_to_pdf, set_font, set_times_new_roman, set_latex_serif,
    latex, latex_base10, base10,
    axis_add_ticks, set_ticklabel_props, π_frac_string, π_labels, frac_ticks, sci_ticks, colorbar_sci_ticks,
    square_axs, axes_labels_opposite, no_tick_labels,
    # pyslice,
    savefig_f, reltext, disp,
    # GridSpec, 
    next_color, hline, vline

end

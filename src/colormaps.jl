# New matplotlib colormaps by Nathaniel J. Smith, Stefan van der Walt,
# and (in the case of viridis) Eric Firing.
#
# This file and the colormaps in it are released under the CC0 license /
# public domain dedication. We would appreciate credit if you use or
# redistribute these colormaps, but do not impose any legal restrictions.
#
# To the extent possible under law, the persons who associated CC0 with
# mpl-colormaps have waived all copyright and related or neighboring rights
# to mpl-colormaps.
#
# You should have received a copy of the CC0 legalcode along with this
# work.  If not, see <http://creativecommons.org/publicdomain/zero/1.0/>.
#
# Adapted for Julia by Jagot

module colormaps
using PythonCall
using PythonPlot

const COL = PythonCall.pynew()

function __init__()
    PythonCall.pycopy!(COL, pyimport("matplotlib.colors"))
end

function (cmap::PythonPlot.ColorMap)(i::Int)
    if :colors in propertynames(cmap)
        pyconvert(Vector, cmap.colors[i])
    else
        pycall(cmap, PyAny, i) # This is probably wrong
    end
end

function lerp(a,b,t)
    (1-t)*a + t*b
end

function lerp(a::Tuple, b::Tuple, t)
    tuple([lerp(a[i],b[i],t) for i = 1:length(a)]...)
end

function (cmap::PythonPlot.ColorMap)(f::Real)
    N = pyconvert(Int, cmap.N)
    i = clamp(f*(N-1),0,N-1)
    fl,ce = floor(Int,i), ceil(Int, i)
    lerp(cmap(fl),cmap(ce),f)
end

(cmap::PythonPlot.ColorMap)(i::Integer, l::Integer) =
    l > 1 ? cmap((i-1)/(l-1)) : cmap(0.5)

end

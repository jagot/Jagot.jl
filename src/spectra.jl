module spectra
using FFTW

get_freq(t, fs) = fftshift(fftfreq(round(Int,fs*(maximum(t)-minimum(t)))+1, fs))
get_freq(t) = get_freq(t, 1.0/(t[2]-t[1]))

get_spectrum(v::Vector) = fftshift(fft((v.*hanning(length(v)))))
get_spectrum(v::Matrix) = fftshift(fft((v'.*hanning(size(v,2)))',2),2)
function get_spectrum(v::Matrix, tlim, ndt)
    j = tlim*ndt
    j = j[1]:j[2]
    get_spectrum(v[:,j])
end
function get_spectrum(v::Vector, tlim, ndt)
    j = tlim*ndt
    j = j[1]:j[2]
    get_spectrum(v[j])
end

export get_freq, get_spectrum
end

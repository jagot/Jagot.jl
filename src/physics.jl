module physics
using Physical

I_au(I_si) = I_si/(3.5094452e16*Watt/(Centi*Meter)^2)

export I_au
end

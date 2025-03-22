#.....SYNCHROTRON EMISSION WITH APPROXIMATED KERNEL 
#.....BASED ON FOUKA & OUICHAOUI, https://arxiv.org/abs/1301.6908
#MODULES
@everywhere using Base
@everywhere using SpecialFunctions

@everywhere function fpxeta(xx, plocal, a1, a2, a3, Cp, kp, b1)
  if b1 >= -1e-2
    b1 = -1e-2
  end
  fp = kp * xx^0.333 * exp(a1 * xx^2 + a2 * xx + a3 * xx^0.666) + Cp * xx^(-0.5 * (plocal - 1)) * (1 - exp(b1 * xx^2))^(plocal * 0.2 + 0.5)
  return fp
end

@everywhere function ap(plocal)
  ap = -0.033 - 0.104 * plocal + 0.115 * plocal^2
  return ap
end


@everywhere function fouka_kernel(xx, etag, plocal)
  fp_fit = 0.0

  a1 = -0.14602 + 3.62307 * 1e-2 * plocal - 5.76507 * 1e-3 * plocal^2 + 3.46926 * 1e-4 * plocal^3

  a2 = -0.366648 + 0.18031 * plocal - 7.30773 * 1e-2 * plocal^2 + 1.12484 * 1e-2 * plocal^3 - 6.17683 * 1e-4 * plocal^4

  a3 = 9.69376 * 1e-2 - 0.48892 * plocal + 0.14024 * plocal^2 - 1.93678 * 1e-2 * plocal^3 + 1.01582 * 1e-3 * plocal^4

  b1 = -0.20250 + 5.43462 * 1e-2 * plocal - 8.44171 * 1e-3 * plocal^2 + 5.21281 * 1e-4 * plocal^3

  kp = (pi * 2^(8 / 3)) / (sqrt(3) * (plocal - 0.333) * gamma(0.3333))

  Cp = 2^((plocal + 1) * 0.5) / (plocal + 1) * gamma(plocal * 0.25 + 19 / (12)) * gamma(plocal * 0.25 - 1 / (12))

  xc = (2.028 - 1.187 * plocal + 0.240 * plocal^2) * etag^2
  if xx < xc
    fp_fit = fpxeta(xx, plocal, a1, a2, a3, Cp, kp, b1) + etag^(-plocal + 1) * fpxeta(xx / (etag^2), plocal, a1, a2, a3, Cp, kp, b1)
  end
  if xx >= xc
    fp_fit = sqrt(pi * 0.5) * etag^(-plocal + 2) * xx^(-0.5) * exp(-xx / (etag^2)) * (1 + ap(plocal) * etag^2 / (xx))
  end
  
  return fp_fit

end



@everywhere function synchrotronK(pval::Vector{Float64}, b0::Float64, pe::Vector{Float64}, freq::Float64)

  local g2 = np - 1
  local g1 = 1

  @inbounds for gg in eachindex(pval)

    if pe[gg] <= 1 && pe[gg+1] >= 1
      g1 = gg
    end
    if pe[gg] >= 1 && pe[gg+1] <= 1
      g2 = gg
    end
  end
  #.,,,,compute syncrothron emissivity

  local integ = 0.0
  local freqc = 0.0
  local emission = 0.0

  vL = qe * b0 / (2 * pi * me * vc)

  emission = 1.0
  totg = 0.0
  @inbounds for gg in g1:g2 #...momenta

    @fastmath gam1 = sqrt(1 + (10^pval[gg])^2)
    @fastmath gam2 = sqrt(1 + (10^pval[gg+1])^2)
    @fastmath v1 = 1.5 * gam1^2 * vL
    @fastmath eta = gam2 / gam1
    @fastmath degam = gam2 - gam1
    xx = freq / v1
    plocal = 0.1

    if pe[gg] > pe[gg+1] && pe[gg+1] >= 1

      @fastmath plocal = (log10(pe[gg] / 10^pval[gg]) - log10(pe[gg+1] / 10^pval[gg+1])) / (log10(gam2) - log10(gam1))
      if plocal <= 1.001
        plocal = 1.001
      end
 
      if isfinite(plocal) == false || plocal > 6.0
        continue
      end

      Cnorm = pe[gg] * (plocal - 1) / (gam1^(-plocal + 1) - gam2^(-plocal + 1))
      P1 = pi * sqrt(3) * qe^2 * vL * gam1^(-plocal + 1) * Cnorm / (vc)
      fp = fouka_kernel(xx, eta, plocal)
      psync = P1 * fp

      if isfinite(psync)
        emission += psync
      end
    end
  end

  return emission
end


@everywhere function bessel_integ1(vv::Float64)
  totf = 0.0
  dx = 0.01 #0.2e-3
  ymin = convert(Int64, trunc(vv / dx))
  #    totf=Array{Float64}(1)
  #    totf=Array{Float64}(1)
  #  totf[1]=0.

  @simd for y in ymin:1e2  #....K53 modified bessel function
    @fastmath argb = y * dx#(1.+y*dx
    #if argb >= 30
    @fastmath bbes = besselk(1.666, argb)

    totf += (vv * bbes * dx)
    # end
  end

  return totf
end


@everywhere function bessel_integ3(vv::Float64)
  #   ak1=[-0.97947838884478688,-0.83333239129525072,0.15541796026816246]
  #   ak2=[-4.69247165562628882e-2,-0.70055018056462881,1.03876297841949544e-2]
  a1 = -0.14602
  a2 = -0.36648
  a3 = 9.69376e-2
  b1 = -0.20250
  p = 3.0
  F1 = π * 2.0^1.666 / (sqrt(3.0) * gamma(0.333))
  κp = 2.0 * F1 / (p - 0.333)
  Cp = 2.0 * ((p + 1.0) * 0.5) / (p + 1.0) * gamma(p * 0.25 + 19.0 / 12.0) * gamma(p * 0.25 - 1.0 / 12.0)
  ffn = κp * vv^0.333 / (1.0 + (κp / Cp) * vv^((p - 0.333) * 0.5))

  return ffn

end

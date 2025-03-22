@everywhere function plot_spectra(pe_tot, pe_tot2, pe_tot3, pval)

   plot_nmin = 1e40
   plot_nmax = 1e70
   lfs = 14
   xfs = 12
   xtfs = 12
   pshow = similar(pval)
   for i in eachindex(pval)
      pshow[i] = 10.0^pval[i]
   end
   model1 = "model 1"
   model2 = "model 2"
   model3 = "model 3"
   plot(pshow, pe_tot[:], label=model1, color="black", line=:solid, linewidth=3, dpi=500, alpha=0.7, grid=false, legendfontsize=lfs, yguidefontsize=xfs, xguidefontsize=xfs, xtickfontsize=xtfs, ytickfontsize=xtfs)
   plot!(pshow, pe_tot2[:], color="green", label=model2, line=:dash, linewidth=2, alpha=0.6)
   plot!(pshow, pe_tot3[:], color="red", label=model3, line=:solid, linewidth=3, alpha=0.5)
   yaxis!(L"pN(p)", :log10, (plot_nmin, plot_nmax), fonts=20)
   xaxis!(L"P/(m_ec)", :log10, (10^p_min, 10^p_max), fonts=20, xticks=[1, 10, 100, 1000, 10000, 1e5, 1e6])
   title!(string("z=", zz), fonts=20)
   filep1 = string(root_out, "/SPECTRA/spectra_new_momenta_radial", snapn, "_initial", test, ".png")

   savefig(filep1)

end

@everywhere function plot_spectra_ini(pe_tot, pe_tot2, pe_tot3, pval)

   plot_nmin = 1e40
   plot_nmax = 1e70
   lfs = 14
   xfs = 12
   xtfs = 12
   pshow = similar(pval)
   for i in eachindex(pval)
      pshow[i] = 10.0^pval[i]
   end
   model1 = "model 1"
   model2 = "model 2"
   model3 = "model 3"
   plot(pshow, pe_tot[:], label=model1, color="black", line=:solid, linewidth=2.5, dpi=500, alpha=0.7, grid=false, legendfontsize=lfs, yguidefontsize=xfs, xguidefontsize=xfs, xtickfontsize=xtfs, ytickfontsize=xtfs)
   plot!(pshow, pe_tot2[:], color="green", label=model2, line=:dash, linewidth=2, alpha=0.6)
   plot!(pshow, pe_tot3[:], color="red", label=model3, line=:solid, linewidth=1.5, alpha=0.5)
   yaxis!(L"pN(p)", :log10, (plot_nmin, plot_nmax), fonts=20)
   xaxis!(L"P/(m_ec)", :log10, (1.1 * 10^p_min, 10^p_max), fonts=20)
   title!(string("z=", zz, "input spectra"), fonts=20)
   filep1 = string(root_out, "/SPECTRA/spectra_new_momenta_radial", snapn, "_initial", test, ".png")
   savefig(filep1)

end


@everywhere function total_spectra(pval::Array{Float64}, pe::Array{Float64}, pe_tot::Array{Float64}, nn1::Int64, nn2::Int64, x::Array{Float64}, y::Array{Float64}, z::Array{Float64}, rspec::Float64, xc::Float64)

   @inbounds @simd for i in nn1:nn2
      @inbounds @simd for g in 1:np
         if pe[g, i] > 0.0 && pe[g, i] < 1e70
            pe_tot[g] += (pe[g, i] * 10^pval[g])
         end
      end
   end

   return pe_tot
end



@everywhere function plot_radio_spectra(total_radio1, total_radio2, total_radio3, freqf)
   plot_nmin = 1e40
   plot_nmax = 1e70
   lfs = 14
   xfs = 12
   xtfs = 12
   model1 = "model 1"
   model2 = "model 2"
   model3 = "model 3"

   plot(freqf, total_radio1, label=model1, color="black", line=:solid, dpi=500, linewidth=2.5, alpha=0.7, grid=false, legendfontsize=lfs, yguidefontsize=xfs, xguidefontsize=xfs, xtickfontsize=xtfs, ytickfontsize=xtfs)
   plot!(freqf, total_radio2, color="green", label=model2, line=:dash, linewidth=2, alpha=0.6)
   plot!(freqf, total_radio3, color="red", label=model3, line=:solid, linewidth=1.5, alpha=0.5)

   #....VERY APPROXIMATE estimate of what can be detected a low redshift
   visib = similar(freqf)
   for i in eachindex(freqf)
      visib[i] = 1e27 * (freqf[1] / freqf[i])
   end
   plot!(freqf, visib, color="grey", label="detectable emission", alpha=0.5, linewidth=1)
   title!(string("z=", zz), fonts=20)
   yaxis!(L"[erg/s/Hz]", :log10, (1e15, 1e35), fonts=20)
   xaxis!(L"\nu[Hz]", :log10, (freqf[1] * 0.9, freqf[nfreq[1]] * 1.1), fonts=20)
   model = ""

   filep1 = string(root_out, "/RADIO/RADIO_spectra_new_momenta_radial", snapn, "_", test, ".png")
   savefig(filep1)

end


@everywhere function map2d(pradio::Array{Float64,2}, x::Array{Float64}, y::Array{Float64}, z::Array{Float64}, nn1::Int64, nn2::Int64, freqf::Array{Float64}, map::Array{Float64,3}, xc::Float64, pee::Array{Float64,2}, dezoom::Float64)#,total_radio::Array{Float64},rspec::Float64)
   xc = 1
   #....maps of all emission frequencies + maps with all projected number of CRs in the tracers (last level of the fits map)

   @inbounds @simd for i in nn1:nn2

      @fastmath x1 = convert(Int32, trunc(dezoom * xc * x[i]))  #...
      @fastmath y1 = convert(Int32, trunc(dezoom * xc * y[i]))
      @fastmath z1 = convert(Int32, trunc(dezoom * xc * z[i]))
      im = findall(x -> (x > 0), pee[:, i])
      @inbounds map[x1, z1, 1] += pradio[i, 1]
      @inbounds map[x1, z1, 2] += pradio[i, 2]
      @inbounds map[x1, z1, 3] += pradio[i, 3]
      @inbounds map[x1, z1, 4] += pradio[i, 4]
      @inbounds map[x1, z1, 5] += pradio[i, 5]
      @inbounds map[x1, z1, 6] += sum(pee[im, i])

   end

   return map
end


@everywhere function chunk(arr, nw)
   n = Int(ceil(length(arr) / nw))
   chunk_out = [arr[i:min(i + n - 1, end)] for i in 1:n:length(arr)]
   return chunk_out
end

@everywhere function write_maps_fits(map, model)
   filep1 = string(root_out, "/MAPS/", model, "/map_radio_tracers", runm, "_", snapn, model, "Y_momenta_", test, ".fits")
   @inbounds for i in eachindex(map)
      map[i] = log10(map[i] + 1)
   end
   f = FITS(filep1, "w")
   write(f, map)
   close(f)
end


@everywhere function do_contours(map1, map2, map3, so)

   filep1 = string(root_out, "/PNG/", "contours_", snapn, "_", test, "_CRe3.png")
   filep2 = string(root_out, "/PNG/", "contours_", snapn, "_", test, "_33.png")
   filed = string(root_map, "/fits/Dmovie_4C", snapn, ".fits")   #.....existing fits maps with projected denisty, X-ray etc... (must be generated externally of ROGER, with the simulation data)

   f1 = FITS(string(filed), "r")
   dens = read(f1[1])
   close(f1)

   dmap = log10.(dens[:, :])

   #...rescale for plots
   dmin = minimum(dmap)
   dmax = maximum(dmap)
   @inbounds for i in eachindex(dmap)
      dmap[i] = (dmap[i] - dmin) / (dmax - dmin)
   end
   i1 = 1
   i2 = n
   j1 = 1
   j2 = n

   heatmap(60 * dmap[i1:i2, j1:j2], colorscale="Plasma", aspect_ratio=:equal, xlims=(1, i2 - i1), ylims=(1, j2 - j1), dpi=300, colorbar=true, title=string("z=", zz), xformatter=_ -> "", yformatter=_ -> "") #xtickfontcolor=:white,ytickfontcolor=:white)  
   contour!(map3[i1:i2, j1:j2, 6], seriescolor=:greens, levels=[55, 55.5, 56, 56.5, 57, 57.5, 58, 58.5, 59, 59.5, 60, 60.5, 61.5, 62, 62.5, 63, 63.5])    #levels=[24,25,26,27,28,29])
   savefig(filep1)

   heatmap(25 * dmap[i1:i2, j1:j2], colorscale="Plasma", aspect_ratio=:equal, xlims=(1, i2 - i1), ylims=(1, j2 - j1), dpi=300, colorbar=true, title=string("z=", zz), xformatter=_ -> "", yformatter=_ -> "") #xtickfontcolor=:white,ytickfontcolor=:white)
   contour!(map3[i1:i2, j1:j2, 2], seriescolor=:reds, levels=[20, 21, 22, 23, 24, 25, 26, 27, 28, 29])
   savefig(filep2)


end



function write_hdf5(pe, pe2, pe3, pval, tr, sel, pradio1, pradio2, pradio3, z, freqf, nn1, nn2)

   filep1 = string(root_out, "_radio_tracer_", runm, "_", snapn, "_pradio_", test, "_3D.hdf5")
   h5write(filep1, "N(p)_model1", pe[:, :])
   h5write(filep1, "N(p)_model2", pe2[:, :])
   h5write(filep1, "N(p)_model3", pe3[:, :])
   h5write(filep1, "P", pval[:])
   h5write(filep1, "Density", tr[7, sel[nn1:nn2]])
   h5write(filep1, "Temperature", tr[8, sel[nn1:nn2]])
   h5write(filep1, "mach", tr[9, sel[nn1:nn2]])
   h5write(filep1, "div_v", tr[11, sel[nn1:nn2]])
   h5write(filep1, "curl_v", tr[12, sel[nn1:nn2]])
   h5write(filep1, "Bfield", tr[10, sel[nn1:nn2]])
   h5write(filep1, "Radio_Power_A", pradio1[:, :])
   h5write(filep1, "Radio_Power_B", pradio2[:, :])
   h5write(filep1, "Radio_Power", pradio3[:, :])
   h5write(filep1, "Redshift", z)
   h5write(filep1, "Freqf", freqf[:])
   h5write(filep1, "x", tr[1, sel[nn1:nn2]])
   h5write(filep1, "y", tr[2, sel[nn1:nn2]])
   h5write(filep1, "z", tr[3, sel[nn1:nn2]])


end

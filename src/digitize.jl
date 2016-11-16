using StatsBase

"""
    digitize(data, cal0, cal1, mode=:equal)

  `data`  --  iterable holding floating point values to be digitized

  `cal0`  --  vectors of floating point values corresponding to samples of
  `cal1`      voltages from measurements of the 0 and 1 states

  `mode`  --  mode of operation for thresholding. Supportted options are `:equal` and `:max`

  Translate floating point votlates in `data` to binary outcomes based on
calibration samples for 0 and 1 outcomes. It is assumed that the number of samples
in `cal0` and `cal1` is large enough that the lower and upper quartile estimate are
accurate.

  The `:equal` mode of operation tries to equalize the conditional error probabilities for the
two possible outcomes, and it is generally less sensitive to statistical fluctuations. The
`:max` mode of operation tries to minimize the overall error probability, but is is much more
sensitive to statistical fluctuations, and is therefore deprecated.

"""
function digitize(data, cal0::Vector{Float64}, cal1::Vector{Float64}; mode = :equal)
    minA, maxA = extrema(vcat(quantile(cal0,[.25,.75]),quantile(cal1,[.25,.75])))
    Arange = linspace(minA,maxA,1000)
    thr = if mode==:max
            Arange[indmax(abs(ecdf(cal0)(Arange)-ecdf(cal1)(Arange)))] 
          elseif mode==:equal
            Arange[indmin(abs(ecdf(cal0)(Arange)-1+ecdf(cal1)(Arange)))]
          else
            error("Unsupported digitization mode")
          end
    
    ddata = data .< thr

    if sum(cal0 .> thr)/length(cal0) > 1/2
        ddata = ~ddata
    end
        
    return round(Int,ddata)
end

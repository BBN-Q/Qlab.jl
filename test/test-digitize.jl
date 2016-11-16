using Base.Test

function test_digitize(sep; C=1000, K=1000)
    cal0 = randn(C)
    cal1 = randn(C)+sep

    # flip 1000 coins, and add noise
    rbits  = round(Int, rand(K) .> 0.5)
    rvolts = sep*rbits + randn(K)

    # estimate of prob of being correct
    pc_est_m = sum(rbits .== digitize(rvolts, cal0, cal1, mode=:max))/K
    pc_est_e = sum(rbits .== digitize(rvolts, cal0, cal1, mode=:equal))/K
    # compute theoretical prob of being correct
    pc     = 1/2 - 1/2*erf(sep/(2*sqrt(2)))
    
    return pc, pc_est_m, pc_est_e
end

@test all(floor(mean(reduce(hcat,[[test_digitize(0.75)...] for _ in 1:1000]),2)*100)/100 .== 0.35)
@test all(floor(mean(reduce(hcat,[[test_digitize(0.5,C=10_000,K=10_000)...] for _ in 1:1000]),2)*100)/100 .== 0.40)

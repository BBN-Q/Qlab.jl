using Base.Test

function test_digitize(sep; C=1000, K=1000)
    cal0 = randn(C)
    cal1 = randn(C)+sep

    # flip 1000 coins, and add noise
    rbits  = round(Int, rand(K) .> 0.5)
    rvolts = sep*rbits + randn(K)

    # estimate of prob of being correct
    dm, pe0m, pe1m = digitize(rvolts, cal0, cal1, mode=:max)
    de, pe0e, pe1e = digitize(rvolts, cal0, cal1, mode=:equal)
    # compute theoretical prob of being correct
    pc = 1/2 - 1/2*erf(sep/(2*sqrt(2)))
    
    return pc, 
           sum(abs(dm-rbits))/length(rbits), 
           sum(abs(de-rbits))/length(rbits), 
           pe0m, 
           pe1m, 
           pe0e, 
           pe1e 

end

@test all(map(x->isapprox(x,0.35,atol=5e-2),mean(reduce(hcat,[[test_digitize(0.75)...] for _ in 1:1000]),2)))
@test all(map(x->isapprox(x,0.40,atol=5e-2),mean(reduce(hcat,[[test_digitize(0.50)...] for _ in 1:1000]),2)))

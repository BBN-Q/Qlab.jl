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
    pe = 1/2 - 1/2*erf(sep/(2*sqrt(2)))
    
    return pe, 
           sum(abs(dm-rbits))/length(rbits), 
           sum(abs(de-rbits))/length(rbits), 
           pe0m, 
           pe1m, 
           pe0e, 
           pe1e 

end

# run a bunch of tests, store each as a column of a matrix, take the mean
# the rows are: 
#   - predicted error rate (asymptotic)
#   - observed error rate for :max mode
#   - observed error rate for :equal mode
#   - estimated error rate based on cals for :max mode (0 preps)
#   - estimated error rate based on cals for :max mode (1 preps)
#   - estimated error rate based on cals for :equal mode (0 preps)
#   - estimated error rate based on cals for :equal mode (1 preps)
test1 = mean(reduce(hcat,[[test_digitize(0.75)...] for _ in 1:1000]),2)
test2 = mean(reduce(hcat,[[test_digitize(0.50)...] for _ in 1:1000]),2)
# make sure these are all close to the prediction
@test all(map(x->isapprox(x,0.35,atol=5e-2),test1))
@test all(map(x->isapprox(x,0.40,atol=5e-2),test2))

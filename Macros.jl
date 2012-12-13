module Macros
using Base

export @timeit

macro timeit(ex)
    quote
        # run once to make sure it passes through JIT
        $(esc(ex))
        singleRunTime = @elapsed $(esc(ex))
        # repeat enough times such that it takes about a second, but at least 3 times
        numRepeats = convert(Int32, max(3, floor(1.0/singleRunTime)))
        times = zeros(numRepeats)
        for ct = 1:numRepeats
            times[ct] = @elapsed $(esc(ex))
        end
        @printf("Mean: %fs, Best: %fs (%d runs)\n", mean(times), min(times), numRepeats)
        mean(times)
    end
end

end
# Copyright 2017 Raytheon BBN Technologies
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the license at
#
#		http://www.apache.org/licenses/LICENSE-2.0
using Test

function in_range(x, minx, maxx)
  return (x >= minx) && (x <= maxx)
end

function check_fit_result(true_values::Dict{String, Float64}, result::Qlab.FitResult)
    @test Set(keys(true_values)) == Set(keys(result.fit_params))
    @test Set(keys(true_values)) == Set(keys(result.errors))

    for k in keys(true_values)
        println("$(k): $(result.fit_params[k]) => $(true_values[k]) ± $(result.errors[k])")
        @test in_range(result.fit_params[k], true_values[k]-10*result.errors[k], true_values[k]+10*result.errors[k])
    end
end

#y = a exp(-t/T₁) + b
function test_T1_fit()
    #generate some fake data
    vals = Dict("a" => 0.87, "T" => 20., "b" => 0.1)
    t = collect(range(0., stop=100., length=201))
    y = vals["a"] * exp.(-t/vals["T"]) .+ vals["b"] .+ 0.01 * randn((201,))
    result = Qlab.fit_t1(t, y)
    check_fit_result(vals, result)
end

#y = ax + b.
function test_line_fit()
    vals = Dict("a" => -0.4, "b" => 2.1)
    x = range(-4, stop=4, length=201)
    y = vals["a"] * x .+ vals["b"] .+ 0.05 * randn((201,))
    result = Qlab.fit_line(x,y)
    check_fit_result(vals, result)
end

#a*exp(-t ./ T).*cos(2πf .* t + ϕ) + b
function test_ramsey_fit()
    p = Dict("a"=>0.8, "T"=>12.3, "f"=>0.55, "ϕ"=>π/12, "b"=>-0.15)
    t = range(0, stop=25, length=201)
    y = p["a"]*exp.(-t ./ p["T"]).*cos.(2π*p["f"] .* t .+ p["ϕ"]) .+ p["b"] + 0.01 * randn((201,))
    result = Qlab.fit_ramsey(t,y)
    check_fit_result(p, result)
end

test_T1_fit()
test_line_fit()
test_ramsey_fit()

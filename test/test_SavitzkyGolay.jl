using Test
using Qlab: savitzky_golay_filter

# test the innput validity checks
# polynomial order greater than window size
@test_throws AssertionError savitzky_golay_filter(1:10, 5, 5)
# window size must be odd
@test_throws AssertionError savitzky_golay_filter(1:10, 10, 2)
# boundary_mode must be one :nearest or :interpolation
@test_throws AssertionError savitzky_golay_filter(1:10, 11, 2; boundary_mode=:silly)  

# a window size of 1 with a zeroth order polynomial should return the data
@test savitzky_golay_filter([1.2], 1, 0) == [1.2]
@test savitzky_golay_filter([1.2, 3.4, 5.6], 1, 0) == [1.2, 3.4, 5.6]

# a line smoothed with a polynomial order of 1 should return itself or `boundary_mode==:interpolation`
@test savitzky_golay_filter(1:10, 3, 1; boundary_mode=:interpolation) ≈ 1:10

# take a low order polynomial and smooth with a higher order polynomial to return the exact value
xpts = range(1,10; length=81)
data = 1.2*xpts.^3 + 3.4*xpts.^2 + 5.6*xpts .+ 7.8

@test savitzky_golay_filter(data, 5, 4) ≈ data
@test savitzky_golay_filter(data, 11, 4) ≈ data

# derivatives are easy to analytically calculate
# we need to scale savitzky_golay_filter by the spacing
deriv1 = 3*1.2*xpts.^2 + 2*3.4*xpts .+ 5.6
@test savitzky_golay_filter(data, 5, 4; deriv_order=1) / (xpts[2] - xpts[1]) ≈ deriv1
deriv2 = 2*3*1.2*xpts .+ 2*3.4
@test savitzky_golay_filter(data, 5, 4; deriv_order=2) / (xpts[2] - xpts[1])^2 ≈ deriv2

# test the :nearest boundary condition
# take 2 points and use a 3 point constant smoothing fit
# first point will be average of [1, 1, 2] = 4/3 and second point will be of [1, 2, 2] = 5/3
@test savitzky_golay_filter([1, 2], 3, 0; boundary_mode=:nearest) ≈ [4/3, 5/3]

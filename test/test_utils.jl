using Test,
      RandomQuantum,
      QuantumInfo,
      Cliffords,
      Qlab

import SchattenNorms.trnorm

function do_unitary2choi()
    choi = unitary2choi(str2unitary("1QX90p"))
    ideal = [0.25+0.0im    0.0+0.25im   0.0+0.25im  0.25+0.0im;
             0.0-0.25im  0.25+0.0im   0.25+0.0im    0.0-0.25im;
             0.0-0.25im  0.25+0.0im   0.25+0.0im    0.0-0.25im;
             0.25+0.0im    0.0+0.25im   0.0+0.25im  0.25+0.0im]
    @test trnorm(choi-ideal) < 3e-2
end

function do_choi2pauliMap()
    choi = unitary2choi(str2unitary("1QX90p"))
    p_map = choi2pauliMap(choi)
    ideal = [1.0 0.0 0.0 0.0;
             0.0 1.0 0.0 0.0;
             0.0 0.0 0.0 1.0;
             0.0 0.0 -1.0 0.0]
    @test isapprox(choi, ideal, atol=1e-14)
end

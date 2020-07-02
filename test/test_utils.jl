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

function do_choi2chi()
    choi = unitary2choi(str2unitary("1QX90p"))
    chi = choi2chi(choi)
    ideal = [0.5 0.0-0.5im 0.0 0.0;
             0.0-0.5im 0.5 0.0 0.0;
             0.0 0.0 0.0 0.0;
             0.0 0.0 0.0 0.0]
    @test isapprox(chi, ideal, atol=1e-14)
end

function test_processFidelity()
    gates = Dict("CNOT12" => Matrix{ComplexF64}([[1,0,0,0] [0,1,0,0] [0,0,0,1] [0,0,1,0]]),
                "InvCNOT12" => Matrix{ComplexF64}([[0,1,0,0] [1,0,0,0] [0,0,1,0] [0,0,0,1]]),
                "CNOT21" => Matrix{ComplexF64}([[1,0,0,0] [0,0,0,1] [0,0,1,0] [0,1,0,0]]),
                "1QId" => I,
                "1QX90p" => exp(-1im*pi*X/4),
                "1QX90m" => exp(1im*pi*X/4),
                "1QY90p" => exp(-1im*pi*Y/4),
                "1QY90m" => exp(1im*pi*Y/4),
                "1QXp" => exp(-1im*pi*X/2),
                "1QYp" => exp(-1im*pi*Y/2),
                "1QX45p" => exp(-1im*pi*X/8),
                "1QX22p" => exp(-1im*pi*X/16),
                "1QHad" => exp(-1im*(pi/2)*(1/sqrt(2))*(X+Z)),
                "1QZ90" => exp(-1im*(pi/4)*Z),
                "1QT" => exp(-1im*(pi/8)*Z),
                "II" => kron(I,I),
                "XI" => exp(-1im*kron(X,I)*pi/2),
                "IX" => exp(-1im*kron(I,X)*pi/2),
                "YI" => exp(-1im*kron(Y,I)*pi/2),
                "IY" => exp(-1im*kron(I,Y)*pi/2),
                "X_2I" => exp(-1im*kron(X,I)*pi/4),
                "IX_2" => exp(-1im*kron(I,X)*pi/4),
                "Y_2I" => exp(-1im*kron(Y,I)*pi/4),
                "IY_2" => exp(-1im*kron(I,Y)*pi/4),
                "XZ90" => exp(-1im*kron(X,Z)*pi/4),
                "ZX90" => exp(-1im*kron(Z,X)*pi/4),
                "XZm90" => exp(1im*kron(X,Z)*pi/4),
                "ZXm90" => exp(1im*kron(Z,X)*pi/4),
                "ZXm90_2" => exp(1im*kron(Z,I)*pi/2)*exp(1im*kron(I,X)*pi/4),
                "X_8I" => exp(-1im*kron(X,I)*pi/16),
                "X_4I" => exp(-1im*kron(X,I)*pi/8),
                "IX_8" => exp(-1im*kron(I,X)*pi/16),
                "IX_4" => exp(-1im*kron(I,X)*pi/8),
                "CNOT12_y" => Matrix{ComplexF64}([[1,0,0,0] [0,1,0,0] [0,0,-1im,0] [0,0,0,1im]]),
                "X90ZX90" => Matrix{ComplexF64}([[1,0,0,0] [0,0,0,-1im] [0,0,1,0] [0,-1im,0,0]]),
                "XX" => exp(-1im*kron(X,I)*pi/2)*exp(-1im*kron(I,X)*pi/2),
                "X90X90" => exp(-1im*kron(X,I)*pi/4)*exp(-1im*kron(I,X)*pi/4),
                "Y90Y90" => exp(-1im*kron(Y,I)*pi/4)*exp(-1im*kron(I,Y)*pi/4)
            );
    for gate in keys(gates)
        proces_fid, gate_fid = getProcessFidelity(str2unitary(gate), gate)
        @test isapprox(proces_fid, 0.0, atol=1e-14)
end

do_unitary2choi()
do_choi2pauliMap()
do_choi2chi()
test_processFidelity()

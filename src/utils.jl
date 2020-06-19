using Cliffords
import LinearAlgebra

"""
    unitary2choi(U::AbstractArray)

Returns the Choi super operator representation of a unitary.
"""
function unitary2choi(U::AbstractArray)

    d=size(U,1);
    choi=zeros(ComplexF64,d^2,d^2);
    for ii in 1:d
        for jj in 1:d
            proj = zeros(ComplexF64,d,d);
            proj[ii,jj] = 1;
            choi = choi + kron(proj,U*proj*U')/d;
        end
    end
    return choi
end

"""
    choi2pauliMap(choi::AbstractArray)

Converts a Choi representation to a Pauli Map representation.
"""
function choi2pauliMap(choi::AbstractArray)

    #Dimension of superoperator
    d2 = size(choi,1);

    #Create the Pauli opearators for n qubits
    pauliOps = allpaulis(sqrt(d2));

    pauliMap = zeros(d2,d2);
    for ct1 in 1:d2
        for ct2 in 1:d2
            pauliMap[ct2,ct1] = real(LinearAlgebra.tr(choi*pauliOps[ct1, ct2]));
        end
    end
    return pauliMap
end

"""
    unitary2pauli(U:AbstractArray)
"""
function unitary2pauli(U::AbstractArray)
    return choi2pauliMap(unitary2choi(U));
end

"""
    str2unitary

Converts a string description of a gate into a unitary matrix.
"""
function str2unitary(strIn::String)

    #Basic Pauli operators
    X = complex(Pauli(1));
    Y = complex(Pauli(3));
    Z = complex(Pauli(2));
    I = complex(Pauli(4));

    Uout = Dict("CNOT12" => Matrix{ComplexF64}([[1,0,0,0] [0,1,0,0] [0,0,0,1] [0,0,1,0]]),
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
                "Id" => kron(I,I),
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
    return Uout[strIn]
end

"""
    choi2chi(choi::AbstractArray)

Converts a Choi superoperator to a chi represtation in Pauli basis.
"""
function choi2chi(choi::AbstractArray)

    # Some dimensions
    d2 = size(choi,1);
    d = Int(sqrt(d2));

    # Get the Kraus operators from the eigen decomposition of the Choi
    vals, vecs = LinearAlgebra.eigen(choi);
    vals = LinearAlgebra.diagm(vals);

    chi = zeros(ComplexF64,d2,d2);

    pauliOps = allpaulis(log2(d));

    # Transform from the Krauss basis to the Pauli basis
    for kraussct in 1:size(vals,1)
        tmpKrauss = reshape(vecs[:,kraussct],d,d)*sqrt(d); # Krauss operator should have norm of d
        for paulict1 in 1:d2
            pauliLeft = LinearAlgebra.tr(complex(pauliOps[paulict1])*tmpKrauss)/d;
            for paulict2 in 1:d2
                pauliRight = LinearAlgebra.tr(complex(pauliOps[paulict2])*tmpKrauss')/d;
                chi[paulict1, paulict2] = chi[paulict1, paulict2] + vals[kraussct]*pauliLeft*pauliRight;
            end
        end
    end
    return chi
end

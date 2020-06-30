using Cliffords
import LinearAlgebra

###############################################################################
## Useful tools for manipulation of matrix representations in a tomography ####
###############################################################################

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
    num_qubits = Int(log2(sqrt(d2)))
    pauliOps = map(x -> complex(x), vec(permutedims(allpaulis(num_qubits),
                                    reverse(collect(1:num_qubits)))));

    pauliMap = zeros(d2,d2);
    for ct1 in 1:d2
        p1 = pauliOps[ct1]
        for ct2 in 1:d2
            p2 = pauliOps[ct2]
            pauliMap[ct2,ct1] = real(LinearAlgebra.tr(choi*kron(transpose(p2), p1)));
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

    # note the ordering diff between Matlab
    vals = LinearAlgebra.diagm(reverse(vals));
    vecs = reverse(vecs, dims=2);

    chi = zeros(ComplexF64,d2,d2);

    pauliOps = map(x -> complex(x), allpaulis(log2(d)));

    # Transform from the Krauss basis to the Pauli basis
    for kraussct in 1:size(vals,1)
        tmpKrauss = reshape(vecs[:,kraussct],d,d)*sqrt(d); # Krauss operator should have norm of d
        for paulict1 in 1:d2
            pauliLeft = LinearAlgebra.tr(pauliOps[paulict1]*tmpKrauss)/d;
            for paulict2 in 1:d2
                pauliRight = LinearAlgebra.tr(pauliOps[paulict2]*transpose(tmpKrauss))/d;
                chi[paulict1, paulict2] = chi[paulict1, paulict2] + vals[kraussct]*pauliLeft*pauliRight;
            end
        end
    end
    return chi
end

###############################################################################
## Fidelity tools for tomography ##############################################
###############################################################################

function getProcessFidelity(choi::AbstractArray, idealProcess::String)
    nbrQubits = log2(sqrt(size(choi, 1)))
    # Calculate the overlaps with the ideal gate
    if typeof(idealProcess) == String
        unitaryIdeal = str2unitary(idealProcess);
    else
        # assume we are passed a matrix otherwise
        unitaryIdeal = idealProcess;
    end
    choiIdeal = unitary2choi(unitaryIdeal);

    # Convert to chi representation to compute fidelity metrics
    chiExp = choi2chi(choi);
    chiIdeal = choi2chi(choiIdeal);

    processFidelity = real(LinearAlgebra.tr(chiExp*chiIdeal));
    gateFidelity = (2^nbrQubits*processFidelity+1)/(2^nbrQubits+1);

    return processFidelity, gateFidelity
end

# TO-DO: implement plotting

# Create the pauli operator strings
# pauliStrs = allpaulis(nbrQubits);

# Create the pauli map for plotting
# pauliMapIdeal = choi2pauliMap(choiIdeal);
# pauliMapLSQ = choi2pauliMap(choiLSQ);
# pauliMapExp = choi2pauliMap(choiSDP);

"""
    PauliPlot(U::AbstractMatrix; save_fig=false)

Plot the Pauli transfer matrix of an operator U.  Only one and two qubit
operators are supported.

Ex:
  julia> PauliPlot(str2unitary("1QY90p"), save_fig=true)
"""
function PauliPlot(U::AbstractMatrix; save_fig::Bool=false)
    num_qubits = Int(log2(size(U,1)))
    op_num = 0:(4^num_qubits - 1)
    # This code snippet works for n-qubit unitaries, it's just not that
    # verbose in what its doing...  It creates the Pauli string labels.
    snip(s::String) = s[nextind(s,1):end];
    labels = map(p -> snip(string(p)),
                 vec(permutedims(allpaulis(num_qubits),
                 reverse(collect(1:num_qubits)))))
    # if num_qubits == 1
    #     labels = ["I", "X", "Y", "Z"]
    # elseif num_qubits == 2
    #     labels = ["II","IX","IY","IZ",
    #               "XI","XX","XY","XZ",
    #               "YI","YX","YY","YZ",
    #               "ZI","ZX","ZY","ZZ"]
    # else
    #     labels = []
    # end
    fig = figure()
    ax = fig.add_subplot(111)
    cax = ax.matshow(unitary2pauli(U), cmap="bwr")
    fig.colorbar(cax)

    ax.set_xticklabels(labels)
    ax.set_yticklabels(labels)
    ax.set_xticks(op_num)
    ax.set_yticks(op_num)
    ax.xaxis.set_ticks_position("bottom")
    ax.grid(color="grey", linestyle="--", alpha=0.1)
    ax.spines["right"].set_visible(false)
    ax.spines["top"].set_visible(false)
    ax.spines["left"].set_visible(false)
    ax.spines["bottom"].set_visible(false)
    ax.set_xlabel("Input Pauli Operator")
    ax.set_ylabel("Output Pauli Operator")

    if save_fig
        fig.savefig("PauliPlot.png")
        println("Please rename the saved figure!")
    end
    fig.show()
end

using Cliffords
import LinearAlgebra

# helper to remove the first "+/-" char from Pauli string labels
snip(s::String) = s[nextind(s,1):end];

###############################################################################
## Useful tools for manipulation of matrix representations in a tomography ####
###############################################################################

"""
    unitary2choi(U::AbstractArray)

Returns the Choi super operator representation of a unitary.

# Arguments
- `U`: the unitary matrix representation of the process

# Returns
- `choi`: a choi matrix representation

# Examples
julia> choi = unitary2choi(str2unitary("1QX90p"))
4×4 Array{Complex{Float64},2}:
 0.25+0.0im    0.0+0.25im   0.0+0.25im  0.25+0.0im
  0.0-0.25im  0.25+0.0im   0.25+0.0im    0.0-0.25im
  0.0-0.25im  0.25+0.0im   0.25+0.0im    0.0-0.25im
 0.25+0.0im    0.0+0.25im   0.0+0.25im  0.25+0.0im
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

# Arguments
- `choi`: the matrix representation of the process

# Returns
- `pauli_map`: a Pauli representation of the choi matrix

# Examples
julia> choi = unitary2choi(str2unitary("1QX90p"))
julia> choi2pauliMap(choi)
4×4 Array{Float64,2}:
  1.0          0.0   0.0          0.0
  0.0          1.0   0.0          0.0
  0.0          0.0   3.33067e-16  1.0
 -5.55112e-17  0.0  -1.0          3.33067e-16
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

Compute the Pauli representation of an arbitrary unitary

# Arguments
- `U`: the matrix representation of the unitary.

# Returns
- pauli_map: Matrix{ComplexF64}

# Examples
julia> unitary2pauli(str2unitary("1QY90p"))
4×4 Array{Float64,2}:
  1.0          0.0          0.0   0.0
  0.0          3.33067e-16  0.0  -1.0
  0.0          0.0          1.0   0.0
 -5.55112e-17  1.0          0.0   3.33067e-16
"""
function unitary2pauli(U::AbstractArray)
    return choi2pauliMap(unitary2choi(U));
end

"""
    str2unitary(strIn::String)

Converts a string description of a gate into a unitary matrix.  Supported
unitaries are: X_4I, IX_4, IX_2, InvCNOT12, CNOT12, CNOT21, 1QXp, 1QX22p,
IY_2, IY, ZXm90_2, 1QHad, XX, 1QId, IX_8, 1QT, YI, X90ZX90, CNOT12_y, X90X90,
1QX45p, ZXm90, 1QX90p, ZX90, 1QY90p, Y90Y90, 1QYp, X_8I, 1QY90m, X_2I, XZm90,
XI, XZ90, 1QX90m, 1QZ90, II, IX, Y_2I

# Arguments
- `strIn`: string id of the unitary

# Returns
- `U`: Matrix{ComplexF64} representing the unitary transformation

# Examples
julia> str2unitary("1QX90p")
2×2 Array{Complex{Float64},2}:
 0.707107+0.0im            0.0-0.707107im
      0.0-0.707107im  0.707107+0.0im
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

# Arguments
- `choi`: the matrix representation of the process

# Returns
- `chi`: chi matrix representation of the process

# Examples
julia> using RandomQuantum
julia> U = rand(RandomClosedEvolution(2,2))
julia> choi = unitary2choi(U)
julia> choi2chi(choi)
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

    #Create the Pauli opearators for n qubits
    num_qubits = Int(log2(sqrt(d2)));
    pauliOps = map(x -> complex(x), vec(permutedims(allpaulis(num_qubits),
                                    reverse(collect(1:num_qubits)))));

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
## Fidelity and plotting tools for tomography #################################
###############################################################################
"""
    getProcessFidelity(choi::AbstractArray, idealProcess::String)

Compute the gate and process fidelity for a give Choi matrix and a given ideal
process matrix.

# Arguments
- `choi`: the matrix representation of the process
- `idealProcess`: string representing the ideal process in str2unitary()

# Returns
- `processFidelity`: -> real(tr(choi*choiIdeal))
- `gateFidelity`: -> (d*processFidelity +1)/(d + 1)

# Examples
julia> getProcessFidelity(choi, "XI")
"""
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

    processFidelity = real(LinearAlgebra.tr(choi*choiIdeal));
    gateFidelity = (2^nbrQubits*processFidelity+1)/(2^nbrQubits+1);

    return processFidelity, gateFidelity
end

"""
    getProcessFidelity(choi::AbstractArray, idealProcess::AbstractMatrix)

Compute the gate and process fidelity for a give Choi matrix and a given ideal
unitary.

# Arguments
- `choi`: the matrix representation of the process
- `idealProcess`: matrix representing the ideal process unitary

# Returns
- `processFidelity`: -> real(tr(choi*choiIdeal))
- `gateFidelity`: -> (d*processFidelity +1)/(d + 1)

# Examples
julia> getProcessFidelity(choi, U)
"""
function getProcessFidelity(choi::AbstractArray, idealProcess::AbstractMatrix)
    nbrQubits = log2(sqrt(size(choi, 1)))
    # Calculate the overlaps with the ideal gate
    choiIdeal = unitary2choi(idealProcess);

    processFidelity = real(LinearAlgebra.tr(choi*choiIdeal));
    gateFidelity = (2^nbrQubits*processFidelity+1)/(2^nbrQubits+1);

    return processFidelity, gateFidelity
end

"""
    PlotProcess(U::AbstractMatrix; p_map::String="pauli", save_fig::Bool=false)

Plot the Pauli transfer matrix of an operator U.  Only one and two qubit
operators are supported.

# Arguments
- `U`: the matrix representation of the unitary you want to plot.
- `p_map`: which process map to use {"pauli", "chi"}.  "pauli" will plot the
           Pauli transfer matrix and "chi" will plot the chi matrix.
- `save_fig`: set to true if you want a .png of the figure saved in the current
              directory.

# Returns
- void: The function outputs a plot

# Examples
julia> PauliPlot(str2unitary("1QY90p"), save_fig=true)
"""
function PlotProcess(U::AbstractMatrix;
                     p_map::String="pauli",
                     save_fig::Bool=false)

    num_qubits = Int(log2(size(U,1)))
    op_num = 0:(4^num_qubits - 1)
    # This code snippet works for n-qubit unitaries, it's just not that
    # verbose...
    labels = map(p -> snip(string(p)),
                 vec(permutedims(allpaulis(num_qubits),
                 reverse(collect(1:num_qubits)))))
    # a more verbose and brittle option would be something like:
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
    if p_map == "pauli"
        cax = ax.matshow(unitary2pauli(U), cmap="bwr")
        ax.set_title("Pauli transfer matrix")
    elseif p_map == "chi"
        cax = ax.matshow(real(choi2chi(unitary2choi(U))), cmap="bwr")
        ax.set_title(L"$\chi$ process matrix")
    end

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
        fig.savefig("ProcessPlot.png")
        println("Please rename the saved figure!")
    end
    fig.show()
end

"""
    PauliToCounts(Z::AbstractArray,nbrQubits)

Convert a set of 2^nbrQubits-1 of SigmaZ_j and correlator SigmaZ_iSigmaZ_j... expectation values into state populations (expectation value of projection operators)

# Arguments
- `Z`: array of SigmaZ expectation values. Z[k] is the expectation value of a product of I and Z operators corresponding to the binary value of k.
        So for example if k=0110 , Z[k]=<IZZI> . Z needs to have size of 2^nbrQubits-1
- `nbrQubits`: number of Qubits.

# Returns
- Array of projector expectation values

# Examples
julia> PauliToCounts([1,1,1],2)
"""
function PauliToCounts(Z::AbstractArray,nbrQubits)
  P = ones(1,2^nbrQubits)

  for k in 0:2^nbrQubits-1
       kbin = digits!(zeros(Int64,nbrQubits),k,base=2)
       kbin = -2*kbin .+ 1

       for h in 1:2^nbrQubits-1
           hbin = digits!(zeros(Int64,nbrQubits),h,base=2)
           P[k+1] += prod(sign.(kbin.*hbin .+ 0.5))*Z[h]
       end

   end
   return P/2^nbrQubits
end

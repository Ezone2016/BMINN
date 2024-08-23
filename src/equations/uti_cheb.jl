function chebdif(N::Int64, M::Int64)
    eye = LinearAlgebra.I(N);
    n1 = Int(floor(N/2)); n2  = Int(ceil(N/2));     # Indices used for flipping trick.

    k = collect(0:N-1);                        # Compute theta vector.
    th = k.*pi./(N-1)

    x = sin.(pi.*collect(N-1:-2:1-N)./(2*(N-1))); # Compute Chebyshev points.

    T = repeat(th/2,1,N);
    DX = 2*sin.(T'+T).*sin.(T'-T);          # Trigonometric identity.
    DX = [DX[1:n1,:]; -reverse(reverse(DX[1:n2,:], dims = 2), dims=1)];   # Flipping trick.
    DX[eye.==true] = ones(N,1);                      # Put 1's on the main diagonal of DX.

    C = convert(Array{Float64},Matrix(Toeplitz((-1).^vec(k),(-1).^vec(k))));            # C is the matrix with
    C[1,:] = C[1,:].*2; C[N,:] = C[N,:].*2;     # entries c[k]/c[j]
    C[:,1] = C[:,1]/2; C[:,N] = C[:,N]/2

    Z = 1 ./DX;                              # Z contains entries 1/(x[k]-x[j])
    Z[eye.==true] = zeros(N,1);              # with zeros on the diagonal.

    D  = eye;                                 # D contains diff(). matrices.
    DM = zeros(N,N,M);

    for ell = 1:M
        D = ell*Z.*(C.*repeat(diag(D),1,N) - D); # Off-diagonals
        D[eye.==true] = -sum(D', dims = 1);                            # Correct main diagonal of D
        DM[:,:,ell] = D;                                   # Store current D in DM
    end

    return x,DM
end

function cumsummat(N::Int64)
    N = N-1;
    # Matrix mapping coeffs -> values.
    T = cp2cdm(N);
    # Matrix mapping values -> coeffs.
    Tinv = cd2cpm(N);

    # Matrix mapping coeffs -> integral coeffs. Note that the highest order
    # term is truncated.
    k = 1:N;
    k2 = collect(2*(k .-1));  k2[1] = 1;  # avoid divide by zero
    B = zeros(N+1,N+1);
    B[diagind(B, -1)] = 1 ./(2*k);
    B[diagind(B, 1)] .= -1 ./k2
    v = ones(N); v[2:2:end] .= -1;
    B[1,:] = sum(Diagonal(v)*B[2:N+1,:], dims=1);
    B[:,1] = 2*B[:,1];

    Q = T*B*Tinv;
end

function cp2cdm(N::Int64)
    # Values of Cheb. polys at Cheb nodes, x(n)=-cos(pi*n/N).
    T = cos.( pi*(N:-1:0)/N*(0:N)' );
end

function cd2cpm(N::Int64)
    # Three steps: Double the data around the circle, apply the DFT matrix,
    # and then take half the result with 0.5 factor at the ends.
    F = exp.( -1im*((pi/N)*(0:2*N-1))*(0:2*N-1)' );  # DFT matrix
    rows = 1:N+1;  # output upper half only
    # Impose symmetries on data and coeffs.
    C = real( [ F[rows,N+1] F[rows,N:-1:2]+F[rows,N+2:2*N] F[rows,1] ] );
    C = C/N;  C[[1, N+1],:] = 0.5*C[[1, N+1],:];
    return C
end

function clencurt(N)
  thet = pi*(0:N)/N; x = cos.(thet);
  w = zeros(N+1); ii = 2:N; v = ones(N-1);

  if mod(N,2)==0
    w[1] = 1/(N^2-1); w[N+1] = w[1];
    for k in 1:N/2-1
        v -= 2*cos.(2*k*thet[ii])/(4*k^2-1);
    end
    v = v - cos.(N*thet[ii])/(N^2-1);
  else
    w[1] = 1/N^2; w[N+1] = w[1];
    for k in 1:(N-1)/2
        v -= 2*cos.(2*k*thet[ii])/(4*k^2-1);
    end
  end

  w[ii] = 2*v/N;
  return x,w
end



function sigma(k, l, number)
    s = spzeros(number, number)  
    s[k, l] = 1
    return s
end



function zero_except_one!(matrix, row, col)
    # Check if the specified row and column are within the bounds of the matrix
    if row > size(matrix, 1) || col > size(matrix, 2) || row < 1 || col < 1
        error("Row or column is out of bounds.")
    end

    # Get the value of the element to keep
    value_to_keep = matrix[row, col]

    # Set all elements to zero
    fill!(matrix, 0)

    # Set the specified element to its original value
    matrix[row, col] = value_to_keep

    return matrix
end


function element_non_zero(line, column, value )
    global  dim
    A=zeros(ComplexF64, dim, dim)
    for j=1:Int(length(line))
        A[line[j], column[j]]=value[j]
    end
    return A
end

function keep_column(matrix, column_number)
 
    new_matrix = zeros(size(matrix))
    # Copy the specified column from the original matrix to the new matrix
    new_matrix[:, column_number] = matrix[:, column_number]
    return new_matrix
end


function column_non_zero(matrix, column )
    A=zeros(ComplexF64, dim, dim)
    A[:, column]=vector
    return A
end

function greens(r1, r2)
    if r1 == r2
        return 1im / 2
    else
        # Calculate r and its magnitude
        r = r1 - r2
        r_mag = norm(r) # Magnitude of r vector
        # Calculate e^(ik0*r)
        exponential_term = exp(1im * r_mag)

        # For linear polarization along the z-direction, driving in x diraction
        p1 = [0, 0, 1]

        # For circular polarization
        #p1 = 1/(sqrt(2))*[1, -1im, 0]

        # Calculate outer product of r
        outer_product = dot(p1 , r)^2 / r_mag^2

        # Calculate Green's function
        greens_value = 3 * (exponential_term / (4 * r_mag^3)) *
                       (( r_mag^2 + 1im  * r_mag - 1) +
                        (-r_mag^2 - 3im  * r_mag + 3) * outer_product)
        # greens_value = 3 * (exponential_term / (4 * r_mag^3)) *
        #                 (( r_mag^2 ) + (-r_mag^2 ) * outer_product)  #only far field
        return greens_value
    end
end


function subspaceG(number)
    global N
    A=zeros(ComplexF64, 2*N-1,2*N-1)
    if number==1
        A[1:N, 1:N].=1
    end
    if number==2
        A[N+1:2*N-1, N+1:2*N-1].=1
    end
    if number==3
        A[1:N, N+1:2*N-1].=1
    end
    if number==4
        A[N+1:2*N-1, 1:N].=1
    end
    return A
end

function subspace(number)
    global N
    A=zeros(ComplexF64, 2*N-1, 2*N-1)
    if number==1
        A[1:N, 1:N].=1
    end
    if number==2
        A[N+1:2*N-1, N+1:2*N-1].=1
    end
    return A
end

function create_g_combined(g)
    global dim,N

    g_global_middle=[ g g[:,1:N-1]]
    g_global_low=[g[1:N-1,:] g[1:N-1,1:N-1] ]
    g_global =[ g_global_middle; g_global_low]
    return g_global
end
    

function create_block_diagonal(matrices::AbstractMatrix...)
    total_rows = sum(size(m, 1) for m in matrices)
    total_cols = sum(size(m, 2) for m in matrices)

    block_diag = spzeros(ComplexF64, total_rows, total_cols)
    row_start = 1
    col_start = 1

    for m in matrices
        rows, cols = size(m)
        block_diag[row_start:(row_start + rows - 1), col_start:(col_start + cols - 1)] = sparse(m)
        row_start += rows
        col_start += cols
    end

    return block_diag
end



#function evolution_rho(N,rho_matrix, gamma,  x)
function evolution_rho(rho_matrix, N, x, G)
    #the parameters of the problem
    dim=2*N+1
    omega=0.5/sqrt(N)
    Gamma = 2 * imag(G)
    H2 = spzeros(ComplexF64, dim, dim)
    H3 = -sparse(create_block_diagonal(zeros(2,2), G, G[1:N-1, 1:N-1]))
    for i = 3:N+2
        H2[i, 1] = omega / 2 * exp(1im * x[i-2])
    end
    for i = N+3:2*N+1
        H2[i, 2] = omega / 2 * exp(1im * x[i-(N+2)])
    end
    H = copy(H3)  
    H.+=H2+H2'

    
    # Define the Hamiltonian H
    #excited= sum(sigma(2+i, 2+i, dim) for i in 1:2N-1)
    #H1 = -Delta0 * excited 
    #vector1 = [zeros(Complex{Float64}, 2); exp.(1im .* x); zeros(Complex{Float64}, N - 1)]
    #vector2 = [zeros(Complex{Float64}, N + 2); exp.(1im .* x[1:N-1])]
    #H2[:, 1] = omega / 2 * vector1
    #H2[:, 2] = omega / 2 * vector2
    #H2=omega/2 *(element_non_zero(range1, ones(Int,N), [exp(1im *x[i]) for i=1:N] )+ element_non_zero(range2, 2*(ones(Int,N-1)), [exp(1im *x[i]) for i=1:N-1] ))
    #H2=omega/2 *sum(exp(1im *x[i])*(sigma(2+i, 1, dim)+sigma(2+N+i, 2, dim)) for i=1:N-1) +omega/2 *exp(1im * x[N]) * sigma(2 + N, 1, dim)
    #H3=-create_block_diagonal(zeros(2,2), G, G[1:N-1, 1:N-1])

    #H=H1+H2+H3+H2'
    #H3= sum(-(gamma[i,j])* sigma(2+i,1,dim)*sigma(1,2+j,dim) for i=1:N, j=1:N)+sum(-(gamma[i,j] )*sigma(2+N+i,2,dim)*sigma(2,2+N+j, dim) for i=1:N-1, j=1:N-1)
    #H3= sum(Delta[i,j]*sigma(2+i,1,dim)*sigma(1,2+j,dim) for i=1:N, j=1:N)+sum(Delta[i,j]*sigma(2+N+i,2,dim)*sigma(2,2+N+j, dim) for i=1:N-1, j=1:N-1)

    commutator = -1im * (H * rho_matrix - rho_matrix * H')
    
    dissipator = spzeros(Complex{Float64}, dim, dim)
    #the repopulation terms. Can i write it as a matrix multiplication?
    #dissipator[1,1] =sum( Gamma[i, j]*  rho_matrix[2+j,2+i] for i=1:N, j=1:N )
    #dissipator =sum( Gamma[i, j]*  sigma(1, 2+j, dim) * rho_matrix * sigma( 2+i,1, dim) for i=1:N, j=1:N )+sum( Gamma[i, j]*  sigma(2, 2+N+j, dim) * rho_matrix * sigma( 2+N+i,2, dim) for i=1:N-1, j=1:N-1 ) 
    #dissipator +=sum( Gamma[i, j]*  sigma(2, 2+j, dim) * rho_matrix * sigma( 2+N+i,1, dim) for j=1:N, i=1:N-1 )+sum( Gamma[i, j]*  sigma(1, 2+N+j, dim) * rho_matrix * sigma( 2+i,2, dim) for j=1:N-1, i=1:N )
  
    dissipator[1,1]=sum( Gamma[i, j]*  rho_matrix[2+j,2+i] for i=1:N, j=1:N )
    dissipator[2,2]=sum( Gamma[i, j]*  rho_matrix[2+N+j,2+N+i] for i=1:N-1, j=1:N-1 ) 
    dissipator[2,1]=sum( Gamma[i, j]* rho_matrix[2+j,2+N+i] for j=1:N, i=1:N-1 )
    dissipator[1,2]=conj(dissipator[2,1])
  
    #ensure hermitian result
    #commutator .= 0.5 .* (commutator + commutator') 
    #dissipator .= 0.5 .* (dissipator + dissipator')
    return  commutator + dissipator
end 



function create_total_G(G)
    # Assuming G is a predefined N x N matrix
    N = size(G, 1)

    # Top-right block: G with its last column removed
    top_right = G[:, 1:N]

    # Bottom-left block: G with its last row removed
    bottom_left =G[1:N, :]

    # Bottom-right block: G with both its last row and column removed
    bottom_right =  G[1:N, 1:N]

    # Concatenate the blocks
    top_row = hcat(G, top_right)
    bottom_row = hcat(bottom_left, bottom_right)

    # Final matrix
    new_matrix = vcat(top_row, bottom_row)
    return new_matrix
end


function diag_G(N, x, G, omeg, gamma0)
    imG=imag(sparse(create_block_diagonal(zeros(2,2),gamma0*G,gamma0* G[1:N-1, 1:N-1]))   )
    eig_result=eigen(imG)
    V = eig_result.vectors
    D = Diagonal(eig_result.values)
    return V,D
end


function h_nh(N, F, omeg)
    dim=2*N+1
    omega=omeg/sqrt(N)
    H2 = spzeros(ComplexF64, dim, dim)
    for i = 3:N+2
        H2[i, 1] = omega / 2 * exp(1im * x[i-2])
    end
    for i = N+3:2*N+1
        H2[i, 2] = omega / 2 * exp(1im * x[i-(N+2)])
    end
    H = copy(H3)  
    H.+=H2+H2'

end

function evolution_rho_effective_good(rho_matrix, N, Heff, Leff, Leff_dagger)
    commutator = -1im * (Heff * rho_matrix - rho_matrix * Heff)
    dissipator = zeros(ComplexF64,2,2) 
    for k = 1:N
        Leffk_dagger= Leff_dagger[k,:,:]
        Leffk= Leff[k,:,:]
        dissipator +=Leffk*rho_matrix*Leffk_dagger -1/2* Leffk_dagger*Leffk*rho_matrix-1/2* rho_matrix* Leffk_dagger*Leffk
    end
    #println(commutator[1,1], dissipator[1,1])
    res=commutator .+ dissipator
    return  res
end 

function evolution_rho_effective(rho_matrix, N, Heff, Leff, Leff_dagger)
    #commutator = -1im * (Heff * rho_matrix - rho_matrix * Heff)
    dissipator = zeros(ComplexF64,2,2) 
    for k = 1:N
        Leffk_dagger= Leff_dagger[k,:,:]
        Leffk= Leff[k,:,:]
        dissipator +=Leffk*rho_matrix*Leffk_dagger -1/2* Leffk_dagger*Leffk*rho_matrix-1/2* rho_matrix* Leffk_dagger*Leffk
    end
    #println(commutator[1,1], dissipator[1,1])
    #res=commutator .+ dissipator
    return  dissipator
end 


function evolution_rho_omega(rho_matrix, N, x, G, omeg, gamma0)
    #the parameters of the problem
    dim=2*N+1
    omega=omeg/sqrt(N)
    Gamma = 2 * imag(G)*gamma0
    H2 = spzeros(ComplexF64, dim, dim)
    
    #H3 = -sparse(create_block_diagonal(zeros(2,2),gamma0*create_total_G(G)))
    H3 = -sparse(create_block_diagonal(zeros(2,2),gamma0*G,gamma0* G[1:N-1, 1:N-1]))
    for i = 3:N+2
        H2[i, 1] = omega / 2 * exp(1im * x[i-2])
        #H2[i, 2] = omega / 2 * exp(1im * x[i-2])
    end
    for i = N+3:2*N+1
        H2[i, 2] = omega / 2 * exp(1im * x[i-(N+2)])
        #H2[i, 1] = omega / 2 * exp(1im * x[i-(N+2)])
    end
    H = copy(H3)  
    H.+=H2+H2'


    commutator = -1im * (H * rho_matrix - rho_matrix * H')
    
    dissipator = spzeros(Complex{Float64}, dim, dim)
    #the repopulation terms. Can i write it as a matrix multiplication?
    #dissipator[1,1] =sum( Gamma[i, j]*  rho_matrix[2+j,2+i] for i=1:N, j=1:N )
    #dissipator =sum( Gamma[i, j]*  sigma(1, 2+j, dim) * rho_matrix * sigma( 2+i,1, dim) for i=1:N, j=1:N )+sum( Gamma[i, j]*  sigma(2, 2+N+j, dim) * rho_matrix * sigma( 2+N+i,2, dim) for i=1:N-1, j=1:N-1 ) 
    #dissipator +=sum( Gamma[i, j]*  sigma(2, 2+j, dim) * rho_matrix * sigma( 2+N+i,1, dim) for j=1:N, i=1:N-1 )+sum( Gamma[i, j]*  sigma(1, 2+N+j, dim) * rho_matrix * sigma( 2+i,2, dim) for j=1:N-1, i=1:N )
  
    dissipator[1,1]=sum( Gamma[i, j]*  rho_matrix[2+j,2+i] for i=1:N, j=1:N )
    dissipator[2,2]=sum( Gamma[i, j]*  rho_matrix[2+N+j,2+N+i] for i=1:N-1, j=1:N-1 ) 
    dissipator[2,1]=sum( Gamma[i, j]* rho_matrix[2+j,2+N+i] for j=1:N, i=1:N-1 )
    dissipator[1,2]=conj(dissipator[2,1])
  
    #ensure hermitian result
    #commutator .= 0.5 .* (commutator + commutator') 
    #dissipator .= 0.5 .* (dissipator + dissipator')
    return  commutator + dissipator
end 

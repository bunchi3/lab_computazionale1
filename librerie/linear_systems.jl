#File's index
#1. Error handling section: necessario per messaggi di errore
#2. fw_sub(matrix, b) function: risolve un sistema lineare con matrice triangolare inferiore
#3. bw_sub(matrix, b) function: risolve un sistema lineare con matrice triangolare superiore
#4. outer_product(v, w) function: esegue il prodotto esterno tra due vettori
#5. LU_dec(A) function: esegue la decomposizione di una matrice in trinagolare inferiore e superiore
#6. chlsky_dec(A) function: esegue la decomposizione di una matrice definita positiva in trinagolare inferiore e superiore
#7. diag_prod(squared_matr) function: calcola il prodotto degli elementi diagonali di una matrice n*n qualunque
#8. tri_det(tri_matr) function: calcola il determinante di una matrice n*n triangolare tramite diag_prod()
#9. my_det(A) function: calcola il determinante di una matrice tramite LU_dec() e tri_det()
#10. solve_linear_system(A, b) function: risolve un sistema lineare Ax = b tramite LU_dec(), fw_sub() e bw_sub()
#11. sls_chlsky(A, b) function: risolve un sistema lineare Ax = b tramite chlsky_dec(), fw_sub() e bw_sub()
#12. least_sq(A, b) function: calcola i coefficienti di un fit con metodo dei minimi quadrati
#-----------------------------------------------------------------------------------------
#Import section
using LinearAlgebra
#-----------------------------------------------------------------------------------------
#Error handling section
# Definisci un wrapper per la matrice
struct MatriceSummary{T}
    A::T
end
# Sovrascrivi il metodo show per MatriceSummary in modo che stampi solo le dimensioni
Base.show(io::IO, s::MatriceSummary) = print(io, "matrix of dimension ", size(s.A))
#-----------------------------------------------------------------------------------------
#This function operates a forward substitution for lower diagonal matrices (vector of vectors and matrices)
function fw_sub(matrix, b)
    if !(matrix isa AbstractMatrix) 
        throw(ArgumentError("The first parameter must be an AbstractMatrix."))
    end
    if !(b isa AbstractVector)
        if b isa AbstractMatrix
            b = vec(b)
        else
            throw(ArgumentError("The second parameter must be an AbstractVector or an AbstractMatrix."))
        end    
    end
    
    #Mi accerto che la matrice e il vettore siano di tipo Float64
    matrix = convert(Matrix{Float64}, matrix)
    b = convert(Vector{Float64}, b)

    if !istril(matrix) 
        throw(ArgumentError("The matrix must be lower triangular"))
    end 

    n = size(matrix, 1)
    if n != length(b)
        throw(DomainError(MatriceSummary(b), "Dimension mismatch between the matrix and the vector of known terms"))
    end

    #Creo un vettore di zeri con lo stesso tipo degli elementi di b
    solution = zeros(eltype(b), n)
    # i indice di riga, j indice di colonna
    for i in 1:n
        solution[i] = b[i] / matrix[i, i]
        for j in 1:i-1
            solution[i] -= matrix[i, j] * solution[j] / matrix[i, i]
        end
    end
    return solution
end
#------------------------------------------------------------------------------------------------------------------------
#This function operates a backward substitution for upper diagonal matrices (vector of vectors and matrices)
function bw_sub(matrix, b)
    if !(matrix isa AbstractMatrix) 
        throw(ArgumentError("The first parameter must be an AbstractMatrix."))
    end
    if !(b isa AbstractVector)
        if b isa AbstractMatrix
            b = vec(b)
        else
            throw(ArgumentError("The second parameter must be an AbstractVector or an AbstractMatrix."))
        end    
    end

    #Mi accerto che la matrice e il vettore siano di tipo Float64
    matrix = convert(Matrix{Float64}, matrix)
    b = convert(Vector{Float64}, b)
    if !istriu(matrix) 
        throw(ArgumentError("The matrix must be upper triangular"))
    end
    n = size(matrix, 1)
    if n != length(b)
        throw(DomainError(MatriceSummary(b), "Dimension mismatch between the matrix and the vector of known terms"))
    end

    #Creo un vettore di zeri con lo stesso tipo degli elementi di b
    solution = zeros(eltype(b), n)
    # i indice di riga, j indice di colonna
    for i in n:-1:1
        solution[i] = b[i] / matrix[i, i]
        for j in i+1:n
            solution[i] -= matrix[i, j] * solution[j] / matrix[i, i]
        end
    end
    return solution
end
#------------------------------------------------------------------------------------------------------------------------
#Funzione per il prodotto esterno. Restituisce una matrice a partire da due vettori
function outer_product(v, w)
    if !(v isa AbstractVector) 
        if v isa AbstractMatrix
            v = vec(v)
        else
            throw(ArgumentError("The first parameter must be an AbstractVector or an AbstractMatrix"))
        end
    end 
    if !(w isa AbstractVector)
        if w isa AbstractMatrix
            w = vec(w)
        else
            throw(ArgumentError("The second parameter must be an AbstractVector or an AbstractMatrix"))
        end    
    end
    #Mi accerto che i vettori siano di tipo Float64
    v = convert(Vector{Float64}, v)
    w = convert(Vector{Float64}, w)
    #Mi accerto che i vettori siano della stessa lunghezza
    if length(v) != length(w)
        throw(DomainError(MatriceSummary(v), "Dimension mismatch between the two arrays"))
    end
    n = length(v)
    matrix = zeros(n, n)
    for i in 1:n
        for j in 1:n
            matrix[i, j] = v[i]*w[j]
        end
    end
    return matrix
end
#---------------------------------------------------------------------------------------------------------------------
#Esegue la decomposizione di una matrice in due matrici, una triangolare inferiore e l'altra triangolare superiore
#NON INCLUDE ROW-PIVOTING
function LU_dec(A)
    if !(A isa AbstractMatrix)
        throw(ArgumentError("The parameter must be an AbstractMatrix"))
    end
    #Mi accerto che la matrice sia quadrata
    if size(A,1) != size(A,2)
        throw(DomainError(MatriceSummary(A), "The given matrix isn't squared."))
    end
    #Mi accerto che la matrice sia di tipo Float64
    A = convert(Matrix{Float64}, A)

    n = size(A, 1)
    L = diagm(ones(n))
    U = zeros(n, n)
    #B matrice da usare per i calcoli intermedi
    B = A
    #i indice dell'iterazione. Faccio in modo di escludere gli zeri
    for i in 1:n
        U[i, i:n] = B[i, i:n]
        L[i:n, i] = B[i:n, i]/U[i,i]
        B -= outer_product(L[:,i], U[i,:])       
    end
    return L, U
end
#-------------------------------------------------------------------------------------------------------------------------
#Esegue la decomposizione di una matrice definita positiva in due matrici, triangolare inf e sup (Cholesky)
#Restituisce solo la matrice triangolare superiore, l'altra è la trasposta 
function chlsky_dec(A)
    if !(A isa AbstractMatrix)
        throw(ArgumentError("The parameter must be an AbstractMatrix"))
    end
    #Mi accerto che la matrice sia quadrata
    if size(A,1) != size(A,2)
        throw(DomainError(MatriceSummary(A), "The given matrix isn't squared."))
    end
    #Mi accerto che la matrice sia di tipo Float64
    A = convert(Matrix{Float64}, A)

    n = size(A, 1)
    RT = zeros(n, n)    #analogo di RT
    R = zeros(n, n)     #analogo di U
    #B matrice da usare per i calcoli intermedi
    B = A
    #i indice dell'iterazione. Faccio in modo di escludere gli zeri
    for i in 1:n
        if B[i,i] <= 0
            throw(DomainError(MatriceSummary(A), "La matrice non è definita positiva."))
        else
            R[i, i:n] = B[i, i:n]/sqrt(B[i,i])
            RT[i:n, i] = R[i, i:n]
            B -= outer_product(RT[:,i], R[i,:])
        end      
    end
    return R
end
#---------------------------------------------------------------------------------------------------------------------
#Calcola il prodotto degli elementi diagonali di una matrice n*n 
function diag_prod(squared_matr)
    if !(squared_matr isa AbstractMatrix)
        throw(ArgumentError("The parameter must be an AbstractMatrix"))
    end
    #Mi accerto che la matrice sia quadrata
    if size(squared_matr, 1) != size(squared_matr, 2)
        throw(DomainError(MatriceSummary(squared_matr), "The given matrix isn't squared"))
    end
    #Mi accerto che la matrice sia di tipo Float64
    squared_matr = convert(Matrix{Float64}, squared_matr)

    n = size(squared_matr, 1)
    prod = 1.0
    for i in 1:n
        prod *= squared_matr[i, i]
    end
    
    return prod
end
#---------------------------------------------------------------------------------------------------------------------
#Calcola il determinante di una matrice n*n tringolare (superiore o inferiore è indifferente, il det è il prodotto degli elementi diagonali)
function tri_det(tri_matr)
    if !(tri_matr isa AbstractMatrix)
        throw(ArgumentError("The parameter must be an AbstractMatrix"))
    end
    #I want to be sure that the given matrix is squared
    if size(tri_matr, 1) != size(tri_matr, 2)
        throw(DomainError(MatriceSummary(tri_matr), "The given matrix isn't squared"))
    end
    #Mi accerto che la matrice sia di tipo Float64
    tri_matr = convert(Matrix{Float64}, tri_matr)

    det = 0
    #I want to be sure that the matrix given is triangular
    #istriu e istril restituiscono un booleano se la matrice è lower o upper triangular
    if istriu(tri_matr) == true || istril(tri_matr) == true
        det = diag_prod(tri_matr)
    
    #The matrix surely isn't either lower or upper triangular
    else
        throw(DomainError(MatriceSummary(tri_matr), "The given matrix isn't triangular. \n
        To perform the product of diagonal elements in a generic matrix use diag_prod() instead."))
    end

    return det
end
#---------------------------------------------------------------------------------------------------------------------
#Calcola il determinante di una matrice n*n
function my_det(A)
    if !(A isa AbstractMatrix)
        throw(ArgumentError("The parameter must be an AbstractMatrix"))
    end
    #Mi accerto che la matrice sia quadrata
    if size(A, 1) != size(A, 2)
        throw(DomainError(MatriceSummary(A), "The given matrix isn't squared"))
    end
    #Mi accerto che la matrice sia di tipo Float64
    A = convert(Matrix{Float64}, A)

    L, U = LU_dec(A)
    det = tri_det(L) * tri_det(U)
    return det
end
#---------------------------------------------------------------------------------------------------------------------
#Risolve un sistema lineare Ax = b
function solve_linear_system(A, b)
    if !(A isa AbstractMatrix)
        throw(ArgumentError("The first parameter must be an AbstractMatrix."))
    end
    if !(b isa AbstractVector)
        if b isa AbstractMatrix
            b = vec(b)
        else
            throw(ArgumentError("The second parameter must be an AbstractVector or an AbstractMatrix."))
        end    
    end
    #Mi accerto che la matrice e il vettore siano della stessa lunghezza
    if size(A, 1) != length(b)
        throw(DomainError(MatriceSummary(b), "Dimension mismatch between the matrix and the vector of known terms"))
    end
    #Mi accerto che la matrice e il vettore siano di tipo Float64
    A = convert(Matrix{Float64}, A)
    b = convert(Vector{Float64}, b)

    L, U = LU_dec(A)
    y = fw_sub(L, b)
    x = bw_sub(U, y)
    return x
end
#---------------------------------------------------------------------------------------------------------------------
#Risolve un sistema lineare Ax = b usando Cholesky decomposition
function sls_chlsky(A, b)
    if !(A isa AbstractMatrix)
        throw(ArgumentError("The first parameter must be an AbstractMatrix."))
    end
    if !(b isa AbstractVector)
        if b isa AbstractMatrix
            b = vec(b)
        else
            throw(ArgumentError("The second parameter must be an AbstractVector or an AbstractMatrix."))
        end    
    end
    #Mi accerto che la matrice e il vettore siano della stessa lunghezza
    if size(A, 1) != length(b)
        throw(DomainError(MatriceSummary(b), "Dimension mismatch between the matrix and the vector of known terms"))
    end
    #Mi accerto che la matrice e il vettore siano di tipo Float64
    A = convert(Matrix{Float64}, A)
    b = convert(Vector{Float64}, b)

    R = chlsky_dec(A)
    Rt = transpose(R)
    y = fw_sub(Rt, b)
    x = bw_sub(R, y)
    return x
end
#---------------------------------------------------------------------------------------------------------------------
#Calcola i coefficienti di un fit con metodo dei minimi quadrati
#la matrice A è m*n, b è di dimensione n
function least_sq(A, b)
    if !(A isa AbstractMatrix)
        throw(ArgumentError("The first parameter must be an AbstractMatrix."))
    end
    if !(b isa AbstractVector)
        if b isa AbstractMatrix
            b = vec(b)
        else
            throw(ArgumentError("The second parameter must be an AbstractVector or an AbstractMatrix."))
        end    
    end
    
    #Mi accerto che la matrice e il vettore siano di tipo Float64
    A = convert(Matrix{Float64}, A)
    b = convert(Vector{Float64}, b)
    
    m = size(A, 1)
    n = size(A, 2)
    m_b = length(b)

    if m != m_b
        throw(DomainError(MatriceSummary(A), "The system is not well defined. The matrix's size is $m * $n while the vector is $m_b long."))
    end

    At = transpose(A)
    N = At*A   
    z = At*b
    x = sls_chlsky(N, z)
    return x
end

#---------------------------------------------------------------------------------------------------------------------
println("linear_systems.jl loaded correctly")
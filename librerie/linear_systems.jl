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
    #Mi accerto che la matrice e il vettore siano di tipo Float64
    matrix = convert(Matrix{Float64}, matrix)
    b = convert(Matrix{Float64}, b)

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

    #Mi accerto che la matrice e il vettore siano di tipo Float64
    matrix = convert(Matrix{Float64}, matrix)
    b = convert(Matrix{Float64}, b)

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
#---------------------------------------------------------------------------------------------------------------------

#Funzione per il prodotto esterno. Restituisce una matrice a partire da due vettori
function outer_product(v, w)
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

#Esegue la decomposizione di una matrice in due matrici, una triangolare inferiore e l'altra triangolare superiore
#NON INCLUDE ROW-PIVOTING
function LU_dec(A)
    if size(A,1) != size(A,2)
        throw(DomainError(MatriceSummary(A), "The given matrix isn't squared."))
    end

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

#Calcola il prodotto degli elementi diagonali di una matrice n*n 
function diag_prod(squared_matr)
    #I want to be sure that the given matrix is squared
    if size(squared_matr, 1) != size(squared_matr, 2)
        throw(DomainError(MatriceSummary(squared_matr), "The given matrix isn't squared"))
    end
    
    n = size(squared_matr, 1)
    prod = 1.0
    for i in 1:n
        prod *= squared_matr[i, i]
    end
    
    return prod
end

#Calcola il determinante di una matrice n*n tringolare (superiore o inferiore è indifferente, il det è il prodotto degli elementi diagonali)
function tri_det(tri_matr)
    #I want to be sure that the given matrix is squared
    if size(tri_matr, 1) != size(tri_matr, 2)
        throw(DomainError(MatriceSummary(tri_matr), "The given matrix isn't squared"))
    end

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
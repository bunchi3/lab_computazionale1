#Import section
using LinearAlgebra

#This function operates a forward substitution for lower diagonal matrices (vector of vectors and matrices)
function fw_sub(matrix::Vector{Vector{Int64}}, b)
    n = length(matrix)
    if n != length(b)
        throw(DomainError(b, "Dimension mismatch between the matrix and the vector of known terms"))
    end
    solution = zeros(n)
    #i indice di riga, j indice di colonna
    for i in 1:n
        solution[i] = b[i]/matrix[i][i]
        for j in 1:i-1
            solution[i] += -matrix[i][j]*solution[j]/matrix[i][i] 
        end
    end
    return solution
end

function fw_sub(matrix::Vector{Vector{Float32}}, b)
    n = length(matrix)
    if n != length(b)
        throw(DomainError(b, "Dimension mismatch between the matrix and the vector of known terms"))
    end
    solution = zeros(n)
    #i indice di riga, j indice di colonna
    for i in 1:n
        solution[i] = b[i]/matrix[i][i]
        for j in 1:i-1
            solution[i] += -matrix[i][j]*solution[j]/matrix[i][i] 
        end
    end
    return solution
end

function fw_sub(matrix::Vector{Vector{Float64}}, b)
    n = length(matrix)
    if n != length(b)
        throw(DomainError(b, "Dimension mismatch between the matrix and the vector of known terms"))
    end
    solution = zeros(n)
    #i indice di riga, j indice di colonna
    for i in 1:n
        solution[i] = b[i]/matrix[i][i]
        for j in 1:i-1
            solution[i] += -matrix[i][j]*solution[j]/matrix[i][i] 
        end
    end
    return solution
end

function fw_sub(matrix::Matrix{Int64}, b)
    n = size(matrix, 1)
    if n != length(b)
        throw(DomainError(b, "Dimension mismatch between the matrix and the vector of known terms"))
    end
    solution = zeros(n)
    #i indice di riga, j indice di colonna
    for i in 1:n
        solution[i] = b[i]/matrix[i, i]
        for j in 1:i-1
            solution[i] += -matrix[i, j]*solution[j]/matrix[i, i] 
        end
    end
    return solution
end

function fw_sub(matrix::Matrix{Float32}, b)
    n = size(matrix, 1)
    if n != length(b)
        throw(DomainError(b, "Dimension mismatch between the matrix and the vector of known terms"))
    end
    solution = zeros(n)
    #i indice di riga, j indice di colonna
    for i in 1:n
        solution[i] = b[i]/matrix[i, i]
        for j in 1:i-1
            solution[i] += -matrix[i, j]*solution[j]/matrix[i, i] 
        end
    end
    return solution
end

function fw_sub(matrix::Matrix{Float64}, b)
    n = size(matrix, 1)
    if n != length(b)
        throw(DomainError(b, "Dimension mismatch between the matrix and the vector of known terms"))
    end
    solution = zeros(n)
    #i indice di riga, j indice di colonna
    for i in 1:n
        solution[i] = b[i]/matrix[i, i]
        for j in 1:i-1
            solution[i] += -matrix[i, j]*solution[j]/matrix[i, i] 
        end
    end
    return solution
end

#This function operates a backward substitution for upper diagonal matrices (vector of vectors and matrices)
function bw_sub(matrix::Vector{Vector{Int64}}, b)
    n = size(matrix, 1)
    if n != length(b)
        throw(DomainError(b, "Dimension mismatch between the matrix and the vector of known terms"))
    end
    solution = zeros(n)
    #i indice di riga, j indice di colonna
    for i in n:-1:1
        solution[i] = b[i]/matrix[i][i]
        for j in i+1:n
            solution[i] += -matrix[i][j]*solution[j]/matrix[i][i] 
        end
    end
    return solution
end

function bw_sub(matrix::Vector{Vector{Float32}}, b)
    n = size(matrix, 1)
    if n != length(b)
        throw(DomainError(b, "Dimension mismatch between the matrix and the vector of known terms"))
    end
    solution = zeros(n)
    #i indice di riga, j indice di colonna
    for i in n:-1:1
        solution[i] = b[i]/matrix[i][i]
        for j in i+1:n
            solution[i] += -matrix[i][j]*solution[j]/matrix[i][i] 
        end
    end
    return solution
end

function bw_sub(matrix::Vector{Vector{Float64}}, b)
    n = size(matrix, 1)
    if n != length(b)
        throw(DomainError(b, "Dimension mismatch between the matrix and the vector of known terms"))
    end
    solution = zeros(n)
    #i indice di riga, j indice di colonna
    for i in n:-1:1
        solution[i] = b[i]/matrix[i][i]
        for j in i+1:n
            solution[i] += -matrix[i][j]*solution[j]/matrix[i][i] 
        end
    end
    return solution
end

function bw_sub(matrix::Matrix{Int64}, b)
    n = size(matrix, 1)
    if n != length(b)
        throw(DomainError(b, "Dimension mismatch between the matrix and the vector of known terms"))
    end
    solution = zeros(n)
    #i indice di riga, j indice di colonna
    for i in n:-1:1
        solution[i] = b[i]/matrix[i, i]
        for j in i+1:n
            solution[i] += -matrix[i, j]*solution[j]/matrix[i, i] 
        end
    end
    return solution
end

function bw_sub(matrix::Matrix{Float32}, b)
    n = size(matrix, 1)
    if n != length(b)
        throw(DomainError(b, "Dimension mismatch between the matrix and the vector of known terms"))
    end
    solution = zeros(n)
    #i indice di riga, j indice di colonna
    for i in n:-1:1
        solution[i] = b[i]/matrix[i, i]
        for j in i+1:n
            solution[i] += -matrix[i, j]*solution[j]/matrix[i, i] 
        end
    end
    return solution
end

function bw_sub(matrix::Matrix{Float64}, b)
    n = size(matrix, 1)
    if n != length(b)
        throw(DomainError(b, "Dimension mismatch between the matrix and the vector of known terms"))
    end
    solution = zeros(n)
    #i indice di riga, j indice di colonna
    for i in n:-1:1
        solution[i] = b[i]/matrix[i, i]
        for j in i+1:n
            solution[i] += -matrix[i, j]*solution[j]/matrix[i, i] 
        end
    end
    return solution
end

#Funzione per il prodotto esterno. Restituisce una matrice a partire da due vettori
function outer_product(v, w)
    if length(v) != length(w)
        throw(DomainError(v, "Dimension mismatch between the two arrays"))
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

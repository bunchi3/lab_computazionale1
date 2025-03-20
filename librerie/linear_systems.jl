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


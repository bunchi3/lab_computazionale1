#File's index:
#1. IMPORT SECTION
#2. EULER'S METHOD
#3. RUNGE-KUTTA METHOD
#-------------------------------------------------------------------------------------------------------------------------------
#IMPORT SECTION

#-------------------------------------------------------------------------------------------------------------------------------
#EULER'S METHOD
function euler_ode(f::Function, u0::Number, n::Number, a::Number, b::Number)
    a, b = min(a, b), max(a, b)
    h = (b - a)/n
    u = Vector{Number}(undef, n+1)
    u[1] = u0
    
    for i in 1:1:(n)
        ti = a + h*(i-1)
        u[i+1] = u[i] + h*f(ti, u[i])
    end

    return u
end

function euler_ode(f::Vector, u0::Vector, n::Int, a::Number, b::Number)
    m = length(f)
    if m != length(u0)
        throw(DomainError(u0, "euler_ode says: the two parameters must be of the same length."))
    end

    h = (b-a)/n
    u = [Vector{Number}(undef, n+1) for _ in 1:1:m]

    #j will be used as index for the diff_equations
    #i will be used as index for the points of the single equation
    for j in 1:1:m
        u[j][1] = u0[j]        
    end

    for i in 1:1:n
        ti = a + (i-1)*h  
        #Defining a vector containing all the entries for the known function      
        v = [u[j][i] for j in 1:1:m]
        for j in 1:1:m
            u[j][i+1] = u[j][i] + h*f[j](ti, v...)             
        end
    end

    return u
end
#-------------------------------------------------------------------------------------------------------------------------------
#RUNGE-KUTTA METHOD
function IE2(f::Vector, u0::Vector, n::Int, a::Number, b::Number)
    m = length(f)
    if m != length(u0)
        throw(DomainError(u0, "euler_ode says: the two parameters must be of the same length."))
    end

    h = (b-a)/n
    u = [Vector{Number}(undef, n+1) for _ in 1:1:m]

    #j will be used as index for the diff_equations
    #i will be used as index for the points of the single equation
    for j in 1:1:m
        u[j][1] = u0[j]        
    end

    for i in 1:1:n
        ti = a + (i-1)*h  
        #Defining a vector containing all the entries for the known function      
        v = [u[j][i] for j in 1:1:m]
        k1 = [h*f[j](ti, v...) for j in 1:1:m]
        v_stage = [u[j][i] + 0.5*k1[j] for j in 1:1:m]  
        k2 = [h*f[j](ti + 0.5*h, v_stage...) for j in 1:1:m]   
        for j in 1:1:m         
            u[j][i+1] = u[j][i] + k2[j]    
        end
    end

    return u
end

function RK4(f::Vector, u0::Vector, n::Int, a::Number, b::Number)
    m = length(f)
    if m != length(u0)
        throw(DomainError(u0, "euler_ode says: the two parameters must be of the same length."))
    end

    h = (b-a)/n
    u = [Vector{Number}(undef, n+1) for _ in 1:1:m]

    #j will be used as index for the diff_equations
    #i will be used as index for the points of the single equation
    for j in 1:1:m
        u[j][1] = u0[j]        
    end

    for i in 1:1:n
        ti = a + (i-1)*h  
        #Defining a vector containing all the entries for the known function      
        v = [u[j][i] for j in 1:1:m]
        k1 = [h*f[j](ti, v...) for j in 1:1:m]
        v_stage1 = [u[j][i] + 0.5*k1[j] for j in 1:1:m]  
        k2 = [h*f[j](ti + 0.5*h, v_stage1...) for j in 1:1:m]   
        v_stage2 = [u[j][i] + 0.5*k2[j] for j in 1:1:m]
        k3 = [h*f[j](ti + 0.5*h, v_stage2...) for j in 1:1:m]
        v_stage3 = [u[j][i] + k3[j] for j in 1:1:m]
        k4 = [h*f[j](ti + h, v_stage3...) for j in 1:1:m]
        for j in 1:1:m         
            u[j][i+1] = u[j][i] + (k1[j] + 2*k2[j] + 2*k3[j] + k4[j]) / 6
        end
    end

    return u
end
#-------------------------------------------------------------------------------------------------------------------------------
println("diff_equations.jl loaded correctly")
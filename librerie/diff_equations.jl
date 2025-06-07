#File's index:
#1. IMPORT SECTION
#2. EULER'S METHOD
#3. RUNGE-KUTTA METHOD
#-------------------------------------------------------------------------------------------------------------------------------
#IMPORT SECTION
#-------------------------------------------------------------------------------------------------------------------------------
#STEP SIZE FINDING
function h_bool(u1::Vector, u2::Vector, err_tol::Number)
    if length(u1) != length(u2)
        throw(DomainError(u1, "h_bool says: the two parameters must be of the same length."))
    end
    norm_inf = maximum(abs.(u1 .- u2))
    if norm_inf > err_tol
        println("Norma infinito: ", norm_inf)
        return false
    else
        println("Norma infinito: ", norm_inf)
        return true
    end
end

function h_search(ODE_func1::Function, ODE_func2::Function, f::Vector, u0::Vector, n::Number, a::Number, b::Number; err_tol::Number = 1e-6)
    m = length(f)
    if m != length(u0)
        throw(DomainError(u0, "h_search says: f and u must be of the same length."))
    end
    stability = fill(false, m)
    iter = 1
    u1 = ODE_func1(f, u0, n, a, b)
    while all(!x for x in stability) && iter < 10
        u2 = ODE_func2(f, u0, 2^(iter)*n, a, b)
        for j in 1:1:m
            u2_comparison = [u2[j][i] for i in 1:2:length(u2[j])]
            println("$j ° solution: ")
            stability[j] = h_bool(u1[j], u2_comparison, err_tol)
        end
        
        iter += 1
        u1 = u2         #I use u2 as the new u1 for the next iteration
    end
    println("Iterations: ", iter)
    if all(!x for x in stability) && iter == 10
        throw(DomainError(n, "h_search says: the two methods are not stable."))
    end
    return abs((b - a)/(2^(iter)*n))
end
#-------------------------------------------------------------------------------------------------------------------------------
#SOLUTION LONG STEADY STATE
function l_t_steady(t::Vector, u::Vector; err_tol::Number = 10^-9)
    n = length(t)
    if n != length(u)
        throw(DomainError(u, "l_t_steady says: the two parameters must be of the same length."))
    end
    if n < 2
        throw(DomainError(t, "l_t_steady says: the two parameters must have at least two elements."))
    end
    #Check if the time vector is sorted
    if !issorted(t)
        throw(DomainError(t, "l_t_steady says: the time vector must be sorted."))
    end

    for i in Int(floor(n/2)):1:n-1
        if abs(u[i] - u[i+1]) < err_tol
            return t[i]
        else
            throw(DomainError(t, "l_t_steady says: the steady state was not reached."))
        end
    end
end
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
    u = [Vector{Number}(undef, n+1) for _ in 1:1:m]     #n+1 because we need to include the initial point

    #j will be used as index for the diff_equations
    #i will be used as index for the points of the single equation
    #for loop for the initial point
    for j in 1:1:m
        u[j][1] = u0[j]        
    end

    #for loop for the rest of the points, starting from the first one, excluding the known value (u0)
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

    #u contains m vectors, each of length n+1
    return u
end
#-------------------------------------------------------------------------------------------------------------------------------
println("diff_equations.jl loaded correctly")
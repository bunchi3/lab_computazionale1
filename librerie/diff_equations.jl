#File's index:
#1. IMPORT SECTION
#2. STEP SIZE FINDING
#3. SOLUTION LONG STEADY STATE
#4. EULER'S METHOD
#5. RUNGE-KUTTA METHOD
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

#This function computes the next solution point given h, for two different methods, BS2 and BS3
function BS23(f::Vector, ui::Vector, ti::Number, h::Float64)
    m = length(f)
    if m != length(ui)
        throw(DomainError(ui, "euler_ode says: the two parameters must be of the same length."))
    end

    u2, u3 = zeros(m), zeros(m)     #2 and 3 because of the solution order

    #j will be used as index for the diff_equations  
    #Defining a vector containing all the entries for the known function      
    v = ui
    k1 = [h*f[j](ti, v...) for j in 1:m]
    k2 = [h*f[j](ti + 0.5*h, (ui .+ 0.5 .* k1)...) for j in 1:m]
    k3 = [h*f[j](ti + 0.75*h, (ui .+ 0.75 .* k2)...) for j in 1:m]
    k4 = [h*f[j](ti + h, (ui .+ (2/9).*k1 .+ (1/3).*k2 .+ (4/9).*k3)...) for j in 1:m]
    for j in 1:1:m         
        u3[j] = ui[j] + (2.0/9.0)*k1[j] + (1.0/3.0)*k2[j] + (4.0/9.0)*k3[j]
    end
    for j in 1:1:m    
        u2[j] = ui[j] + (7/24)*k1[j] + (1/4)*k2[j] + (1/3)*k3[j] + (1/8)*k4[j]
    end 

    #u3, u2 contains m vectors, each of length n+1
    return u3, u2
end

function adaptive_bs23(f::Vector, ui::Vector, ti::Number, h::Float64, δ::Float64, k1_prec = nothing)
    m = length(f)
    k4 = []
    ti_1 = 0.0
    if m != length(ui)
        throw(DomainError(ui, "euler_ode says: the two parameters must be of the same length."))
    end
    q = 1.0
    H = q*h
    u_num2 = zeros(m)  #u_num2 is u_{i+1} for the second order method
    u_num3 = zeros(m)
    iter = 0
    v = ui
    #if k1_prec is nothing k1 will be the first expression, else will be k1_prec
    k1 = isnothing(k1_prec) ? [H*f[j](ti, v...) for j in 1:m] : k1_prec
    while true 
        #u_num is u_{i+1}
        k2 = [H*f[j](ti + 0.5*H, (ui .+ 0.5 .* k1)...) for j in 1:m]
        k3 = [H*f[j](ti + 0.75*H, (ui .+ 0.75 .* k2)...) for j in 1:m]
        k4 = [H*f[j](ti + H, (ui .+ (2.0/9.0).*k1 .+ (1.0/3.0).*k2 .+ (4.0/9.0).*k3)...) for j in 1:m]
        
        for j in 1:1:m         
            u_num3[j] = ui[j] + (2.0/9.0)*k1[j] + (1.0/3.0)*k2[j] + (4.0/9.0)*k3[j]
        end
        for j in 1:1:m    
            u_num2[j] = ui[j] + (7.0/24.0)*k1[j] + (1.0/4.0)*k2[j] + (1.0/3.0)*k3[j] + (1.0/8.0)*k4[j]
        end
        
        ϵi = δ*(1 + maximum(abs.(ui)))
        Ei = maximum(abs.(u_num2 .- u_num3))
        
        if Ei < ϵi
            ti_1 = ti+H
            q = 0.8 * (ϵi/Ei)^(1.0/3.0)
            q = min(q, 4)
            f4 = k4/H
            H = q*H         #this is the H which will be the guess for the next step 
            k4 = H*f4       #this means that k4 won't be exactly the k1 for the next step, but needs to be re-escalate
            return u_num3, ti_1, H, k4
        end

        q = 0.8 * (ϵi/Ei)^(1.0/3.0)
        q = min(q, 4)
        H = q*H
        k1 = [H*f[j](ti, v...) for j in 1:m]
        
        iter += 1
        if iter > 1000
            println("Max iterations reached")
            break
        end
    end
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

#= function BS2(f::Vector, u0::Vector, n::Int, a::Number, b::Number)
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
        v_stage2 = [u[j][i] + (3.0/4.0)*k2[j] for j in 1:1:m]
        k3 = [h*f[j](ti + (3.0/4.0)*h, v_stage2...) for j in 1:1:m]
        v_stage3 = [u[j][i] + (2.0/9.0)*k1[j] + (1.0/3.0)*k2[j] + (4.0/9.0)*k3[j] for j in 1:1:m]
        k4 = [h*f[j](ti + h, v_stage3...) for j in 1:1:m]
        for j in 1:1:m         
            u[j][i+1] = u[j][i] + ((7/24)*k1[j] + (1/4)*k2[j] + (1/3)*k3[j] + (1/8)*k4[j])
        end
    end

    #u contains m vectors, each of length n+1
    return u
end

function BS3(f::Vector, u0::Vector, n::Int, a::Number, b::Number)
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
        v_stage2 = [u[j][i] + (3.0/4.0)*k2[j] for j in 1:1:m]
        k3 = [h*f[j](ti + (3.0/4.0)*h, v_stage2...) for j in 1:1:m]
        v_stage3 = [u[j][i] + (2.0/9.0)*k1[j] + (1.0/3.0)*k2[j] + (4.0/9.0)*k3[j] for j in 1:1:m]
        #k4 = [h*f[j](ti + h, v_stage3...) for j in 1:1:m]
        #first same as last
        for j in 1:1:m         
            u[j][i+1] = v_stage3[j]
        end
    end

    #u contains m vectors, each of length n+1
    return u
end =#
#-------------------------------------------------------------------------------------------------------------------------------
#The key in understanding the process was to separate h of the step from the h of the guess
function solution_adapt_bs23(f::Vector, u0::Vector, a::Number, b::Number, δ::Float64)
    m = length(f)
    u = [[] for _ in 1:1:m]     #n+1 because we need to include the initial point
    t = [a]
    h = []
    h_guess = δ^(1.0/3.0) / 2.0
    #i index of the solution Number
    #j index of the equations
    for j in 1:1:m
        push!(u[j],  u0[j])        
    end

    k1 = nothing

    #ui_1 is u_{i+1} and hi1 is the time step to reach t_{i+1}
    while t[end] < b
        ui_1, ti_1, h_guess, k1 = adaptive_bs23(f, [u[j][end] for j in 1:1:m], t[end], h_guess, δ, k1) 
        
        ti = t[end]
        if ti_1 > b
            for j in 1:1:m
                u_last_forced = BS23(f, [u[j][end] for j in 1:1:m], t[end], (b-ti))[1]
                push!(u[j], u_last_forced[j])
                
            end
            push!(t, b)
            push!(h, b-ti)
            return t, u, h
        end
        push!(h, ti_1 - ti) 
        #Interrupting the routine if h doesn't make the time grow.
        if ti == ti + h_guess
            return t, u, h
        end
        
        for j in 1:1:m
            push!(u[j], ui_1[j])            
        end  
        push!(t, ti_1)

    end   
end
#-------------------------------------------------------------------------------------------------------------------------------
println("diff_equations.jl loaded correctly")
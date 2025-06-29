#File's index
# 1.  #IMPORT SECTION
# 2.  #NEWTON-COTES INTEGRATION RULES
# 3.  #GAUSS-LEGENDRE QUADRATURE FIRST VERSION
# 4.  #GAUSS-LEGENDRE QUADRATURE
# 5.  #CLENSHAW-CURTIS QUADRATURE
# 6.  #DOUBLE EXPONENTIAL QUADRATURE
#-------------------------------------------------------------------------------------------------------------------------------
#IMPORT SECTION
using Dates
include("c:\\ALL\\Stefano\\Bicocca\\3terzo_anno\\lab_comp\\lab_computazionale1\\librerie\\non_linear_roots.jl")
include("c:\\ALL\\Stefano\\Bicocca\\3terzo_anno\\lab_comp\\lab_computazionale1\\librerie\\polynomials.jl")
#-------------------------------------------------------------------------------------------------------------------------------
#NEWTON-COTES INTEGRATION RULES
#=   
-IntegralTrap: returns the integral with the trapezoidal rule.
               Requires the function f, the integration interval extrema (a, b), and the number of intervals 
-IntegralSim: returns the integral with Simpson's rule.
              Requires the function f, the integration interval extrema (a, b), and the number of intervals 
=#

function IntegralTrap(f::Function, a::Number, b::Number, m::Integer)
    if a>b
        shift = b
        b = a 
        a = shift        
    end
    if m <= 0 
        throw(DomainError(m, "IntegralTrap says: the last parameter is the intervals' number, can't be negative or zero"))
    end

    h = (b-a)/m
    sum = 0.0

    for j in 1:1:(m-1)
        sum+=f(a+j*h)
    end
    
    return h/2 * (f(a) + f(b) + 2*sum)
end 
function IntegralSim(f::Function, a::Number, b::Number, m::Integer)
    if a>b
        shift = b
        b = a 
        a = shift        
    end
    if m <= 0 
        throw(DomainError(m, "IntegralTrap says: the last parameter is the intervals' number, can't be negative or zero"))
    end

    h = (b-a) / (2*m)
    sum1 = 0.0
    sum2 = 0.0

    for j in 1:1:m
        sum1 += f(a + (2*j - 1)*h)        
    end
    for j in 1:1:(m-1)
        sum2 += f(a + 2*j*h)
    end

    return (h/3)*(f(a) + f(b) + 4*sum1 + 2*sum2)
end
#---------------------------------------------------------------------------------------------------------------------
#GAUSS-LEGENDRE QUADRATURE FIRST VERSION
#NB: this version actually works, but is too slow.
#= 
-w_slow(Int): returns the weights for every root of a legendre polynomial
         requires n, the polynomial degree
         even though it returns every weight, since the roots are symmetric it just computes half of them and then
         duplicates them.
-IntegralGaussLegendre_slow(Function, Int): returns the value of the integral between -1, 1
                                       requires f (the function to be integrated), n, the integration precision
=#
#= function newt_met_steps(f::Function, df::Function, x_start::Number)
    x_tol::Float64 = 1.0e-15
    y_tol::Float64 = 1.0e-15
    iter_max::Int = 1000
    iter::Int = 0
    x = Vector{Float64}(undef, iter_max)
    x[1] = x_start

    t2 = 0.0
    result = nothing
    t1 = @elapsed begin
        for i in 2:iter_max
            t2 += @elapsed x[i] = x[i-1] - f(x[i-1])/df(x[i-1])
            if abs(f(x[i])) < y_tol
                result = (x[i], x)
                break
            end
            if abs(x[i] - x[i-1]) < x_tol
                result = (x[i], x)
                break
            end
            if i >= iter_max
                println("newt_met_steps says: iterations' maximum number reached. Found zero could be inaccurate.")
                result = (x[i], x)
                break
            end
            iter += i
        end
    end
    #println(t1)
    println("newt_met_steps step: ", t2)
    println("Step number:         ", iter)
    return result
end

function w_slow(n::Int, pl::Vector{Function}, dpl::Vector{Function})
    P = pl[n+1]
    dP = dpl[n+1]

    #calculating legendre roots----------------------------------------
    roots = ones(n)
    if n%2 == 0
        t1 = @elapsed for k in 1:1:Int(n/2)
            x0_k = cos(pi*(4k-1)/(4n+2))
            r, r_steps = newt_met_steps(P, dP, x0_k)
            println("Root found:          ", r)
            println()
            roots[k] = r
            roots[end-k+1] = -r
        end
    else
        for k in 1:1:Int((n+1)/2)
            x0_k = cos(pi*(4k-1)/(4n+2))
            r, r_steps = newt_met_steps(P, dP, x0_k)
            roots[k] = r
            if k != length(roots)
                roots[end-k+1] = -r
            end
        end
    end
    
    #calculating legendre roots----------------------------------------
    
    #calculating the weights-------------------------------------------
    
    w = ones(n)
    
    if n%2 == 0
        t2 = @elapsed for k in 1:1:Int((length(roots))/2)
            xk = roots[k]
            wk = 2/((1-xk^2)*(dP(xk))^2)
            w[k] = wk
            w[end-k+1] = wk
        end
    
    else
        for k in 1:1:Int((length(roots)+1)/2)
            xk = roots[k]
            wk = 2/((1-xk^2)*(dP(xk))^2)
            if k == length(roots)
                w[k] = wk
            else
                w[k] = wk
                w[end-k+1] = wk
            end
        end
    end
    
    #calculating the weights-------------------------------------------
    println("Inside w: time sum of every call of newt_met_steps: ", t1)
    println("Inside w: time sum of calculating the weights: ", t2)
    return w, roots
end

function IntegralGaussLegendre_slow(f::Function, n::Int)
    t1 = @elapsed pl, dpl = all_leg_pol(n)
    t2 = @elapsed weight_arr, xk_arr = w_slow(n, pl, dpl)
    sum::Float64 = 0.0
    for i in 1:1:n
        sum += f(xk_arr[i])*weight_arr[i]
    end
    println("Function all_leg_pol: ", t1)
    println("Function w: ", t2)
    return sum
end =#
#---------------------------------------------------------------------------------------------------------------------
#GAUSS-LEGENDRE QUADRATURE
#= 
-w(Int): returns the weights for every root of a legendre polynomial
         requires n, the polynomial degree
         even though it returns every weight, since the roots are symmetric it just computes half of them and then
         duplicates them.
-IntegralGaussLegendre(Function, Int): returns the value of the integral between -1, 1
                                       requires f (the function to be integrated), n, the integration precision
=#
function wGL(n::Int)
    x_tol::Float64 = 1.0e-15
    y_tol::Float64 = 1.0e-15
    iter_max::Int = 1000
    roots = ones(n)
    weight = ones(n)

    for k in 1:1:Int(floor(n/2))
        iter::Int = 1
        xk::Float64 = cos(pi*(4k-1)/(4n+2))
        dx::Float64 = 1.0
        p0::Float64 = 1.0
        while abs(dx) > x_tol && abs(p0) > y_tol
            p0, p0der = legforintegration(n, xk) 
            dx = p0/p0der
            xk = xk - dx  
            iter += 1
            if iter > iter_max
                println("w says: iterations' maximum number reached. Found zero could be inaccurate.")
                if abs(p0) > y_tol
                    throw(DomainError(n, "w says: the function calculated in the estimated root is not close enoughto zero."))  
                elseif abs(dx > x_tol)
                    throw(DomainError(n, "w says: the root could not be found with the required precision."))
                else
                    println("w says: the root is: ", xk)
                end
            end 
        end
        roots[k] = xk
        roots[end-k+1] = -xk

        wk = 2/((1-xk^2)*(p0der)^2)
        weight[k] = wk
        weight[end-k+1] = wk
    end

    #Adding the middle root if n is odd
    if n%2 == 1
        xk = 0.0
        middle_index::Int = Int((n+1)/2)
        roots[middle_index] = xk
        p0, p0der = legforintegration(n, xk)
        weight[middle_index] = 2/((1-xk^2)*(p0der)^2)
    end

    return roots, weight
end

#n is the polynomial order
function IntegralGaussLegendre(f::Function, n::Int; a::Number = -1.0, b::Number = 1.0)
    a, b = min(a, b), max(a, b)
    sum::Float64 = 0.0
    norm::Float64 = (b-a)/2
    root, weight = wGL(n)
    
    g = x -> f(0.5*((b-a)*x + (a+b)))

    for k in 1:1:n
        sum += g(root[k])*weight[k]
    end
    return norm*sum
end

#---------------------------------------------------------------------------------------------------------------------
#CLENSHAW-CURTIS QUADRATURE
#= 
-w_CC_even: returns the weights for Clenshaw-Curtis integration, with even n
            requires n, the number of weights, and k, the weight's index
-IntegralCC_even: returns the value of the integral between -1, 1, having even n.
                  requires f (the function to be integrated), n, the number of weights.
=#

function w_CC_even(n::Int64, k::Int)
    if n%2 == 1
        throw(DomainError(n, "w_CC_even says: n must be even."))
    end

    wk = 0.0
    n1 = Int(n/2)
    b = [2 for j in 1:1:(n1-1)]
    push!(b, 1)
    ck = 2.0
    if k == 0 || k == n
        ck = 1.0        
    end
    
    sum = 0.0
    for j in 1:1:n1
        sum += b[j]*(cos(2*j*k*pi/n))/(4*j^2 - 1)
    end
    wk = (ck/n)*(1-sum)

    return wk
end

function IntegralCC_even(f::Function, n::Int; a::Number = -1.0, b::Number = 1.0)
    if n%2 == 1
        throw(DomainError(n, "IntegralCC_even says: n must be even."))
    end
    a, b = min(a, b), max(a, b)
    norm::Float64 = (b-a)/2
    sum = 0.0

    for k in 0:1:n
        xk = cos(k*pi/n)
        xk = 0.5*((b-a)*xk + (a+b))
        sum += w_CC_even(n, k)*f(xk)        
    end

    return norm*sum
end
#---------------------------------------------------------------------------------------------------------------------
#DOUBLE EXPONENTIAL QUADRATURE
#= 
-Integral_DE: has two methods
    1: returns the value of the integral of a function f in a finite interval, using the double exponential method
       requires the function f to be integrated; N which is half of the points of the polynomial approximation;
       a and b which are the interval boundaries.
    2: returns the value of the integral of a function f in a infinite interval and in [-1, 1], using the double 
       exponential method
       requires the function f to be integrated; N which is half of the points of the polynomial approximation;
       a and b which specify the boundaries (l is for left, r is for right, e.g.: lone = left one = -1 ); 
       e specifies if in the function there is e^-x. 
=#
function Integral_DE(f::Function, N::Int, a::Number, b::Number)
    a, b = min(a, b), max(a, b)
    norm = (b-a)/2.0

    g = x -> f(0.5*((b-a)*x + a + b))
    I = Integral_DE(g, N)
    
    return norm*I
end

function Integral_DE(f::Function, N::Int; a::Symbol=:lone, b::Symbol=:rone, e::Bool=false)
    if (a != :zero && a != :linf && a != :lone)
        throw(DomainError(a, "Integral_DE says: the first parameter must be one of the following: :zero, :linf, :rone"))
    end
    if (b != :rinf && b != :rone)
        throw(DomainError(b, "Integral_DE says: the second parameter must be the following: :rinf, :rone"))
    end
    if N <= 0
        throw(DomainError(N, "Integral_DE says: the last parameter is the intervals' number, can't be negative or zero"))
    end

    sum = 0.0
    
    if a == :lone && b == :rone
        Φ = t -> tanh(pi/2 * sinh(t))
        Φ_der = t -> (pi/2 * cosh(t))/(cosh(pi/2 * sinh(t)))^2    
    end

    if a == :zero && b == :rinf
        Φ = t -> exp(pi/2 * sinh(t))
        Φ_der = t -> pi/2 * exp(pi/2 * sinh(t)) * cosh(t)
    end

    if a == :zero && b == :rinf && e == true
        Φ = t -> exp(t - exp(-t))
        Φ_der = t -> exp(-exp(-t)) * (exp(t) + 1)
    end

    if a == :linf && b == :rinf
        Φ = t -> sinh((pi/2) * sinh(t))
        Φ_der = t -> pi/2 * cosh(t) * cosh(pi/2 * sinh(t))
    end

    g = t -> f(Φ(t)) * Φ_der(t)

    #I start with a big tM already (and not tM=0.0)
    tM = 2.0
    t_prec = 10.0^(-5)
    iter::Int = 0
    iter_max::Int = 10^6

    while (abs(g(-tM)) > eps() || abs(g(tM)) > eps()) && iter < iter_max
        tM +=  t_prec
        iter += 1
    end
    
    if iter >= iter_max
        display("Integral_DE says: finding tM was too long. Integral found could be inaccurate.")     
    end

    h = tM/N

    for k in -N:1:N
        sum += g(k*h)      
    end

    return h*sum
end
#---------------------------------------------------------------------------------------------------------------------
println("integration.jl loaded correctly")
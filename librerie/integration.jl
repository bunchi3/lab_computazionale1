#File's index
# 1.  #IMPORT SECTION
# 2.  #NEWTON-COTES INTEGRATION RULES
#-------------------------------------------------------------------------------------------------------------------------------
#IMPORT SECTION
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
#GAUSS-LEGENDRE QUADRATURE
#= 
-w: returns the weights for every root of a legendre polynomial
     requires n, the polynomial degree
     even though it returns every weight, since the roots are symmetric it just computes half of them and then
     duplicates them.
-IntegralGaussLegendre: returns the value of the integral between -1, 1
                        requires f (the function to be integrated), n, the integration precision
=#

function w(n::Int)
    pl, dpl = all_leg_pol(n) 

    P = pl[n+1]
    dP = dpl[n+1]

    #calculating legendre roots----------------------------------------
    roots = Float64[]
    if n%2 == 0
        for k in 1:1:(n/2)
            x0_k = cos(pi*(4k-1)/(4n+2))
            r, r_steps = newt_met_steps(P, dP, x0_k, 1)
            push!(roots, r)
        end
    else
        for k in 1:1:((n+1)/2)
            x0_k = cos(pi*(4k-1)/(4n+2))
            r, r_steps = newt_met_steps(P, dP, x0_k, 1)
            push!(roots, r)
        end
    end
    #calculating legendre roots----------------------------------------
    
    #calculating the weights-------------------------------------------
    w = ones(n)
    if n%2 == 0
        for k in 1:1:length(roots)
            xk = roots[k]
            wk = 2/((1-xk^2)*(dP(xk))^2)
            w[k] = wk
            w[end-k+1] = wk
        end
    else
        for k in 1:1:length(roots)
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
    
    return w
end

function IntegralGaussLegendre(f::Function, n::Int)
    weight_arr = w(n)
    xk_arr = leg_roots(n)
    sum = 0.0

    for i in 1:1:n
        sum += f(xk_arr[i])*weight_arr[i]
    end
    return sum
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

function IntegralCC_even(f::Function, n::Int)
    if n%2 == 1
        throw(DomainError(n, "IntegralCC_even says: n must be even."))
    end

    sum = 0.0
    for k in 0:1:n
        xk = cos(k*pi/n)
        sum += w_CC_even(n, k)*f(xk)        
    end

    return sum
end
#---------------------------------------------------------------------------------------------------------------------
println("integration.jl loaded correctly")
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
-wk: returns the weight for a certain root of a legendre polynomial
     requires n, the polynomial degree
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
#---------------------------------------------------------------------------------------------------------------------
println("integration.jl loaded correctly")
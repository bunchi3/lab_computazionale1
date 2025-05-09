#File's index
# 1.  #IMPORT SECTION
# 2.  #NEWTON-COTES INTEGRATION RULES
#-------------------------------------------------------------------------------------------------------------------------------
#IMPORT SECTION

#-------------------------------------------------------------------------------------------------------------------------------
#NEWTON-COTES INTEGRATION RULES
#=  
-IntegralTrap calculates the integral with the trapezoidal rule.
Requires the function f, the integration interval extrema (a, b), and the number of intervals 
-IntegralSim calculates the integral with Simpson's rule.
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
println("integration.jl loaded correctly")
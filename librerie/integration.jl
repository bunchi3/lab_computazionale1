#File's index
# 
#-------------------------------------------------------------------------------------------------------------------------------
#import section

#-------------------------------------------------------------------------------------------------------------------------------
#=  
IntegralTrap calculates the integral with the trapezoidal rule.
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

#---------------------------------------------------------------------------------------------------------------------
println("integration.jl loaded correctly")
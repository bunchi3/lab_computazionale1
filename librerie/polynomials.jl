#File's index
# 1.  #IMPORT SECTION
# 2.  #LEGENDRE POLYNOMIALS
#-------------------------------------------------------------------------------------------------------------------------------
#IMPORT SECTION

#-------------------------------------------------------------------------------------------------------------------------------
#LEGENDRE POLYNOMIALS
#= 
-leg_pol: returns a vector containing all legendre polynomials as functions, up to n-th order
          requires n, which is the maximum polynomial order
          useful when you just need the polynomials
-leg_pol_der: returns a vector containing all first derivatives of legendre polynomials as functions, up to n-th order
              requires n, which is the maximum polynomial order
              useful when you just need the first derivatives' polynomials
-all_leg_pol: returns two vector containing both legendre polynomials and their derivative, as functions, up to n-th order
              requires n, which is the maximum polynomial order
              useful when you need the everything
=#

function leg_pol(n::Int)
    #this is a vector which will contain the polynomial in function format
    pol_vec = [x -> 1, x -> x]
    
    #= for m we calcualate order P_{m+1}, this means that for m = 1 we calculate P_2 (we already have P_0 and P_1).
    We want P_n, so m must be n-1 at best.  
    Moreover, Julia language starts numerating arrays from 1, so every polynomial is shifted by 1 in position and degree=#
    for m in 1:1:(n-1) 
        func = x -> ((2m+1)*x*pol_vec[m+1](x) - m*pol_vec[m](x))/(m+1)
        push!(pol_vec, func)
    end
    
    return pol_vec
end
function leg_pol_der(n::Int)
    pol = leg_pol(n)
    der_vec = [x -> 0, x -> 1]

    for m in 2:1:n
        func = x -> m*(x*pol[m+1](x) - pol[m](x))/(x^2-1)    
        push!(der_vec, func)     
    end

    return der_vec
end
function all_leg_pol(n::Int)
    pol_vec = [x -> 1, x -> x]
    der_vec = [x -> 0, x -> 1]

    for m in 1:1:(n-1) 
        func = x -> ((2m+1)*x*pol_vec[m+1](x) - m*pol_vec[m](x))/(m+1)
        push!(pol_vec, func)
    end

    for m in 2:1:n
        func = x -> m*(x*pol_vec[m+1](x) - pol_vec[m](x))/(x^2-1)    
        push!(der_vec, func)     
    end

    return pol_vec, der_vec
end
#---------------------------------------------------------------------------------------------------------------------
println("polynomials.jl loaded correctly")
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
              useful when you need both legendre polynomials and their derivative.
-leg_roots: return a vector containing the zeros of legendre polynomial of n-th order.
            requires n, which is the polynomial degree.
            just return the positive zeros, the others are symmetric.

Please note: using a certain n will return a vector of lenght n+1 because of polynomial of 0-th order. Polynomial of
             n-th order will be stored in n+1 position because of Julia arrays ordering
-legforintegration: returns two numbers, the value of legendre polynomial of n-th order evaluated in a specific x,
                    and the value of its derivative evaluated in the same x.
                    requires n, which is the polynomial degree, and the point in which evaluate the polynomial and 
                    its derivative.
=#

function leg_pol(n::Int)
    #this is a vector which will contain the polynomial in function format
    pol_vec = Function[x -> 1, x -> x]
    
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
    der_vec = Function[x -> 0, x -> 1]

    for m in 2:1:n
        func = x -> m*(x*pol[m+1](x) - pol[m](x))/(x^2-1)    
        push!(der_vec, func)     
    end

    return der_vec
end
function all_leg_pol(n::Int)
    pol_vec = Function[x -> 1, x -> x]
    der_vec = Function[x -> 0, x -> 1]

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
function leg_roots(n::Int, pl::Vector{Function}, dpl::Vector{Function})
    P = pl[n+1]
    dP = dpl[n+1]

    roots = ones(n)
    
    if n%2 == 0
        for k in 1:1:Int(n/2)
            x0_k = cos(pi*(4k-1)/(4n+2))
            r, r_steps = newt_met_steps(P, dP, x0_k, 1)
            roots[k] = r
            roots[end-k+1] = -r
        end
    else
        for k in 1:1:Int((n+1)/2)
            x0_k = cos(pi*(4k-1)/(4n+2))
            r, r_steps = newt_met_steps(P, dP, x0_k, 1)
            roots[k] = r
            if k != length(roots)
                roots[end-k+1] = -r
            end
        end
    end
    return roots
end
function legforintegration(n::Int, x::Float64)
    p0 = 1
    p1 = x
    ptemp = 0.0
    pder = 0.0

    for m in 1:1:(n-1)
        ptemp = p1
        p1 = ((2m+1)*x*ptemp - m*p0)/(m+1)
        p0 = ptemp
    end

    pder = n*(x*p1 - p0)/(x^2-1)

    return p1, pder
end
#---------------------------------------------------------------------------------------------------------------------
println("polynomials.jl loaded correctly")
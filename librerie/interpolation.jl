#File's index:
# 1. POLYNOMIAL FIT  WITH LEAST SQUARES METHOD
# 2. FOURIER FIT  WITH LEAST SQUARES METHOD
# 3. CHANGE OF VARIABLES
# 4. SUPPORT COEFFICIENTS FOR LAGRANGIAN INTERPOLATION
# 5. LAGRANGIAN INTERPOLATION
# 6. ERROR INDICATION FUNCTION

#-------------------------------------------------------------------------------------------------------------------------------
#import section
using Markdown     #markdown visualization in display
using LaTeXStrings
include("c:\\ALL\\Stefano\\Bicocca\\3terzo_anno\\lab_comp\\lab_computazionale1\\librerie\\linear_systems.jl")
#-------------------------------------------------------------------------------------------------------------------------------
#POLYNOMIAL FIT  WITH LEAST SQUARES METHOD
#Polynomial series. Requires x and the series coefficients
#Returns an array with the polynomial evaluated at the x values
function pol_expansion(x, c)
    # Verifica che x e c siano vettori
    if !(x isa AbstractVector)
        throw(ArgumentError("Il parametro x deve essere un AbstractVector, ma è di tipo $(typeof(x))"))
    end
    if !(c isa AbstractVector)
        throw(ArgumentError("Il parametro c deve essere un AbstractVector, ma è di tipo $(typeof(c))"))
    end
    # Verifica che la lunghezza di c sia minore o uguale alla lunghezza di x
    if length(c) > length(x)
        throw(ArgumentError("La lunghezza dei coefficienti deve essere minore o uguale alla lunghezza dei valori in x"))
    end

    sum = zeros(length(x))
    for j in 1:1:length(x)
        for i in 1:1:length(c)
            sum[j] += c[i]*(x[j]^(i-1))
        end
    end
    return sum
end
#Polynomial fit. Requires x data, the y data, the number of coefficients you want from the fit
function pol_fit(x, y, n_coeff)
    #Sezione di test iniziali
    # Verifica che x e y siano vettori
    if !(x isa AbstractVector)
        throw(ArgumentError("Il parametro x deve essere un AbstractVector, ma è di tipo $(typeof(x))"))
    end
    if !(y isa AbstractVector)
        throw(ArgumentError("Il parametro c deve essere un AbstractVector, ma è di tipo $(typeof(c))"))
    end
    # Verifica che la lunghezza di x e y siano uguali
    if length(x) != length(y)
        throw(ArgumentError("La lunghezza di x e y devono essere uguali, ma sono $(length(x)) e $(length(y))"))
    end
    # Verifica che n_coeff sia un numero intero positivo
    if !(n_coeff isa Integer) || n_coeff <= 0
        throw(ArgumentError("Il numero di coefficienti deve essere un intero positivo, ma è $(n_coeff)"))
    end
    # Verifica che n_coeff sia minore o uguale alla lunghezza di x
    if n_coeff > length(x)
        throw(ArgumentError("Il numero di coefficienti deve essere minore o uguale alla lunghezza di x, ma sono $(n_coeff) e $(length(x))"))
    end

    #Creating the design matrix (in italiano matrice dei regressori)
    A = hcat([x.^(i-1) for i in 1:1:n_coeff]...) 
    c = least_sq(A, y)
    return c  
end
#-------------------------------------------------------------------------------------------------------------------------------
#FOURIER FIT  WITH LEAST SQUARES METHOD
#Fourier series. Requires x and the series coefficients
#Returns an array with the truncation of the Fourier series evaluated at the x values
function fourier_expansion(x, c)
    # Verifica che x e c siano vettori
    if !(x isa AbstractVector)
        throw(ArgumentError("Il parametro x deve essere un AbstractVector, ma è di tipo $(typeof(x))"))
    end
    if !(c isa AbstractVector)
        throw(ArgumentError("Il parametro c deve essere un AbstractVector, ma è di tipo $(typeof(c))"))
    end
    # Verifica che la lunghezza di c sia minore o uguale alla lunghezza di x
    if length(c) > length(x)
        throw(ArgumentError("Il numero dei coefficienti deve essere minore o uguale alla lunghezza dei valori in x"))
    end
    # Verifica che la lunghezza di c sia dispari
    if length(c) % 2 == 0
        throw(ArgumentError("Il numero dei coefficienti deve essere dispari, ma è $(length(c)). \nLo sviluppo di Fourier ha sia coseno che seno, e il termine noto."))
    end
    sum = zeros(length(x))
    #j cicle just to evaluate the function in the x values
    for j in 1:1:length(x)
        sum[j] = c[1] #c[1] is the constant term
        for i in 2:2:length(c)
            sum[j] += c[i]*cos((i/2)*x[j]) + c[i+1]*sin((i/2)*x[j])
        end
    end
    return sum
end
#Goniometric functions fit. Requires x data, the y data, the number of coefficients you want from the fit
function fourier_fit(x, y, n_coeff)
    #Sezione di test iniziali
    # Verifica che x e y siano vettori
    if !(x isa AbstractVector)
        throw(ArgumentError("Il parametro x deve essere un AbstractVector, ma è di tipo $(typeof(x))"))
    end
    if !(y isa AbstractVector)
        throw(ArgumentError("Il parametro c deve essere un AbstractVector, ma è di tipo $(typeof(c))"))
    end
    # Verifica che la lunghezza di x e y siano uguali
    if length(x) != length(y)
        throw(ArgumentError("La lunghezza di x e y devono essere uguali, ma sono $(length(x)) e $(length(y))"))
    end
    # Verifica che n_coeff sia un numero intero positivo
    if !(n_coeff isa Integer) || n_coeff <= 0
        throw(ArgumentError("Il numero di coefficienti deve essere un intero positivo, ma è $(n_coeff)"))
    end
    # Verifica che n_coeff sia minore o uguale alla lunghezza di x
    if n_coeff > length(x)
        throw(ArgumentError("Il numero di coefficienti deve essere minore o uguale alla lunghezza di x, ma sono $(n_coeff) e $(length(x))"))
    end

    #Creating the design matrix (in italiano matrice dei regressori)
    A = hcat(ones(length(x), 1))
    for i in 2:1:n_coeff
        if i % 2 == 0
            A = hcat(A, cos.((i/2).*x))
        else
            A = hcat(A, sin.(((i-1)/2).*x))
        end
    end
    c = least_sq(A, y_meas)
    return c  
end
#-------------------------------------------------------------------------------------------------------------------------------
#CHANGE OF VARIABLES
#Change of variable for finite intervals
#Requires the values between [-1, 1] and the interval extrema
function ChangeOfVar(xn::AbstractVector{<:Real}, a::Number, b::Number)
    mn, mx = minimum(xn), maximum(xn) 
    if mx>1 || mn<-1
        throw(DomainError(:xn, "xn must contain values between [-1, 1]. In $(:xn) found min=$mn and max=$mx."))
    end
    if b < a
        change = b
        b = a
        a = change
        display("ChangeOfVar says: the interval extrema were swapped!")
    end
    xn_ab = [a + (b-a)*(x+1)/2 for x in xn]
    return xn_ab
end
#Change of variable for the entire real line
#Requires the values between [-1, 1]
function ChangeOfVar(x::Number)
    if x>=1 || x<=-1
        throw(DomainError(:x, "x must be a value between (-1, 1). $(:x) was $x."))
    end
    return 2x/(1-x^2) 
end
#Inverse:
function ChangeOfVarInvFinite(xn::AbstractVector{<:Real})
    a, b = minimum(xn), maximum(xn) 
    if a == b
        throw(DomainError(:xn, "xn must contain non constant values. In $(:xn) found inferior bound a=$a and upper bound b=$b."))
    end
    x = 2 .* (xn .- a)./(b-a) .- 1
    return x
end
function ChangeOfVarInvInfinite(x::Number)
    return (-1 + sqrt(1 + x^2)) / x
end
#-------------------------------------------------------------------------------------------------------------------------------
#SUPPORT COEFFICIENTS FOR LAGRANGIAN INTERPOLATION
#Has two methods:
# 1. wj(j::Int, xn::AbstractVector) for generic nodes
# 2. wj(j::Int, n::Int; pt_type="unif") for uniform or Chebyshev points

#where xn are the x values
function wj(j::Int, xn::AbstractVector)
    n = length(xn)
    prod=1
    for i in 1:1:n 
        if xn[i] != xn[j]
            prod *= xn[j] -xn[i]
        end
    end

    return 1/prod
end
#where n is the number of nodes, pt_type is the type of points (uniform or Chebyshev)
function wj(j::Int, n::Int, pt_type::Symbol)
    if !(pt_type in (:unif, :c1, :c2))
        throw(DomainError(pt_type, "Must be :unif, :c1, :c2. Inserted $pt_type"))
    end
    if pt_type == :unif
        return (-1)^(j+1)*factorial(big(n))/(factorial(big(j))*factorial(big(n-j))) 
    end
    if pt_type == :c1
        return (-1)^j*sin(pi*(2j-1)/2n)
    end
    if pt_type == :c2
        if j==1 || j==n
            return 0.5(-1)^j
        else 
            return (-1)^j
        end
    end
end
#-------------------------------------------------------------------------------------------------------------------------------
#LAGRANGIAN INTERPOLATION
#lag_fit calculates a function p(x), approximating the function f(x) in the nodes xn, using Lagrangian interpolation.
#Has three methods:
# 1. lag_fit(xn::AbstractVector, f::Function) for generic nodes
# 2. lag_fit(n::Int, a::Number, b::Number, f::Function, pt_type::Symbol) for uniform or Chebyshev points in a finite domain
# 3. lag_fit(n::Int, f::Function) for interpolation on the real line. It is used and explained in 4_4_4.ipynb

#xn is the x values, f is the function to be interpolated, n is the number of nodes, a and b are the extrema of the interval, pt_type is the type of points (uniform or Chebyshev)
function lag_fit(xn::AbstractVector, f::Function)
    yn = f.(xn)

    #nodes' number
    n = length(xn)
    weights = [wj(j, xn) for j in 1:1:n]
    
    function p(x)
        num = 0.0
        den = 0.0
        for j in 1:1:n
            if x != xn[j]
                num += weights[j]*yn[j]/(x-xn[j])
                den += weights[j]/(x-xn[j])
            else
                return yn[j]
            end
        end
        return num/den
    end

    return p
end
#where n is the number of nodes, a and b are the extrema of the interval, f is the function to be interpolated and pt_type is the type of points (uniform or Chebyshev)
function lag_fit(n::Int, a::Number, b::Number, f::Function, pt_type::Symbol)
    if a>b
        c=a
        a=b
        b=c 
        display("lag_fit says: the interval extrema were swapped.")       
    end
    #Case of uniform points 
    if pt_type == :unif
        xn = ChangeOfVar(range(-1, 1, length=n), a, b)
        weights = [wj(j, n, pt_type) for j in 1:1:n]
    #Case of first kind Chebyshev points (endpoints escluded)
    elseif pt_type == :c1
        #change of coordinates inside xn declaration
        xn = ChangeOfVar([cos((2j-1)pi/(2n)) for j in 1:1:n], a, b)
        weights = [wj(j, n, pt_type) for j in 1:1:n]
    #Case of second kind Chebyshev points (including endpoints)
    elseif pt_type == :c2
        #change of coordinates inside xn declaration
        xn = ChangeOfVar([cos((j-1)pi/(n-1)) for j in 1:1:n], a, b)
        weights = [wj(j, n, pt_type) for j in 1:1:n]
    else
        throw(DomainError(pt_type, "Non valid pt_type. Must be :unif or :c1 or :c2. Inserted $pt_type."))
    end

    yn = f.(xn) 
    
    p = function(x)
        num = 0.0
        den = 0.0
        for j in 1:1:n
            if x != xn[j]
                num += weights[j]*yn[j]/(x-xn[j])
                den += weights[j]/(x-xn[j])
            else
                return yn[j]
            end
        end
        return num/den
    end

    return p, xn, yn
end
#Where n is the nodes number, f is the function on the real line
function lag_fit(n::Int, f::Function)
    #Using first kind Chebyshev nodes to avoid possible infinities
    #change of coordinates inside xn declaration
    xn = [cos((2j-1)pi/(2n)) for j in 1:1:n]
    zn = ChangeOfVar(xn)
    weights = [wj(j, n, :c1) for j in 1:1:n]

    yn = f.(zn) 
    
    #Performing the fit exploiting the fact that xn is between [-1, 1]
    #P accepts value on the entire real line, but the function changes them between (-1, 1) and performs the fit
    p = function(z)
        x = ChangeOfVarInvInfinite(z) 
        num = 0.0
        den = 0.0
        for j in 1:1:n
            if x != xn[j]
                num += weights[j]*yn[j]/(x-xn[j])
                den += weights[j]/(x-xn[j])
            else
                return yn[j]
            end
        end
        return num/den
    end

    return p, xn, zn, yn
end
#-------------------------------------------------------------------------------------------------------------------------------
#ERROR ANALYSIS
#Error indication function
function ErrIndFunc(xn::AbstractArray)
    n = length(xn)
    f = function(x)
        prod=1.0
        for i in 1:1:n
            if x!=xn[i]
                prod *= (x - xn[i])
            else
                return 0.0
            end
        end
        return prod
    end
    return f
end
#infinite norm. Computes the infinite distance between two functions in a given x domain
function InftyNorm(approx::Function, f::Function, x::AbstractVector)
    return maximum(abs.(f.(x)-approx.(x)))    
end
#-------------------------------------------------------------------------------------------------------------------------------
#CARDINAL FUNCTIONS
#x is the variable of the cardinal function, xk is the k-th node, N is the nodes' number
function tau_k(x::Number, xk::Number, N::Int)
    return sin(N*pi*(x-xk) / 2) / (N * sin(pi*(x - xk) / 2))
end
#-------------------------------------------------------------------------------------------------------------------------------
#TRIGONOMETRIC INTERPOLATION
#This function interpolates between [-1, 1] the f function. Needs n, from which we can find N = 2n + 1 number of nodes.
function tri_fit(n::Int, f::Function)
    N = 2n+1
    xn = [2k/N for k in -n:1:n]
    yn = f.(xn)

    p = function(x)
        sum = 0.0
        for k in 1:1:N
            sum += yn[k]*tau_k(x, xn[k], N)
        end
        
        return sum
    end

    return p
end
#It's the same function as before but returns various p, each one for a different node's number, contained in n::AbstractVector
function tri_fit(n::Vector{Int}, f::Function)
    
    polynomials = []
    for i in 1:1:length(n)
        N = 2n[i]+1
        xn = [2k/N for k in -n[i]:1:n[i]]
        yn = f.(xn)

        p = function(x)
            sum = 0.0
            for k in 1:1:N
                sum += yn[k]*tau_k(x, xn[k], N)
            end
            
            return sum
        end
        push!(polynomials, p)
    end

    return polynomials
end
#-------------------------------------------------------------------------------------------------------------------------------
println("interpolation.jl loaded correctly")
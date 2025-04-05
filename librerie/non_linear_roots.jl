#File's index
#1. bis_met_rec(func, a, b) function: calculates recursevly function's zero with bisection method. INCOMPLETE
#2. bis_met(func, A, B) function: calculates a function's zero with bisection method.
#3. bis_met_steps(func, A, B) function: returns the function's zero and the steps required 
#4. newt_met_steps(f, x1) function: calculates the steps required to find a function's zero with Newton's method. 
#--------------------------------------------------------------------------------------------------------------------------
#Import section
using ForwardDiff
#--------------------------------------------------------------------------------------------------------------------------
#Calculates a function's zero with bisection method. bis_met_rec is a recursive function. INCOMPLETE
#The problem of the recursive method is in escaping the recursion: the condition needs the intial values of a, and b.
function bis_met_rec(func, a, b)
    if !(func isa Function)
        throw(ArgumentError("The first parameter should be a Function, instead is a $(typeof(func))"))
    end
    if !(a isa Number)
        throw(ArgumentError("The second parameter should be a Number, instead is a $(typeof(a))"))
    end
    if !(b isa Number)
        throw(ArgumentError("The third parameter should be a Number, instead is a $(typeof(b))"))
    end
    #making sure the interval (a, b) is well defined
    if (a==b)
        throw(ArgumentError("a and b must be different"))
    elseif (b < a)
        c = b
        b = a
        a = c
    end
    
    #m is the middle point of the interval, it is defined so because of rounding errors
    m = a + (b-a)/2.0
    func0 = func(m)
    
    #If the interval is small enough the function returns the middle point
    if (b-a) < 2^(-16)*max(abs(a), abs(b))
        return func0
    end

    if func(m)*func(a) > 0
        bis_met(func, m, b)
    else
        bis_met(func, a, m)
    end
end
#--------------------------------------------------------------------------------------------------------------------------
#Calculates a function's zero with bisection method. 
function bis_met(func, A, B)
    if !(func isa Function)
        throw(ArgumentError("The first parameter should be a Function, instead is a $(typeof(func))"))
    end
    if !(A isa Number)
        throw(ArgumentError("The second parameter should be a Number, instead is a $(typeof(A))"))
    end
    if !(B isa Number)
        throw(ArgumentError("The third parameter should be a Number, instead is a $(typeof(B))"))
    end
    #making sure the interval (A, B) is well defined
    if (A==B)
        throw(ArgumentError("a and b must be different"))
    elseif (B < A)
        c = B
        B = A
        A = c
    end

    iter_max = 1000     #maximum number of iterations
    I_length = B-A      #length of the initial interval

    #initial interval endpoints
    b=B
    a=A
    iter = 0

    while I_length > 10^(-16)*max(abs(a), abs(b))
        m = a + abs(b-a)/2.0       #middle point, defined so because of rounding errors
        if func(m)*func(a) > 0
            a = m
        elseif func(m)*func(a) < 0
            b = m
        else
            return m
        end  

        #updating data
        I_length = abs(b -a)
        iter += 1
        
        if iter == iter_max
            println("bis_met says: iterations' maximum number reached. Found zero could be inaccurate.")
            return m, step_vector
        end
    end
    return m
end
#--------------------------------------------------------------------------------------------------------------------------
#Returns the function's zero and the steps required for achieving the function's zero
function bis_met_steps(func, A, B)
    if !(func isa Function)
        throw(ArgumentError("The first parameter should be a Function, instead is a $(typeof(func))"))
    end
    if !(A isa Number)
        throw(ArgumentError("The second parameter should be a Number, instead is a $(typeof(A))"))
    end
    if !(B isa Number)
        throw(ArgumentError("The third parameter should be a Number, instead is a $(typeof(B))"))
    end
    #making sure the interval (A, B) is well defined
    if (A==B)
        throw(ArgumentError("a and b must be different"))
    elseif (B < A)
        c = B
        B = A
        A = c
    end

    iter_max = 1000     #maximum number of iterations
    I_length = B-A      #length of the initial interval

    #initial interval endpoints
    b=B
    a=A
    iter = 0
    step_vector = Float64[]        #will contain all the values of m, the successive approximations of the function's root

    while I_length > 10^(-16)*max(abs(a), abs(b))
        m = a + abs(b-a)/2.0       #middle point, defined so because of rounding errors
        push!(step_vector, m)
        if func(m)*func(a) > 0
            a = m
        elseif func(m)*func(a) < 0
            b = m
        else                        #if i find exactly the function's zero, I return the step_vector
            return m, step_vector
        end  

        #updating data
        I_length = abs(b -a)
        iter += 1
        
        if iter == iter_max
            println("bis_met_steps says: iterations' maximum number reached. Found zero could be inaccurate.")
            return m, step_vector
        end
    end
    
    return m, step_vector
end
#--------------------------------------------------------------------------------------------------------------------------
#Calculates the steps for finding a function's zero with Newton's method. Requires the function, the starting point, the root's multiplicity.
#The first step is x1
function newt_met_steps(f, x1, q)
    if !(f isa Function)
        throw(ArgumentError("The first parameter should be a Function, instead is a $(typeof(f))"))
    end
    if !(x1 isa Number)
        throw(ArgumentError("The second parameter should be a Number, instead is a $(typeof(x1))"))
    end
    if !(q isa Number)
        throw(ArgumentError("The third parameter should be a Number, instead is a $(typeof(q))"))
    end

    x = Float64[x1]
    df = x -> ForwardDiff.derivative(f, x)
    
    x_tol = 1.0e-15
    y_tol = 1.0e-15
    iter_max = 1000
    iter = 1

    while abs(f(x[iter])) >= y_tol
        #This if statement could not be integrated in while because it could access x[0] (undefined)
        if iter>=2
            if abs(x[iter] - x[iter-1]) <=  x_tol
                return x[iter], x
            end
        end
        if iter > iter_max
            println("newt_met_steps says: iterations' maximum number reached. Found zero could be inaccurate.")
            return x[iter], x
        end
        push!(x, x[iter] - q * f(x[iter])/df(x[iter]))
        iter += 1
    end

    return x[iter], x
end
#--------------------------------------------------------------------------------------------------------------------------
#Calculates the steps for finding a function's zero with the secant's method. Requires the function, the starting points, the root's multiplicity.
#The first step is x1
function secant_met_steps(f, x1, x2)
    if x1==x2
        throw(ArgumentError("The initial interval must be non degenerate."))
    end
    #I must be sure that x1 and x2 are in the right order
    if x1 > x2
        a=x2
        x2=x1
        x1=a
    end

    x = [x1, x2]
    x_tol = 1.0e-15
    y_tol = 1.0e-15
    iter = 1
    iter_max = 1000

    while abs(x[iter+1] - x[iter]) > x_tol && abs(f(x[iter+1])) > y_tol 
        if iter > iter_max
            println("secant_met_steps says: iterations' maximum number reached. Found zero could be inaccurate.")
            return x[iter], x
        end
        
        m = (f(x[iter+1]) - f(x[iter])) / (x[iter+1] - x[iter])
        push!(x, x[iter+1] - f(x[iter+1]/m))
        iter+=1
    end
    return x[iter], x
end
#---------------------------------------------------------------------------------------------------------------------------
println("non_linear_roots.jl loaded correctly")
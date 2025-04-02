#File's index:
# 1. pol_expansion: polynomial series
# 2. pol_fit: polynomial fit
# 3. fourier_expansion: Fourier series
# 4. fourier_fit: Fourier fit
#import section
using Markdown     #markdown visualization in display
using LaTeXStrings
include("c:\\ALL\\Stefano\\Bicocca\\3terzo_anno\\lab_comp\\lab_computazionale1\\librerie\\linear_systems.jl")
#-------------------------------------------------------------------------------------------------------------------------------
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
println("interpolation.jl loaded correctly")
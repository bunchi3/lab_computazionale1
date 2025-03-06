using Plots

function approx_exp(x, N)
    sum = 0
    for n = 0:N
        sum += x^n / factorial(n)
    end
    return sum
end

#Creo un array di valori di x
x = collect(0:0.1:1)
y = exp.(x)

N = [x for x in 1:4]

#Plot del grafico x y
fig = plot(x, y, label="exp(x)", title="Approximation of exp(x)", xlabel="x", ylabel="y", lw=2)
savefig(fig, "C:/ALL/Stefano/Bicocca/3terzo_anno/lab_comp/lab_computazionale1/esercizi/1lezione/exp_from_julia.png")

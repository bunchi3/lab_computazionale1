import numpy as np
import matplotlib.pyplot as plt
import math

#Funzione per lo sviluppo in serie dell'esponenziale
def approx_exp(x, N):
    sum = 0
    for n in range(0, N):
        sum += (x**n)/math.factorial(n)
        
    return sum

#Creo la funzione esatta
x = np.arange(0, 1, 0.1)
f = np.exp(x)

#Inizializzo il numero di confronti che farò
N = [1, 2, 3, 4]

#Genero la funzione approssimata con l'espansione in serie
#Ogni elemento di gn è una funzione approssimata con un diverso N
#è di fatto una lista di len(N) array
gn = []
for i in range(0, len(N), 1):
    gn.append(approx_exp(x, N[i]))

#Errore assoluto. E' anch'esso una lista di len(N) array
delta = np.abs(f-gn)

#Genero la funzione andamento.
#La immagazzino in una lista di len(N) array
delta_scale = [] 
for i in range(0, len(N), 1):    
    delta_scale.append(x**(i+1)/math.factorial(i+1))

fig, ax = plt.subplots(len(N))
for i in range(0, len(N), 1):
    #ax[i].title("N=")
    ax[i].plot(x, delta[i])
    ax[i].plot(x, delta_scale[i])
    ax[i].legend(["Absolute error", "Scaling function"])
plt.show()

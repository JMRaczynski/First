from itertools import product, cycle
from random import randrange, seed
from math import floor, ceil, log, factorial, e, pi
from time import time
import matplotlib.pyplot as plt


def main():
    seed(time())
    r = int(input("Wprowadź wartość parametru r "))  # zmienny parametr bazowy
    N = (2**r) ** r  # ilość cząsteczek
    #N = 72000
    R = 3 * r + 1  # maksymalna wartość położenia
    P = 2 * r + (1 - r % 2)  # maksymalna wartość pędu
    #P = 3
    t = 0  # czas
    times = []  # do wykresu
    values = [] # do wykresu
    print("P= ", P," R= ", R)
    deltat = 1/2 / P  # krok czasu
    states = tuple(product(range(-R, R), range(-R, R), range(-P, P), range(-P, P)))
    step = (2 * P) ** 2
    border = (2 * R) * step
    atoms = []
    ns = [0*i for i in range(len(states))]  # licznik ilości atomów w danym stanie
    atomsleft = N
    for i in cycle(range(0, border, step)):
        atoms.append(randrange(i, i + step))
        ns[atoms[len(atoms)-1]] = ns[atoms[len(atoms)-1]]+1
        atomsleft -= 1
        if atomsleft == 0:
            break
    print("Ilość stanów ", len(states))
    print("ilość cząsteczek", len(atoms))
    for i in range(N):
        print(atoms[i])

    def daje_ns_od_tj(j, atoms1, ns1):
        tj = j*deltat  # j to numer kroku
        print("Czas: ", round(tj, 2))
        times.append(tj)
        for i in range(len(atoms1)):
            xatoms = floor(states[atoms1[i]][0] + tj*states[atoms1[i]][2])  # zmiana połozenia X
            yatoms = floor(states[atoms1[i]][1] + tj*states[atoms1[i]][3])  # zmiana połozenia Y
            ns1[atoms1[i]] = ns1[atoms1[i]] - 1
            while xatoms >= R or xatoms < -R:
                if xatoms >= R:
                    xatoms = 2*R-xatoms - 1  # odbicie od ściany
                    atoms1[i] -= (2 * P + 1) * (2 * abs(states[atoms1[i]][2]))  # mnożenie wektora razy -1
                elif xatoms < -R:
                    xatoms = -xatoms-2*R - 1  # odbicie od ściany
                    atoms1[i] += (2 * P + 1) * (2 * abs(states[atoms1[i]][2]))  # mnożenie wektora razy -1
            while yatoms >= R or yatoms < -R:
                if yatoms >= R:
                    yatoms = 2*R-yatoms - 1 # odbicie od ściany
                    atoms1[i] -= 2*abs(states[atoms1[i]][3]) # mnożenie wektora razy -1
                elif yatoms < -R:
                    yatoms = -yatoms-2*R - 1 # odbicie od ściany
                    atoms1[i] += 2*abs(states[atoms1[i]][3]) # mnożenie wektora razy -1

            atoms1[i] += (xatoms - states[atoms1[i]][0])*border + (yatoms - states[atoms1[i]][1])*step  # zamiana miejsca atomu w tuple
            ns1[atoms1[i]] = ns1[atoms1[i]] + 1
            print("Położenie cząstki ", i, ":\t(", states[atoms1[i]][0], ",", states[atoms1[i]][1], ")", sep='')
            print("Pęd cząstki ", i, ":\t\t\t(", states[atoms1[i]][2], ",", states[atoms1[i]][3], ")", sep='')
        return ns1

    def entrop(N1, ns1):
        if N1 < 128:
            return log(prawdopodobienstwo(N1, ns1))
        else:
            ent = N1*log(N1/e)+log(2*pi*N1)/2
            il = 1
            for i in ns1:
                il *= factorial(i)
            return ent - log(il)

    def prawdopodobienstwo(N1, ns1):
        praw = factorial(N1)
        for i in ns1:
            praw /= factorial(i)
        return praw

    maxj = int(input("Wprowadź liczbę kroków czasu "))
    print('\n')
    praw = []
    for i in range(maxj):
        patoms = atoms[:]  # kopia tablicy atoms
        pns = ns[:]  # kopia tablicy
        ns = daje_ns_od_tj(i, atoms, ns) # przeprowadzenie symulacji dla kolejnego kroku czasu
        #for j in range(len(ns)):
            #if ns[j] != 0:
               # print("Liczba cząsteczek w stanie o indeksie ", j, ": ", ns[j], sep='')
        if N < 128:
            praw.append(prawdopodobienstwo(N, ns)) # obliczanie prawdopodobienstwa
        print('\n')
        entropia = entrop(N, ns) # obliczanie entropii
        values.append(entropia)
        #print(atoms[:])
        print(ns[:])
        atoms = patoms[:]
        ns = pns[:]
        #print(atoms[:])
        print(ns[:])
    print('\n\n')
    if N < 128:
        for i in range(len(times)):
            print("Czas: ", round(times[i], 2), " Wartość prawdopodobieństwa: ", '{:0.2e}'.format(praw[i]), sep='')
    plt.plot(times, values)
    plt.xlabel("Czas")
    plt.grid()
    plt.ylabel("Entropia układu")
    plt.title("Zależność entropii od czasu")
    axes = plt.gca()
    axes.set_xlim(0, floor(max(times) + 1))
    axes.set_ylim(ceil(min(values) - 1), floor(max(values) + 1))
    plt.show()
    del ns, atoms, states, values, times, patoms, pns


if __name__ == "__main__":
    main()

from itertools import product, cycle
from random import randrange
from math import floor, log, factorial
import matplotlib.pyplot as plt


def main():
    r = int(input("Wprowadź wartość parametru r"))  # zmienny parametr bazowy
    N = 2 ** r  # ilość cząsteczek
    R = 2 * r + 1  # maksymalna wartość położenia
    P = r + (1 - r % 2)  # maksymalna wartość pędu
    t = 0  # czas
    times = []  # do wykresu
    values = [] # do wykresu
    print("P= ", P," R= ", R)
    deltat = 1/2 / P  # krok czasu
    states = tuple(product(range(-R, R + 1), range(-R, R + 1), range(-P, P + 1), range(-P, P + 1)))
    step = (2 * P + 1) ** 2
    border = (2 * R + 1) * step
    atoms = []
    ns = [0*i for i in range(len(states))]  # licznik ilości atomów w danym stanie
    atomsleft = N
    for i in cycle(range(0, border, step)):
        atoms.append(randrange(i, i + step))
        ns[atoms[len(atoms)-1]] = ns[atoms[len(atoms)-1]]+1
        atomsleft -= 1
        if atomsleft == 0:
            break
    #print(len(states),states)
    for i in range(N):
        print(states[atoms[i]])

    def daje_ns_od_tj(j, atoms1, ns1):
        tj = j*deltat  # j to numer kroku
        print("Czas: ",tj)
        times.append(tj)
        # print(R, P, tj, len(atoms1))
        for i in range(len(atoms1)):
            xatoms = floor(states[atoms1[i]][0] + tj*states[atoms1[i]][2])  # zmiana połozenia X
            yatoms = floor(states[atoms1[i]][1] + tj*states[atoms1[i]][3])  # zmiana połozenia Y
            ns1[atoms1[i]] = ns1[atoms1[i]] - 1
            # print(xatoms, yatoms, i)
            while xatoms >= R or xatoms < -R:
                # print(xatoms)
                if xatoms >= R:
                    xatoms = 2*R-xatoms - 1  # odbicie od ściany
                    atoms1[i] -= (2 * P + 1) * (2 * abs(states[atoms1[i]][2]))  # mnożenie wektora razy -1
                elif xatoms < -R:
                    xatoms = -xatoms-2*R - 1  # odbicie od ściany
                    atoms1[i] += (2 * P + 1) * (2 * abs(states[atoms1[i]][2]))  # mnożenie wektora razy -1
            while yatoms >= R or yatoms < -R:
                # print(yatoms)
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

    def prawdopodobienstwo(N1, ns1):
        praw = factorial(N1)
        for i in ns1:
            praw /= factorial(i)
        return praw

    maxj = int(input("Wprowadź liczbę kroków czasu"))
    print('\n')
    for i in range(maxj):
        patoms = atoms[:]  # kopia tablicy atoms
        pns = ns[:]  # kopia tablicy
        ns = daje_ns_od_tj(i, atoms, ns)
        print('\n')
        entropia = log(prawdopodobienstwo(N, ns))
        values.append(log(prawdopodobienstwo(N, ns)))
        #print(entropia)
        atoms = patoms[:]
        ns = pns[:]
    #print(values)
    plt.plot(times, values)
    plt.xlabel("Czas")
    plt.grid()
    plt.ylabel("Entropia układu")
    plt.title("Zależność entropii od czasu")
    axes = plt.gca()
    axes.set_xlim(0, times[-1] + 5 - times[-1] % 5)
    axes.set_ylim(values[0] - 5 + values[0] % 5, values[-1] + 5 - values[-1] % 5)
    plt.show()


if __name__ == "__main__":
    main()

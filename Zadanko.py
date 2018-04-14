from itertools import product, cycle
from random import randrange
from math import floor, log


def main():
    r = int(input())  # zmienny parametr bazowy
    N = 2 ** r  # ilość cząsteczek
    R = 2 * r + 1  # maksymalna wartość położenia
    P = r + (1 - r % 2)  # maksymalna wartość pędu
    t = 0  # czas
    deltat = 1/2 / P  # krok czasu
    states = tuple(product(range(-R, R + 1), range(-R, R + 1), range(-P, P + 1), range(-P, P + 1)))
    step = (2 * P + 1) ** 2
    print(step)
    border = (2 * R + 1) * step
    print(border)
    atoms = []
    ns = [0*i for i in range(len(states))]
    atomsleft = N
    for i in cycle(range(0, border, step)):
        atoms.append(randrange(i, i + step))
        ns[atoms[len(atoms)-1]] = ns[atoms[len(atoms)-1]]+1
        atomsleft -= 1
        if atomsleft == 0:
            break
    for i in range(N):
        print(states[atoms[i]])
        # print(ns[atoms[i]])

    def daje_ns_od_tj(j, atoms1, ns1):
        tj = j*deltat
        # print(R, P, tj, len(atoms1))
        for i in range(len(atoms1)):
            xatoms = floor(states[atoms1[i]][0] + tj*states[atoms1[i]][2])
            yatoms = floor(states[atoms1[i]][1] + tj*states[atoms1[i]][3])
            ns1[atoms1[i]] = ns1[atoms1[i]] - 1
            # print(xatoms, yatoms, i)
            while xatoms >= R or xatoms < -R:
                # print(xatoms)
                if xatoms >= R:
                    xatoms = 2*R-xatoms - 1 # odbicie od ściany
                    atoms1[i] -= (2 * P + 1) * (2 * abs(states[atoms1[i]][2])) # mnożenie wektora razy -1
                elif xatoms < -R:
                    xatoms = -xatoms-2*R # odbicie od ściany
                    atoms1[i] += (2 * P + 1) * (2 * abs(states[atoms1[i]][2])) # mnożenie wektora razy -1
            while yatoms >= R or yatoms < -R:
                # print(yatoms)
                if yatoms >= R:
                    yatoms = 2*R-yatoms - 1 # odbicie od ściany
                    atoms1[i] -= 2*abs(states[atoms1[i]][3]) # mnożenie wektora razy -1
                elif yatoms < -R:
                    yatoms = -yatoms-2*R # odbicie od ściany
                    atoms1[i] += 2*abs(states[atoms1[i]][3]) # mnożenie wektora razy -1

            atoms1[i] += (xatoms - states[atoms1[i]][0])*border + (yatoms - states[atoms1[i]][1])*step # zamiana miejsca atomu w tuple
            ns1[atoms1[i]] = ns1[atoms1[i]] + 1
        return ns1

    def silnia(k):
        sil = 1
        i = 1
        while i <= k:
            sil *= i
            i = i + 1
        return sil

    def prawdopodobienstwo(N1, ns1):
        praw = silnia(N1)
        for i in ns1:
            praw /= silnia(i)
        return int(praw)

    #for i in range(7):
       # print(silnia(i))

    maxj = int(input())
    for i in range(maxj):
        patoms = atoms[:]
        pns = ns[:]
        ns = daje_ns_od_tj(i, atoms, ns)
        # !!! tutaj liczenie tych prawdopodobieńst !!! #
        entropia = log(prawdopodobienstwo(N, ns))
        print(entropia)
        # !!! tutaj liczenie tych prawdopodobieńst !!! #
        atoms = patoms[:]
        ns = pns[:]







if __name__ == "__main__":
    main()

from itertools import product, cycle
from random import randrange
from math import floor


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
    print(P)
    print(states[atoms[1]])
    atoms[1] += 5*border + 3*step
    print(states[atoms[1]])

    #czas
    j = int(input())
    tj = j*deltat
    for i in range(len(atoms)):
        xatoms = floor(states[atoms[i]][0] + tj*states[atoms[i]][2])
        yatoms = floor(states[atoms[i]][1] + tj*states[atoms[i]][3])
        ns[atoms[i]] = ns[atoms[i]] - 1
        while xatoms > R or xatoms < -R:
            print("a")
            if xatoms > R:
                xatoms = R - (xatoms - abs(R - states[atoms[i]][0])) # odbicie od ściany
                atoms[i] -= (2 * P + 1) * (2 * abs(states[atoms[i]][2])) # mnożenie wektora razy -1
            elif xatoms < R:
                xatoms = -R + xatoms - abs(-R - states[atoms[i]][0]) # odbicie od ściany
                atoms[i] += (2 * P + 1) * (2 * abs(states[atoms[i]][2])) # mnożenie wektora razy -1
        while yatoms > R or yatoms < -R:
            print("b")
            if yatoms > R:
                yatoms = R - (yatoms - abs(R - states[atoms[i]][1])) # odbicie od ściany
                atoms[i] -= 2*abs(states[atoms[i]][3]) # mnożenie wektora razy -1
            elif yatoms < -R:
                yatoms = -R + yatoms - abs(-R - states[atoms[i]][1]) # odbicie od ściany
                atoms[i] += 2*abs(states[atoms[i]][3]) # mnożenie wektora razy -1
        atoms[i] += (xatoms - states[atoms[i]][0])*border + (yatoms - states[atoms[i]][1])*step # zamiana miejsca atomu w tuple
        ns[atoms[i]] = ns[atoms[i]] + 1





if __name__ == "__main__":
    main()

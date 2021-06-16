import numpy as np
from scipy.spatial.transform import Rotation as R

filename = 'template.xyz'
n_points = 6  # number of rotattion points
A = 16  # fst atom
B = 10  # snd atom

A -= 1
B -= 1

table = {
    1: 'H',
    6: 'C',
    7: 'N',
    8: 'O',
}

with open(filename) as file:
    template = []
    bonds = False
    bonds_pair = []
    for line in file:
        if line == "Bonds\n":
            bonds = True
            continue
        if not bonds:
            splitted = line.split()
            charge = int(splitted[0])
            x = float(splitted[1])
            y = float(splitted[2])
            z = float(splitted[3])
            template.append([charge, np.array([x, y, z])])
        if bonds:
            splitted = line.split()
            fst = int(splitted[0]) - 1
            snd = int(splitted[1]) - 1
            bonds_pair.append([fst, snd])
    adj_list = [[] for i in range(len(template))]
    adj_list[1].append(2)
    for p in bonds_pair:
        fst = p[0]
        snd = p[1]
        if snd not in set(adj_list[fst]):
            adj_list[fst].append(snd)
            adj_list[snd].append(fst)

    # disconnect them
    adj_list[A].remove(B)
    adj_list[B].remove(A)

    comp_Al = []
    comp_Bl = []

    # A components
    used = set()
    queue = [A]
    while len(queue):
        v = queue.pop()
        used.add(v)
        comp_Al.append(v)
        for w in adj_list[v]:
            if w not in used:
                queue.append(w)

    # B components
    used.clear()
    queue.clear()
    queue = [B]
    while len(queue):
        v = queue.pop()
        used.add(v)
        comp_Bl.append(v)
        for w in adj_list[v]:
            if w not in used:
                queue.append(w)

    comp_A = set(comp_Al)
    comp_B = set(comp_Bl)

    # Rc = np.array([0., 0., 0.])
    # M = 0
    # for a in comp_A:
    #     M += 1
    #     Rc += template[a][1]
    # Rc /= M
    Rc = template[A][1].copy()
    # rotate
    n = template[A][1] - template[B][1]
    n /= np.sqrt(np.sum(n**2))
    for k in range(1, n_points + 1):
        phi = 2 * np.pi / n_points
        r = R.from_rotvec(phi * n).as_matrix()
        for a in comp_A:
            template[a][1] -= Rc
            template[a][1] = r @ template[a][1]
            template[a][1] += Rc
        file_k = str(k) + '.xyz'
        with open(file_k, 'w') as res:
            for u in range(len(template)):
                e = template[u]
                if u in comp_A or 1:
                    print(e[0], e[1][0], e[1][1], e[1][2], file=res)
                    #print(table[e[0]], e[1][0], e[1][1], e[1][2], file=res)
            #print(32, Rc[0], Rc[1], Rc[2], file=res)

#    print(adj_list)

import numpy as np
import os


def tri2node(tn, rho):

    nt = tn.shape[0]
    np1 = np.max(tn)

    count = np.zeros(np1, dtype=int)
    Quant = np.zeros(np1)

    for i in range(nt):
        for j in range(3):
            node_index = tn[i, j] - 1
            count[node_index] += 1
            Quant[node_index] += rho[i]

    nrho = Quant / count

    return nrho

def pre_process(case_dir):
    os.chdir(case_dir)

    files = [f for f in os.listdir(case_dir) if f.endswith('.d')]
    stime = [float(f[:-2]) for f in files]
    stime.sort()
    stime[0]=0

    with open('0.d', 'r') as fp01:
        lines = fp01.readlines()

    spt, lpt = int(lines[1].split()[1]), int(lines[1].split()[2])
    sct, lct = int(lines[2].split()[1]), int(lines[2].split()[2])

    sx, lx = None, None
    sy, ly = None, None
    srho, lrho = None, None
    sfdt, lfdt = None, None
    sT, lT = None, None
    sTr, lTr = None, None

    for i, line in enumerate(lines[3:], start=3):
        if line.split()[0]=='x':
            sx, lx = int(line.split()[1]), int(line.split()[2])
        elif line.split()[0]=='y':
            sy, ly = int(line.split()[1]), int(line.split()[2])
        elif line.split()[0]=='rho':
            srho, lrho = int(line.split()[1]), int(line.split()[2])
        elif line.split()[0]=='frac1':
            sfdt, lfdt = int(line.split()[1]), int(line.split()[2])
        elif line.split()[0]=='T':
            sT, lT = int(line.split()[1]), int(line.split()[2])
        elif line.split()[0]=='TR':
            sTr, lTr = int(line.split()[1]), int(line.split()[2])


    with open('0', 'rb') as fp02:
        a = np.fromfile(fp02, dtype=np.float32)
    a = a[1:]
    pt = a[spt-1:spt+lpt-1]


    nt = lct
    tn = np.zeros((nt, 3), dtype=int)
    tn[:, 0] = pt[3 * np.arange(nt)].astype(int)
    tn[:, 1] = pt[3 * np.arange(nt) +1].astype(int)
    tn[:, 2] = pt[3 * np.arange(nt) +2].astype(int)

    T_bar = []

    # Loop over the time steps
    for jj in range(101):
        i = np.where(np.array(stime) >= jj * 1e-11)[0]
        if len(i) == 0:
            continue
        i = i[0]
        istart = (stime[i] > 0) * 4 * nt

        with open(str(stime[i]), 'rb') as fp02:
            a = np.fromfile(fp02, dtype=np.float32)
            a = a[1:]

        x = a[sx - istart-1:sx - istart + lx-1] * 1e4
        y = a[sy - istart-1:sy - istart + ly-1] * 1e4
        rho = a[srho - istart-1:srho - istart + lrho-1]
        T = a[sT - istart-1:sT - istart + lT-1]

        rho_interpolation = tri2node(tn, rho)

        T_bar=get_T_bar(T_bar,x,y,T,rho_interpolation)
        # t_plot = np.arange(0, 101) * 1e-11
    return T_bar





def get_T_bar(T_bar,x,y,T,rho_interpolation):
    range_T = 70
    center_pos = np.where(np.sqrt(x ** 2 + y ** 2) < range_T)[0]
    T_bar.append(rho_interpolation[center_pos] @ (y[center_pos] * T[center_pos]) / (
                rho_interpolation[center_pos] @ y[center_pos]))
    return T_bar

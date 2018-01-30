#!/usr/bin/env python
# Copyright (C) 2017 Sur Herrera Paredes
import os
import pandas as pd
import re


if __name__ == "__main__":
    dirs = os.listdir('./')
    Res = []
    for d in dirs:
        if not os.path.isdir(d):
            continue

        files = os.listdir('./' + d + '/')
        # print(files)
        for f in files:
            if re.match('^mk_results.*\.txt$', f) is not None:
                print(f)
                print("===========")
                dat = pd.read_csv(d + '/' + f, sep="\t", header=0)
                # print(dat.head())

                ngenes = dat.shape[0]
                n_Dn = sum(dat.Dn > 0)
                n_Ds = sum(dat.Ds > 0)
                n_D = sum((dat.Dn + dat.Ds) > 0)
                n_P = sum((dat.Pn + dat.Ps) > 0)
                n_p = sum(dat.hg_p < 0.1)
                n_p2 = sum(dat.hg_p < 0.05)

                res = [d, f, str(ngenes), str(n_Dn),
                       str(n_Ds), str(n_D), str(n_P)
                       str(n_p), str(n_p2)]
                # print(d)
                # print(f)
                # print(ngenes)
                # print(n_dn)
                # print(n_p)
                # print(res)
                Res.append(res)

    print(Res)
    with open('summary.txt', 'w') as o:
        o.write("\t".join(['Strain', 'Comparison', 'ngenes', 'n_Dn',
                           'n_Ds', 'n_D', 'n_P',
                           'p0.1', 'p0.05']) + "\n")
        for r in Res:
            o.write("\t".join(r) + "\n")

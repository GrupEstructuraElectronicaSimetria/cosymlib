periodic_table_info = dict(H=[1, 1], He=[2, 2], Li=[3, 1], Be=[4, 2], B=[5, 3], C=[6, 4], N=[7, 5], O=[8, 6], F=[9, 7],
                           Ne=[10, 8], Na=[11, 1], Mg=[12, 2], Al=[13, 3], Si=[14, 4], P=[15, 5], S=[16, 6], Cl=[17, 7],
                           Ar=[18, 8], K=[19, 1], Ca=[20, 2], Sc=[21, 3], Ti=[22, 4], V=[23, 5], Cr=[24, 6], Mn=[25, 7],
                           Fe=[26, 8], Co=[27, 9], Ni=[28, 10], Cu=[29, 11], Zn=[30, 12], Ga=[31, 13], Ge=[32, 14],
                           As=[33, 15], Se=[34, 16], Br=[35, 17], Kr=[36, 18], Rb=[37, 1], Sr=[38, 2], Y=[39, 3],
                           Zr=[40, 4], Nb=[41, 5], Mo=[42, 6], Tc=[43, 7], Ru=[44, 8], Rh=[45, 9], Pd=[46, 10],
                           Ag=[47, 11], Cd=[48, 12], In=[49, 13], Sn=[50, 14], Sb=[51, 15], Te=[52, 16], I=[53, 17],
                           Xe=[54, 18], Cs=[55, 1], Ba=[56, 2], La=[57, 3], Ce=[58, 4], Pr=[59, 5], Nd=[60, 6],
                           Pm=[61, 7], Sm=[62, 8], Eu=[63, 9], Gd=[64, 10], Tb=[65, 11], Dy=[66, 12], Ho=[67, 13],
                           Er=[68, 14], Tm=[69, 15], Yb=[70, 16], Lu=[71, 17], Hf=[72, 18], Ta=[73, 19], W=[74, 20],
                           Re=[75, 21], Os=[76, 22], Ir=[77, 23], Pt=[78, 24], Au=[79, 25], Hg=[80, 26], Tl=[81, 27],
                           Pb=[82, 28], Bi=[83, 29], Po=[84, 30], At=[85, 31], Rn=[86, 32], Fr=[87, 1], Ra=[88, 2],
                           Ac=[89, 3], Th=[90, 4], Pa=[91, 5], U=[92, 6], Np=[93, 7], Pu=[94, 8], Am=[95, 9],
                           Cm=[96, 10], Bk=[97, 11], Cf=[98, 12], Es=[99, 13], Fm=[100, 14], Md=[101, 15], No=[102, 16],
                           Lr=[103, 17], Rf=[104, 18], Db=[105, 19], Sg=[106, 20], Bh=[107, 21], Hs=[108, 22],
                           Mt=[109, 23], Ds=[110, 24], Rg=[111, 25], Cn=[112, 26], Uut=[113, 27], Uuq=[114, 28],
                           Uup=[115, 29], Uuh=[116, 30], Uus=[117, 31], Uuo=[118, 32])


def atomic_number_to_element(z):
    for element, info in periodic_table_info.items():
        if int(z) == info[0]:
            return element


def element_to_atomic_number(symbol):
    for element, info in periodic_table_info.items():
        if symbol.capitalize() == element:
            return info[0]


def element_valence_electron(symbol):
    for element, info in periodic_table_info.items():
        if symbol.capitalize() == element:
            return info[1]
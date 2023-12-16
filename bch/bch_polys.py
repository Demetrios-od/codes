# create different generator polynomials of the same BCH code

import gfield, bch

D = 5
polys = [0x11d, 0x12b, 0x12d, 0x14d, 0x15f, 0x163, 0x165, 0x169, 0x171, 0x187, 0x18d, 0x1a9, 0x1c3, 0x1cf, 0x1e7, 0x1f5]

for poly in polys:
    gf = gfield.gfield(poly)
    bch_gpoly = bch.ebch.generator_polynomial(gf, D)
    print('poly = ', poly, ', wt = ', len(bch_gpoly['powers']), sep='')
    print('powers =', bch_gpoly['powers'])
    print('zeros =', bch_gpoly['zeros'])
    print()


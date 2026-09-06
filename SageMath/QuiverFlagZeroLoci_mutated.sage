# A sublist of F_list_transverse_ND_fails and its mutation-equivalent replacements

F_sublist_transverse_ND_fails = [F_25, F_102, F_104, F_116, F_126, F_188, F_195, F_211, F_212, \
                                 F_219, F_230, F_238, F_253, F_288, F_294, F_297, F_318, F_321, \
                                 F_327, F_341, F_342, F_361, F_385, F_413, F_421, F_422, F_423, \
                                 F_430, F_439, F_464, F_473, F_478, F_483, F_508, F_509, F_512, \
                                 F_517, F_519, F_521, F_524, F_529, F_549, F_558, F_559, F_562, \
                                 F_582, F_593, F_597, F_600, F_603, F_604, F_616, F_617, F_626, \
                                 F_637, F_662, F_665, F_679, F_699, F_717, F_723, F_734];

F_25_mutated = w*x*y + w*z + x + z + 1/y + z/(x*y) + 1/(w*x) + 1/(w*x*z) + 1/(w*x^2*y)
F_102_mutated = x*z*w + x*z + y*w + x + y + z + y*z^-1 + x^-1*z^-1 + y^-1*z^-1*w^-1 + x^-1*z^-2*w^-1
F_104_mutated = y*w + z*w + x + y + z + w + x^-1*y*z^-1*w + x^-1*w + x^-1 + x^-1*y^-1*z + x^-1*z^-1*w + x^-1*z^-1 + x^-1*y^-1 + x^-1*y^-1*w^-1
F_116_mutated = x*z + y*w + x + y + x*y^-1*z + z + w + x*y^-1 + z^-1*w^-1 + x^-1*z^-1 + y^-1*z^-1*w^-1
F_126_mutated = x*z + x*w + y*w + x^2*y^-1 + 2*x + y + z + w + x*y^-1 + x^-1*z^-1 + y^-1*z^-1*w^-1 + x^-1*z^-1*w^-1
F_188_mutated = x*y*z*w + x*y*w + x + y + z + w + x^-1 + 2*x^-1*y^-1 + x^-1*y^-1*w^-1 + x^-1*y^-1*z^-1 + x^-2*y^-2*z^-1 + x^-2*y^-2*z^-1*w^-1
F_195_mutated = x*z*w + x + y + z + w + y^-1*z*w + x^-1*y*z^-1 + 2*x^-1 + x^-1*y^-1*w + x^-1*w^-1 + x^-1*y^-1 + x^-2*z^-1 + x^-2*z^-1*w^-1
F_211_mutated = x*y*z*w + x*y + y*z + x*w + y + z + w + x^-1 + x^-1*w^-1 + x^-1*y^-1 + x^-1*y^-1*w^-1 + x^-2*y^-1*z^-1*w^-1 + x^-2*y^-2*z^-1*w^-1
F_212_mutated = x*y + z*w + y + z + w + w^-1 + z^-1*w^-1 + x^-1*y^-2*z + 2*x^-1*y^-2 + x^-1*y^-2*w^-1 + x^-1*y^-2*z^-1 + x^-1*y^-2*z^-1*w^-1
F_219_mutated = y*w + x + y + z + y*z^-1*w + w + 2*z^-1*w + y^-1 + x^-1 + y^-1*z^-1*w + x^-1*w^-1 + x^-1*z^-1 + x^-1*y^-1*w^-1 + x^-1*y^-1*z^-1
F_230_mutated = y*w + x + y + z + w + x^-1*y*z^-1 + y^-1 + x^-1 + x^-1*w^-1 + 2*x^-1*z^-1 + x^-1*y^-1*w^-1 + x^-1*y^-1*z^-1 + x^-2*z^-1*w^-1 + x^-2*y^-1*z^-1*w^-1
F_238_mutated = x*y*z*w + x*z*w + y*w + x + y + z + w + x^-1*y*z^-1 + y^-1 + x^-1 + x^-1*w^-1 + 2*x^-1*z^-1 + x^-1*y^-1*w^-1 + x^-1*y^-1*z^-1
F_253_mutated = w*z + w + w*z^2/y + w*z/y + x + y + y/z + z + z/x + 2/x + 1/(x*z) + z/(x*y) + 1/(x*y) + y/(w*x) + 2*y/(w*x*z) + y/(w*x*z^2) + 1/(w*x) + 1/(w*x*z)
F_288_mutated = x*y*z*w + x*z + x*w + y + z + w + z^-1 + 2*x^-1*z^-1 + x^-1*z^-1*w^-1 + x^-1*y^-1*z^-2 + x^-2*y^-1*z^-2 + x^-2*y^-1*z^-2*w^-1
F_294_mutated = x*y*z*w + x*y*w + x*z + y + z + w + z^-1*w + z^-1 + x^-1 + 2*x^-1*z^-1 + x^-1*z^-1*w^-1 + x^-2*y^-1*z^-1 + x^-2*y^-1*z^-2 + x^-2*y^-1*z^-2*w^-1
F_297_mutated = x*y*z*w + x*y + z*w + y + w + y^-1 + x^-1*w^-1 + x^-1*y^-1 + x^-1*z^-1*w^-1 + x^-1*y^-1*w^-1 + x^-1*y^-1*z^-1 + x^-1*y^-1*z^-1*w^-1
F_318_mutated = x*y*w + y*z*w + x*w + x + z + z^-1 + y^-1 + x^-1*w^-1 + y^-1*z^-1 + x^-1*z^-1*w^-1 + x^-1*y^-1*w^-1 + x^-1*y^-1*z^-1*w^-1
F_321_mutated = x*y + y*w + z*w + w + y*z^-1 + y^-1*z + z^-1 + y^-1 + x^-1*z^-1 + 2*x^-1*y^-1 + x^-1*y^-2*z + x^-1*y^-1*w^-1 + x^-1*y^-2*z*w^-1 + 2*x^-1*y^-1*z^-1 + 2*x^-1*y^-2 + x^-1*y^-2*w^-1 + x^-1*y^-2*z^-1
F_327_mutated = x*y^2*z*w + x*y*z*w + y*z*w + x*y + y*z + y*w + z*w + y + z + x^-1*y^-1 + x^-1*y^-1*w^-1 + x^-1*y^-1*z^-1 + x^-1*y^-1*z^-1*w^-1 + x^-1*y^-2*w^-1 + x^-1*y^-2*z^-1*w^-1
F_341_mutated = y*w + z*w + x + w + y^-1*z + w^-1 + y^-1*z*w^-1 + x^-1*y*z^-1 + y^-1 + x^-1 + y^-1*w^-1 + x^-1*w^-1 + x^-1*y^-1*z*w^-1 + x^-1*z^-1 + x^-1*z^-1*w^-1 + 2*x^-1*y^-1*w^-1 + x^-1*y^-1*z^-1*w^-1
F_342_mutated = x*y*w + y*z*w + x*y + y*w + z*w + y + w + x^-1*z + y^-1 + x^-1 + x^-1*y^-1*z + x^-1*y^-1 + x^-1*y^-1*w^-1 + x^-1*y^-1*z^-1*w^-1 + x^-1*y^-2*w^-1 + x^-1*y^-2*z^-1*w^-1
F_361_mutated = w^2*y/z^2 + w^2/(x*z) + 3*w*y/z + 2*w*y/z^2 + w/z + 2*w/x + 2*w/(x*z) + w/(x*y) + x + 3*y + 4*y/z + y/z^2 + 1/z + z/x + 2/x + 1/(x*z) + z/(x*y) + 1/(x*y) + y*z/w + 2*y/w + y/(w*z) + z/w + 1/w
F_385_mutated = x*y*z*w + x*y*w + x*z*w + y*z*w + x*w + y*w + z*w + w + x^-1 + z^-1*w^-1 + y^-1*w^-1 + y^-1*z^-1*w^-1 + x^-1*y^-1*w^-1
F_413_mutated = y*z*w + y*w + z*w + x + y + z + w + x^-1*z + z^-1 + y^-1 + 2*x^-1 + x^-1*y^-1*z + x^-1*w^-1 + x^-1*z^-1 + 2*x^-1*y^-1 + x^-1*z^-1*w^-1 + x^-1*y^-1*z^-1
F_421_mutated = x*y*z*w + x*y*z + x*y*w + x*y + y + z + w + w^-1 + y^-1 + z^-1*w^-1 + y^-1*w^-1 + y^-1*z^-1*w^-1 + x^-1*y^-1*w^-1 + x^-1*y^-1*z^-1*w^-1 + x^-1*y^-2*w^-1 + x^-1*y^-2*z^-1*w^-1
F_422_mutated = x*y^2 + x*y + y*w + x*y*w^-1 + y + w + x*w^-1 + z*w^-1 + w^-1 + y^-1*z*w^-1 + z^-1 + y^-1 + y^-1*w^-1 + x^-1*y^-1 + x^-1*y^-1*z^-1
F_423_mutated = y*z*w + x*y + y*w + x*y*w^-1 + x*y*z^-1 + x + w + x*y*z^-1*w^-1 + x*w^-1 + x*z^-1 + x*z^-1*w^-1 + w^-1 + z^-1 + y^-1 + z^-1*w^-1 + y^-1*w^-1 + y^-1*z^-1 + x^-1*y^-1 + y^-1*z^-1*w^-1
F_430_mutated = x*y^2 + x*y + z*w + y + w + w^-1 + z^-1 + y^-1 + z^-1*w^-1 + y^-1*w^-1 + x^-1*y^-1 + y^-1*z^-1*w^-1 + x^-1*y^-1*w^-1 + x^-1*y^-1*z^-1 + x^-1*y^-1*z^-1*w^-1 + x^-1*y^-2*w^-1 + x^-1*y^-2*z^-1*w^-1
F_439_mutated = w*y*z + w*y + w*z + w + w*y*z/x + w*y/x + w*z/x + w/x + x + x/y + y*z + y + z + y*z/x + y/x + z/x + 1/x + x/(w*y) + x/(w*y*z)
F_464_mutated = y*w + x + y + z + y*z^-1*w + w + z*w^-1 + z^-1*w + w^-1 + y^-1*z*w^-1 + x^-1 + y^-1*w^-1 + x^-1*z^-1 + x^-1*y^-1 + x^-1*y^-1*z^-1
F_473_mutated = x*y*z*w + y^2*w + 2*y*z*w + z^2*w + y*w + z*w + y*z^-1 + z^-1 + z^-1*w^-1 + y^-1*w^-1 + x^-1*z^-1 + x^-1*y^-1 + y^-1*z^-1*w^-1 + x^-1*y^-1*z^-1 + x^-1*y^-1*z^-1*w^-1 + x^-1*y^-2*w^-1 + x^-1*y^-2*z^-1*w^-1
F_478_mutated = w*(z + 1)*(x + y + 1)/y + x + z + 1/x + 1/z + (x + 1)*(z + 1)/(y*z) + (x + y*z + y + 1)/(w*z)
F_483_mutated = w^2*z/(x*y) + 2*w^2/(x*y) + w^2/(x*y*z) + w*z + w + w/y + w/(y*z) + w*z/x + 3*w/x + 2*w/(x*z) + w/(x*y) + w/(x*y*z) + x + y + 1/z + y/x + y/(x*z) + 1/x + 1/(x*z) + x/w + 1/w
F_508_mutated = x*y + x*y*z^-1*w + x + 2*y + z + y*z^-1*w + x^-1*y + z^-1*w + w^-1 + z^-1 + y^-1 + 2*x^-1 + y^-1*w^-1 + x^-1*w^-1 + x^-1*y^-1 + x^-1*y^-1*w^-1
F_509_mutated = w*(x+1)*(y+z+1) + (y+1)/z + (y+z+1)*(x*z+y+z)/(w*x*y*z)
F_512_mutated = x*y + y*w + x*y*w^-1 + x*y*z^-1 + y + z + w + x*w^-1 + z*w^-1 + 2*w^-1 + y^-1*z*w^-1 + z^-1 + y^-1 + 2*y^-1*w^-1 + x^-1*y^-1*z*w^-1 + x^-1*y^-1*w^-1 + x^-1*y^-2*z*w^-1 + x^-1*y^-2*w^-1
F_517_mutated = x*y^2*z*w + x*y*z*w + y*z*w + x*y + y*z + y*w + z*w + x + y + z + w^-1 + y^-1 + z^-1*w^-1 + y^-1*w^-1 + y^-1*z^-1*w^-1 + x^-1*y^-1*w^-1 + x^-1*y^-1*z^-1*w^-1 + x^-1*y^-2*w^-1 + x^-1*y^-2*z^-1*w^-1
F_519_mutated = x*y*w + x*y + z*w + x*y*z^-1 + y + z + w + x*y*z^-1*w^-1 + x*z^-1 + x*z^-1*w^-1 + w^-1 + z^-1 + y^-1 + 2*z^-1*w^-1 + y^-1*w^-1 + y^-1*z^-1 + 2*y^-1*z^-1*w^-1 + x^-1*y^-1*w^-1 + x^-1*y^-1*z^-1*w^-1 + x^-1*y^-2*w^-1 + x^-1*y^-2*z^-1*w^-1
F_521_mutated = x*y + y*w + z*w + y + z + w + y*z^-1 + y*z^-1*w^-1 + 2*w^-1 + y^-1*z*w^-1 + z^-1 + z^-1*w^-1 + y^-1*w^-1 + x^-1*z^-1 + x^-1*y^-1 + x^-1*z^-1*w^-1 + 2*x^-1*y^-1*w^-1 + x^-1*y^-2*z*w^-1 + x^-1*y^-1*z^-1 + x^-1*y^-1*z^-1*w^-1 + x^-1*y^-2*w^-1
F_524_mutated = w*(y + 1)*(y*z + 1)/(y*z) + z*(y + 1)*(x + 1) + x + 2*y + (y + 1)/(y*z) + 1/x + (y + 1)*(x*z + z + 1)/w
F_529_mutated = w*x + w/y + w/(x*z) + x*y + x*z + x + 2*z/y + 1/y + y/(x*z) + 2/x + 1/(x*z) + 1/(x*y) + y/w + 2*z/w + 1/w + z^2/(w*y) + z/(w*y) + y/(w*x) + z/(w*x) + 1/(w*x) + z/(w*x*y)
F_549_mutated = x*y^2*z*w + x*y^2*w + x*y*z*w + x*y*z + x*y*w + x*y + y*w + y + z + w + w^-1 + y^-1 + z^-1*w^-1 + y^-1*w^-1 + y^-1*z^-1*w^-1 + x^-1*y^-1*w^-1 + x^-1*y^-1*z^-1*w^-1 + x^-1*y^-2*w^-1 + x^-1*y^-2*z^-1*w^-1
F_558_mutated = x*y + y*w + z*w + y^2*z^-1 + 2*y + z + w + 2*y*z^-1 + w^-1 + y^-1*z*w^-1 + z^-1 + x^-1*y*z^-1 + 3*x^-1 + 3*x^-1*y^-1*z + x^-1*y^-2*z^2 + y^-1*w^-1 + 2*x^-1*z^-1 + 4*x^-1*y^-1 + 2*x^-1*y^-2*z + x^-1*y^-1*z^-1 + x^-1*y^-2
F_559_mutated = w*x*y + w*x + w*y + w + w/z + w/(x*z) + x*y*z + x*y + x*z + x + y*z + y + z + 1/z + 1/x + 2/(x*z) + 1/(x*y*z) + 1/(w*x) + 1/(w*x*z) + 1/(w*x*y) + 1/(w*x*y*z)
F_562_mutated = w^2/(x*z) + w*y + w + 2*w/x + 2*w/(x*z) + w/(x*y*z) + x + 2*y*z + 2*y + 2*z + 1/y + z/x + 2/x + 1/(x*z) + 1/(x*y) + 1/(x*y*z) + y*z^2/w + 2*y*z/w + y/w + z^2/w + 3*z/w + 2/w + z/(w*y) + 1/(w*y)
F_582_mutated = w*x*y*z + w*x*y + w*x + w + x*y*z + x*y + x + y*z + y + z + 1/z + 1/(x*y) + 1/(x*y*z) + y*z/w + y/w + 1/w + z/(w*x) + 2/(w*x) + 1/(w*x*z) + 1/(w*x*y) + 1/(w*x*y*z)
F_593_mutated = y^2*w + y*z + 2*y*w + x + 2*y + z + y*z^-1*w + w + z*w^-1 + z^-1*w + x^-1*y*z^-1*w + w^-1 + y^-1*z*w^-1 + z^-1 + y^-1 + x^-1 + x^-1*z^-1*w + y^-1*w^-1 + x^-1*z^-1
F_597_mutated = x*y^2*z + x*y*z*w + x*y^2 + x*y*z + x*y*w + x*y + y*z + z*w + y + z + y*w^-1 + y*z^-1*w^-1 + 2*w^-1 + z^-1 + y^-1 + 2*z^-1*w^-1 + y^-1*w^-1 + x^-1*w^-1 + y^-1*z^-1 + y^-1*z^-1*w^-1 + x^-1*z^-1*w^-1 + 2*x^-1*y^-1*w^-1 + 2*x^-1*y^-1*z^-1*w^-1 + x^-1*y^-2*w^-1 + x^-1*y^-2*z^-1*w^-1
F_600_mutated = w*x*y^2 + w*x*y + w*y*z + w*y + w*z + w + x*y + y + y/z + 2/z + 1/y + 1/(y*z) + 1/x + 1/(x*z) + 2/(x*y) + 2/(x*y*z) + 1/(x*y^2) + 1/(x*y^2*z) + 1/w + 1/(w*z) + 1/(w*y) + 1/(w*y*z) + 1/(w*x*y) + 1/(w*x*y*z) + 1/(w*x*y^2) + 1/(w*x*y^2*z)
F_603_mutated = w*x*z + w*x + w*z + 2*w + w/z + w*z/y + w/y + x*z + x + y + y/z + z + 2/z + z/y + 1/y + 1/(x*y) + 1/(x*y*z) + y/(w*z) + 1/w + 1/(w*z) + 1/(w*x*z) + 1/(w*x*y) + 1/(w*x*y*z)
F_604_mutated = w*x/z + w*y + w + w/z + w/(y*z) + x + x/z + y*z + 2*y + z + 1/z + 1/y + 1/(y*z) + 1/x + 1/(x*y) + y*z/w + y/w + z/w + 1/w + z/(w*x) + 1/(w*x) + z/(w*x*y) + 1/(w*x*y)
F_616_mutated = w*x + w*x/z + w*x/(y*z) + w + w/z + w/(y*z) + x*z + x + y*z + 2*y + y/z + 2*z + 2/z + 1/(y*z) + z/x + 1/x + y*z/w + y/w + 1/w + y*z/(w*x) + y/(w*x) + 1/(w*x);
F_617_mutated = w*x^2*y*z + w*x^2*z + 2*w*x*y*z + w*x*z + w*x + w*y*z + w + x*y*z + x*y + x*z + 2*x + x/y + y*z + 2*y + 1/z + 1/y + 1/(y*z) + y/x + 1/x + 1/(x*z) + 1/(x*y*z) + y/w + 1/w + y/(w*x) + 1/(w*x)
F_626_mutated = w*z + w + w*z^2/x + 3*w*z/x + 3*w/x + w/(x*z) + x + x/y + y + z + 1/z + 2*y*z/x + 4*y/x + 2*y/(x*z) + 2*z/x + 4/x + 2/(x*z) + y/w + y/(w*z) + 1/w + 1/(w*z) + y^2/(w*x) + y^2/(w*x*z) + 2*y/(w*x) + 2*y/(w*x*z) + 1/(w*x) + 1/(w*x*z)
F_637_mutated = x^2*w + x*y + x*z + 2*x*w + 2*x + y + z + x*z^-1*w + x*y^-1*w + w + y*w^-1 + z*w^-1 + y*z^-1 + y^-1*z + z^-1*w + y^-1*w + w^-1 + x^-1*y*w^-1 + x^-1*z*w^-1 + z^-1 + y^-1 + x^-1 + x^-1*w^-1
F_662_mutated = w*x^2 + w*x^2*z/y + 3*w*x + w*x/z + 2*w*x*z/y + w*x/y + 3*w + 2*w/z + w*z/y + w/y + w/x + w/(x*z) + x*y + 2*x + x*z/y + 2*y + y/z + 1/z + z/y + y/x + y/(x*z) + 2/x + 1/(x*z) + y/w + 1/w + y/(w*x) + y/(w*x*z) + 1/(w*x)
F_665_mutated = w*x*y + w*x*z + w*y + w*z + w + w*z/y + x*y + x*y/z + x + y + 2*y/z + 2/z + 2/y + y/(x*z) + 1/x + 2/(x*z) + 2/(x*y) + 1/(x*y*z) + 1/(x*y^2) + x*y/(w*z) + x/w + 2*y/(w*z) + 2/w + 2/(w*z) + 2/(w*y) + y/(w*x*z) + 1/(w*x) + 2/(w*x*z) + 2/(w*x*y) + 1/(w*x*y*z) + 1/(w*x*y^2)
F_679_mutated = w*x/(y*z) + w + w/z + w/y + w/(y*z) + w*y/x + w*z/x + w/x + x^2/(y*z) + 2*x + 2*x/z + x/y + 2*x/(y*z) + 3*y + y/z + 2*z + 2/z + 1/y + 1/(y*z) + y^2/x + y*z/x + 2*y/x + z/x + 1/x + x^2/w + 2*x*y/w + x*z/w + 2*x/w + y^2/w + y*z/w + 2*y/w + z/w + 1/w
F_699_mutated = w*x*y/z + w*x + w*x/z + 2*w*y/z + w + 2*w/z + w*y/(x*z) + w/(x*z) + 2*x*y + x*z + 3*x + x*z/y + x/y + 4*y + y/z + z + 1/z + 1/y + 2*y/x + y/(x*z) + 2/x + 1/(x*z) + x*y*z/w + 2*x*z/w + x*z/(w*y) + 2*y*z/w + y/w + 3*z/w + 1/w + z/(w*y) + y*z/(w*x) + y/(w*x) + z/(w*x) + 1/(w*x)
F_717_mutated = w + w/z + w/(y*z) + w/(x*z) + w/(x*y*z) + x*y*z + 2*x*y + x*z + 3*x + x/y + y*z + 4*y + z + 1/z + 2/y + 1/(y*z) + 2*y/x + 3/x + 1/(x*z) + 1/(x*y) + 1/(x*y*z) + x^2*y^2*z/w + 2*x^2*y*z/w + x^2*z/w + 3*x*y^2*z/w + 6*x*y*z/w + x*y/w + 3*x*z/w + x/w + 3*y^2*z/w + 6*y*z/w + 2*y/w + 3*z/w + 2/w + y^2*z/(w*x) + 2*y*z/(w*x) + y/(w*x) + z/(w*x) + 1/(w*x)
F_723_mutated = w^2/(y*z) + 2*w*x/y + w + 3*w/z + 2*w/y + 2*w/(y*z) + w/(x*z) + x^2*z/y + x*z + 4*x + 2*x*z/y + 2*x/y + y + 3*y/z + z + 4/z + z/y + 2/y + 1/(y*z) + 2*y/(x*z) + 2/x + 2/(x*z) + x^2*z/w + 2*x*y/w + 3*x*z/w + 2*x/w + y^2/(w*z) + 4*y/w + 2*y/(w*z) + 3*z/w + 4/w + 1/(w*z) + y^2/(w*x*z) + 2*y/(w*x) + 2*y/(w*x*z) + z/(w*x) + 2/(w*x) + 1/(w*x*z)
F_734_mutated = w*x*z^2 + 3*w*x*z + 3*w*x + w*x/z + 2*w*y*z + 4*w*y + 2*w*y/z + 2*w*z + 4*w + 2*w/z + w*y^2/x + w*y^2/(x*z) + 2*w*y/x + 2*w*y/(x*z) + w/x + w/(x*z) + x*z^2 + 4*x*z + 5*x + 2*x/z + x*z/y + 2*x/y + x/(y*z) + 2*y*z + 6*y + 4*y/z + 2*z + 6/z + 2/y + 2/(y*z) + y^2/x + 2*y^2/(x*z) + 2*y/x + 5*y/(x*z) + 1/x + 4/(x*z) + 1/(x*y*z) + x*z/w + 2*x/w + x/(w*z) + x*z/(w*y) + 2*x/(w*y) + x/(w*y*z) + 2*y/w + 2*y/(w*z) + 4/w + 4/(w*z) + 2/(w*y) + 2/(w*y*z) + y^2/(w*x*z) + 3*y/(w*x*z) + 3/(w*x*z) + 1/(w*x*y*z)

F_list_mutated = [F_25_mutated, F_102_mutated, F_104_mutated, F_116_mutated, F_126_mutated, F_188_mutated, F_195_mutated, F_211_mutated, F_212_mutated, \
                  F_219_mutated, F_230_mutated, F_238_mutated, F_253_mutated, F_288_mutated, F_294_mutated, F_297_mutated, F_318_mutated, F_321_mutated, \
                  F_327_mutated, F_341_mutated, F_342_mutated, F_361_mutated, F_385_mutated, F_413_mutated, F_421_mutated, F_422_mutated, F_423_mutated, \
                  F_430_mutated, F_439_mutated, F_464_mutated, F_473_mutated, F_478_mutated, F_483_mutated, F_508_mutated, F_509_mutated, F_512_mutated, \
                  F_517_mutated, F_519_mutated, F_521_mutated, F_524_mutated, F_529_mutated, F_549_mutated, F_558_mutated, F_559_mutated, F_562_mutated, \
                  F_582_mutated, F_593_mutated, F_597_mutated, F_600_mutated, F_603_mutated, F_604_mutated, F_616_mutated, F_617_mutated, F_626_mutated, \
                  F_637_mutated, F_662_mutated, F_665_mutated, F_679_mutated, F_699_mutated, F_717_mutated, F_723_mutated, F_734_mutated];


# Repaired only to nested non-degeneracy, not transverse

F_669_mutated = x^2*y^-1*z^2*w^2 + x^2*z*w + 2*x*z*w + 2*x*y^-1*z^2*w + 2*x*y^-1*z*w^2 + x*y + 2*x*z + 2*x*w + 2*x*y^-1*z*w + x + y + 2*z + y^-1*z^2 + 2*w + 2*y^-1*z*w + y^-1*w^2 + y*w^-1 + z*w^-1 + y*z^-1 + 2*y^-1*z + z^-1*w + 2*y^-1*w + w^-1 + z^-1 + y^-1 + x^-1*y*z^-1*w^-1 + x^-1*w^-1 + x^-1*z^-1 + x^-1*z^-1*w^-1
F_730_mutated = w/y + w/(y*z) + 3*w/x + 2*w/(x*z) + 3*w/(x*y) + 2*w/(x*y*z) + 3*w*y/x^2 + w*y/(x^2*z) + 6*w/x^2 + 2*w/(x^2*z) + 3*w/(x^2*y) + w/(x^2*y*z) + w*y^2/x^3 + 3*w*y/x^3 + 3*w/x^3 + w/(x^3*y) + x^2/y + 4*x + x/z + 2*x*z/y + 3*x/y + 5*y + y/z + 4*z + 1/z + 4*z/y + 2/y + 2*y^2/x + 2*y*z/x + 4*y/x + 4*z/x + 2/x + 2*z/(x*y) + x^3/w + x^3*z/(w*y) + 2*x^2*y/w + 3*x^2*z/w + x^2/w + x^2*z^2/(w*y) + x^2*z/(w*y) + x*y^2/w + 2*x*y*z/w + x*y/w + x*z^2/w + 2*x*z/w + x*z^2/(w*y)
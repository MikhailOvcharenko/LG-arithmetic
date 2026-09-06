# A sublist of F_list_transverse_ND_fails and its mutation-equivalent replacements

F_sublist_transverse_ND_fails = [F_730, F_729, F_725, F_717, F_712, F_709, F_699, F_679, F_710, F_678, \
                                 F_708, F_689, F_671, F_707, F_666, F_676, F_704, F_660, F_673, F_659, \
                                 F_701, F_653, F_670, F_34, F_159, F_164, F_308, F_524, F_176, F_538, \
                                 F_618, F_570, F_651, F_3, F_17, F_22, F_26, F_30, F_32, F_46, \
                                 F_58, F_59, F_60, F_63, F_68, F_111, F_115, F_124, F_128, F_210, \
                                 F_216, F_228, F_233, F_256, F_261, F_279, F_284, F_302, F_309, F_326, \
                                 F_341, F_343, F_350, F_368, F_387, F_19, F_28, F_33, F_37, F_71, \
                                 F_106, F_108, F_139, F_148, F_168, F_223, F_275, F_317, F_347, F_352, \
                                 F_388, F_431, F_444, F_547, F_551, F_35, F_39, F_140, F_145, F_169, \
                                 F_171, F_221, F_291, F_294, F_300, F_357, F_358, F_363, F_401, F_404, \
                                 F_405, F_465, F_466, F_467, F_502, F_38, F_161, F_178, F_361, F_402, \
                                 F_419, F_519, F_522, F_549, F_572, F_477, F_480, F_486, F_525, F_528, \
                                 F_533, F_545, F_634, F_663, F_672, F_688, F_495, F_531, F_537, F_577, \
                                 F_579, F_604, F_612, F_624, F_629, F_631, F_674, F_693, F_703, F_149, \
                                 F_245, F_443, F_327, F_311, F_158, F_378, F_400, F_448, F_461, F_342, \
                                 F_214, F_248, F_107, F_437, F_445, F_500, F_299, F_188, F_262, F_301];

F_730_mutated = w^2*z*(x*z + y^2)^4/(x^2*y^3) + w*(x*z + y^2)*(3*x^3*z^2 + x^2*y^3*z + 6*x^2*y^2*z + 3*x^2*y*z^2 + 3*x^2*z^2 + 3*x*y^4 + 6*x*y^3*z + 6*x*y^2*z + 3*y^5 + 3*y^4)/(x^2*y^3) + (3*x^4*z^2 + x^3*y^3*z + 6*x^3*y^2*z + 6*x^3*y*z^2 + 6*x^3*z^2 + x^2*y^4*z + 3*x^2*y^4 + 3*x^2*y^2*z^2 + 12*x^2*y^2*z + 6*x^2*y*z^2 + 3*x^2*z^2 + 6*x*y^5 + 6*x*y^4*z + 6*x*y^4 + 12*x*y^3*z + 6*x*y^2*z + 3*y^6 + 6*y^5 + 3*y^4)/(x^2*y^3*z) + (x*z + y^2)*(x + y + 1)^3/(w*x^2*y^3*z^2)
F_729_mutated = w*(y+z)*(x*y*z + x*z + 1)/z + (x^3*y^5*z^3 + 2*x^3*y^4*z^4 + 3*x^3*y^4*z^3 + x^3*y^3*z^5 + 6*x^3*y^3*z^4 + 3*x^3*y^3*z^3 + 3*x^3*y^2*z^5 + 6*x^3*y^2*z^4 + x^3*y^2*z^3 + 3*x^3*y*z^5 + 2*x^3*y*z^4 + x^3*z^5 + 3*x^2*y^4*z^2 + 6*x^2*y^3*z^3 + 6*x^2*y^3*z^2 + 3*x^2*y^2*z^4 + 3*x^2*y^2*z^2 + 6*x^2*y*z^4 + 6*x^2*y*z^3 + 3*x^2*z^4 + 3*x*y^3*z + 6*x*y^2*z^2 + 3*x*y^2*z + 3*x*y*z^3 + 6*x*y*z^2 + 3*x*z^3 + y^2 + 2*y*z + z^2)/(x^2*y^2*z^3) + 2*(y+z)*(x*y*z + x*z + 1)^2/(w*x^2*y^2*z^2) + (x*y*z + x*z + 1)/(w^2*x^2*y^2*z)
F_725_mutated = w*(z + 1)*((y + 1)*(x*y + z))/z + (z + 1)^2*((y + 1)*(x*y + z))^2/(x*y^2*z^2) - 8 + 2*(z + 1)*((y + 1)*(x*y + z))/(w*x*y^2*z) + 1/(w^2*x*y^2)
F_717_mutated = x^4*y^2*z^5*w^2 + 2*x^3*y^2*z^4*w^2 + 2*x^3*y*z^4*w^2 + x^2*y^2*z^3*w^2 + 4*x^2*y*z^3*w^2 + 2*x^2*y*z^3*w + x^2*z^3*w^2 + 2*x^2*y*z^2*w + 2*x*y*z^2*w^2 + x^2*y*z^2 + 2*x*y*z^2*w + 2*x*z^2*w^2 + 2*x*y*z*w + 2*x*z^2*w + 2*x*y*z + 2*x*z*w + z*w^2 + 2*x*z + 2*z*w + y + z + 2*w + 2*w^-1 + z^-1 + y^-1 + 2*z^-1*w^-1 + 2*x^-1*z^-1 + 2*x^-1*z^-1*w^-1 + 2*x^-1*y^-1*z^-1 + 2*x^-1*z^-2*w^-1 + 2*x^-1*y^-1*z^-1*w^-1 + 2*x^-1*y^-1*z^-2*w^-1 + x^-2*y^-1*z^-2 + 2*x^-2*y^-1*z^-2*w^-1 + x^-2*y^-1*z^-2*w^-2 + 2*x^-2*y^-1*z^-3*w^-1 + 2*x^-2*y^-1*z^-3*w^-2 + x^-2*y^-1*z^-4*w^-2
F_712_mutated = y*z^3 + y*z^2*w + 3*y*z^2 + x^-1*y*z^3 + y*z*w + x^-1*y*z^2*w + 3*y*z + 3*x^-1*y*z^2 + x^-1*y*z*w + x*z*w^-1 + x + y + 3*z + 3*x^-1*y*z + w + 2*x*w^-1 + z*w^-1 + x^-1*y + 3*x^-1*z + x^-1*w + x*z^-1*w^-1 + 2*w^-1 + 3*z^-1 + 6*x^-1 + 2*x*y^-1*z^-1*w^-1 + z^-1*w^-1 + 3*y^-1*z^-1 + 3*x^-1*z^-1 + 2*x*y^-1*z^-2*w^-1 + 2*y^-1*z^-1*w^-1 + 3*y^-1*z^-2 + 3*x^-1*y^-1*z^-1 + 2*y^-1*z^-2*w^-1 + 3*x^-1*y^-1*z^-2 + x*y^-2*z^-3*w^-1 + y^-2*z^-3 + y^-2*z^-3*w^-1 + x^-1*y^-2*z^-3
F_709_mutated = x*y*z^2*w + y*z^3 + x*y*z*w + y*z^2*w + 3*y*z^2 + x^-1*y*z^3 + y*z*w + 3*y*z + 3*x^-1*y*z^2 + x*w + x + y + 3*z + 3*x^-1*y*z + w + z*w^-1 + x^-1*y + 3*x^-1*z + 2*w^-1 + 3*z^-1 + 6*x^-1 + z^-1*w^-1 + 3*y^-1*z^-1 + 3*x^-1*z^-1 + 2*y^-1*z^-1*w^-1 + 3*y^-1*z^-2 + 3*x^-1*y^-1*z^-1 + 2*y^-1*z^-2*w^-1 + 3*x^-1*y^-1*z^-2 + y^-2*z^-3 + y^-2*z^-3*w^-1 + x^-1*y^-2*z^-3
F_699_mutated = x - 4 + (x + 1)*(w*(z + 1)*(y*z + 1)/x + (x + w*(z + 1)*(y*z + 1))^2/(w^2*x*y*z^2))
F_679_mutated = w*(z + 1)*(x*y + x*z + y) + x - 4 + (x + w*(z + 1)*(x*y + x*z + y))^2/(w^2*x^2*y*z^2)
F_710_mutated = w*(z + 1)*(x*y^2 + x*y + 1) + x*(y + 1)^3 + 3*y + z + 1/z + 3/y + 3/(x*y) + 3/(x*y^2) + 1/(x^2*y^3) + (z + 1)*(x*y^2 + x*y + 1)^2/(w*x^2*y^3*z)
F_678_mutated = y*z^2*w + y*z^2 + y*z*w + x^-1*y*z^2*w + 2*y*z + x^-1*y*z^2 + z*w + x^-1*y*z*w + x + y + 2*z + 2*x^-1*y*z + w + x^-1*z*w + x*w^-1 + x^-1*y + 2*x^-1*z + x^-1*w + x*z^-1*w^-1 + w^-1 + 2*z^-1 + y^-1 + 4*x^-1 + x*y^-1*z^-1*w^-1 + z^-1*w^-1 + 2*y^-1*z^-1 + 2*x^-1*z^-1 + x^-1*y^-1 + x*y^-1*z^-2*w^-1 + y^-1*z^-1*w^-1 + y^-1*z^-2 + 2*x^-1*y^-1*z^-1 + y^-1*z^-2*w^-1 + x^-1*y^-1*z^-2
F_708_mutated = w*(x^2*y*z+x*y*z+x*z+y)^2 + x^3*z + 3*x^2*z*(y+1)/y + 3*x*(y^2*z+y^2+2*y*z+z)/y^2 + (y^4+y^3*z+3*y^2*z+6*y^2+3*y*z+z)/y^3 + 3*(y^2*z+y^2+2*y*z+z)/(x*y^2*z) + 3*(y+1)/(x^2*y*z) + 1/(x^3*z^2) + (x^2*y*z+x*y*z+x*z+y)/(w*x^3*y^3*z^2)
F_689_mutated = w*(x*y+z)*(x^2*y^2+x^2*y*z+x^2*y+2*x*y*z+z^2)/(x^2*y) + 2*x*y/z + x + x/z + 2*y + 2*y/z + z + 1/z + 4*z/x + 4/x + 2*z/(x*y) + 2*z^2/(x^2*y) + 2*z/(x^2*y) + (x*y+z)*(x+z+1)^2/(w*x^2*y*z^2)
F_671_mutated = x*y*z^2*w + x*y*z*w + y*z^2*w + y*z^2 + x*z*w + y*z*w + 2*y*z + x^-1*y*z^2 + x*w + z*w + x + y + 2*z + 2*x^-1*y*z + w + x^-1*y + 2*x^-1*z + w^-1 + 2*z^-1 + y^-1 + 4*x^-1 + z^-1*w^-1 + 2*y^-1*z^-1 + 2*x^-1*z^-1 + x^-1*y^-1 + y^-1*z^-1*w^-1 + y^-1*z^-2 + 2*x^-1*y^-1*z^-1 + y^-1*z^-2*w^-1 + x^-1*y^-1*z^-2
F_707_mutated = w*x^2*y + w*x^2 + w*x*y*z + 3*w*x*y + 3*w*x + w*y*z + 3*w*y + 3*w + w*y/x + w/x + 3*x + x/z + 3*x/y + y + z + 2/z + 6/y + 3/x + 1/(x*z) + 3/(x*y) + 3/(w*y) + 2/(w*y*z) + 3/(w*y^2) + 3/(w*x*y) + 2/(w*x*y*z) + 3/(w*x*y^2) + 1/(w^2*x*y^2) + 1/(w^2*x*y^2*z) + 1/(w^2*x*y^3)
F_666_mutated = w*(x+z)*(y+1)*(z+1)/x + (z+1)*(x+y)/z + z + 2/z + (z+1)/(y*z) + (z+1)*(y+1)^2/(x*y) + (x+z)*(y+1)*(z+1)/(w*y*z^2)
F_676_mutated = x*z*w + x*y^-1*z^2*w + z^2*w + y^-1*z^3*w + x*w + z*w + x + 2*x*y^-1*z + 2*z + x*y^-2*z^2 + 4*y^-1*z^2 + x^-1*z^2 + 2*y^-2*z^3 + 2*x^-1*y^-1*z^3 + x^-1*y^-2*z^4 + 2*x*z^-1 + y*z^-1 + 2*x*y^-1 + 4*y^-1*z + 2*x^-1*z + 2*x^-1*y^-1*z^2 + x*z^-2 + 2*z^-1 + x^-1 + 2*z^-1*w^-1 + 2*y^-1*w^-1 + 2*x^-1*w^-1 + 2*x^-1*y^-1*z*w^-1 + 2*z^-2*w^-1 + 2*x^-1*z^-1*w^-1 + x^-1*z^-2*w^-2
F_704_mutated = y*z^3 + y*z^2*w + 3*y*z^2 + y*z*w + x^-1*y*z^2*w + 3*y*z + x^-1*y*z*w + x + y + 3*z + w + z*w^-1 + x^-1*w + 2*w^-1 + 3*z^-1 + x^-1 + z^-1*w^-1 + 3*y^-1*z^-1 + 2*y^-1*z^-1*w^-1 + 3*y^-1*z^-2 + 2*y^-1*z^-2*w^-1 + y^-2*z^-3 + y^-2*z^-3*w^-1
F_660_mutated = w^2/(x*z) + w^2/(x*y*z^2) + w + w/z + 2*w/(y*z) + w*y/x + w/x + 2*w/(x*z) + w/(x*y*z) + x + x/y + y + z + 2/y + 1/(y*z) + y*z/x + y/x + 2/x + x*z/w + x/w + x*z/(w*y) + x/(w*y) + y*z/w + 2*z/w + 1/w + y*z/(w*x)
F_673_mutated = w*y*z^2 + w*y*z + w*z + w + w/x + x*y*z^2 + x*y*z + x*z + x + y*z^2 + 2*y*z + y + 2*z + 2/z + 1/y + 2/(y*z) + 1/(y*z^2) + 1/x + 1/(x*z) + 1/(x*y*z) + 1/(x*y*z^2) + x/w + 1/w + 1/(w*z) + 1/(w*y*z) + 1/(w*y*z^2)
F_659_mutated = w^3*x/(y*z) + w^3/(y*z) + w^2*x/y + w^2*x/(y*z) + w^2/y + w^2/(y*z) + w*x/z + w*x/y + w + 2*w/z + w/y + w/(x*z) + x + z + 1/x + y/w + 2*z/w + 1/w + y/(w*x) + y*z/w^2 + y/w^2 + z/w^2 + y*z/(w^2*x) + y/(w^2*x) + y*z/w^3 + y*z/(w^3*x)
F_701_mutated = w^2*y*z^3/x + w^2*z/x + w^2*y*z^3/x^2 + w^2*z/x^2 + w*y*z^3 + 3*w*z + 3*w/(y*z) + w/(y^2*z^3) + w*y*z^2/x + w*y*z^2/x^2 + x + 3*y*z^2 + 3/(y*z^2) + x/w + 2*x/(w*y*z^2) + x/(w*y^2*z^4) + 3*y*z/w + 3/(w*z) + 2*x/(w^2*z) + 2*x/(w^2*y*z^3) + y/w^2 + x/(w^3*z^2)
F_653_mutated = w*z^2/(x*y^2) + w*y*z/x^2 + w*z^3/(x^2*y^2) + x^2/(y^2*z) + x^2/y^3 + 2*x*y/z^2 + 2*x/z + 2*x/y^2 + 2*x*z/y^3 + y^4/z^3 + y^3/z^2 + 2*y/z + z/y + z/y^2 + z^2/y^3 + y^2/x + z^2/(x*y) + x^3/(w*y*z^2) + 2*x^2*y^2/(w*z^3) + 2*x^2/(w*y*z) + x*y^5/(w*z^4) + 2*x*y^2/(w*z^2) + x/(w*y)
F_670_mutated = w*x*z + w*x + w*y*z + w + x*z + x + x*z/y + 3*x/y + 3*x/(y*z) + x/(y*z^2) + 2*z + 2/z + y*z/x + y/x + 2*x*z/(w*y) + 6*x/(w*y) + 6*x/(w*y*z) + 2*x/(w*y*z^2) + 2*z/w + 4/w + 2/(w*z) + 1/(w*y) + 2/(w*y*z) + 1/(w*y*z^2) + 1/(w*x) + 1/(w*x*z) + x*z/(w^2*y) + 3*x/(w^2*y) + 3*x/(w^2*y*z) + x/(w^2*y*z^2) + 1/(w^2*y) + 2/(w^2*y*z) + 1/(w^2*y*z^2)
F_34_mutated =  x*y^2*z*w + x*y + y + z + w + z^-1*w^-1 + 2*x^-1*y^-1*z^-1 + x^-1*y^-1*z^-1*w^-1 + x^-1*y^-2*z^-2*w^-1 + x^-2*y^-3*z^-2 + x^-2*y^-3*z^-2*w^-1
F_159_mutated = x*y*z*w + x + y + z + w + z^-1 + z^-1*w^-1 + 2*x^-1*z^-1 + x^-1*z^-2*w^-1 + x^-1*y^-1*z^-2*w^-1 + x^-2*y^-1*z^-2 + x^-2*y^-1*z^-3*w^-1
F_164_mutated = w + w/(x*y*z) + y + y*z + x*y*z + x*y*z/w + 1/(w*y) + 1/(w*y*z) + 1/(x*y^2*z) + 1/(x*y^2*z^2)
F_308_mutated = w + w/x + w/(x*y*z) + x*y*z + x + y*z + y + 1/(y^2*z) + 1/x + 1/(x*y^2*z) + x*y*z/w + x/(w*y) + 1/(w*y)
F_524_mutated = x*y*z + x + y + z + w*z + w + 2*w/y + w/(y^2*z) + w/(x*y) + w/(x*y^2*z) + 2/y + 1/(y^2*z) + 1/(x*y) + 1/(x*y^2*z) + y/w + 1/w
F_176_mutated = w*x*y + w*z + w + w/(x*z) + x*y + x + 2*z + 2/(x*z) + 1/(x*y*z) + z/w + 1/w + 1/(w*x*z) + 1/(w*x*y*z)
F_538_mutated = w + w/(y*z) + w/(x*z) + w/(x*y*z) + x*z + x + y*z + y + z + 1/(y*z) + 1/(x*z) + 1/(x*y*z) + x*z/w + y*z/w + z/w + 1/w
F_618_mutated = w^3/(x^2*y*z) + w^2/(x*y) + w^2/(x*y*z) + w^2/(x^2*z) + w^2/(x^2*y*z) + w + w/y + 2*w/x + w/(x*y) + w/(x*y*z) + x + y + z + 1/y + x*y*z/w + 2*x/w + y/w + 1/w + x*y*z/w^2 + x/w^2
F_570_mutated = w^2*z/(x^2*y^2) + 2*w^2*z/(x^3*y) + w^2*z/x^4 + 2*w*z/(x*y) + 2*w/(x*y) + 2*w*z/(x*y^2) + 2*w*z/x^2 + 2*w/x^2 + 8*w*z/(x^2*y) + 8*w*z/x^3 + 2*w*y*z/x^4 + x + z + 1/z + 2*z/y + 2/y + z/y^2 + 6*z/x + 6/x + 8*z/(x*y) + 2*y*z/x^2 + 2*y/x^2 + 15*z/x^2 + 8*y*z/x^3 + y^2*z/x^4 + 2*z/w + 2/w + 2*z/(w*y) + 2*y*z/(w*x) + 2*y/(w*x) + 8*z/(w*x) + 8*y*z/(w*x^2) + 2*y^2*z/(w*x^3) + z/w^2 + 2*y*z/(w^2*x) + y^2*z/(w^2*x^2)
F_651_mutated = w + w*z/x + 2*w/x + w/(x*z) + 2*w/(x^2*y*z) + 2*w/(x^2*y*z^2) + w/(x^3*y^2*z^3) + x^2*y*z^2 + x^2*y*z + x + z/x + 2/x + 1/(x*z) + 2/(x^2*y*z) + 2/(x^2*y*z^2) + 1/(x^3*y^2*z^3) + x^2*y*z^2/w + x^2*y*z/w + x/w + 1/w
F_3_mutated = w*x*y*z + w + x + z + 1/(w*x*z) + 1/(w*x^2*y*z^2)
F_17_mutated = w*y + w + x + y + z + y/(w*z) + 1/(w*z) + 1/(w*x*z) + 1/(w*x*y*z)
F_22_mutated = w*x + w + x + y + z + 1/(x*y) + x/(w*z) + 1/(w*z)
F_26_mutated = w*x*y*z + w*x + w*y*z + x + z + 1/(w*x*z) + 1/(w*x*y*z^2) + 1/(w*x^2*z)
F_30_mutated = w*x + w*y + w + x + y + z + 1/(w*y*z) + 1/(w*x*z) + 1/(w*x*y*z)
F_32_mutated = w*y*z + w*y + w*z + w + w*y*z/x + w*z/x + x + y + 1/(w^2*z) + 1/(w^2*y*z)
F_46_mutated = w*x*y + w*z^2 + x + z + 1/z + 1/(w*x) + z^2/(w*x^2*y)
F_58_mutated = w*z + w + x + y + z + 1/z + z/(w*x*y) + 1/(w*x*y)
F_59_mutated = w*y*z + w*y + w*z + w + x + z + 1/z + z/(w^2*x*y) + 1/(w^2*x*y)
F_60_mutated = w*x*y + w*z + x + z + 1/z + 1/(w*x) + z/(w*x^2*y)
F_63_mutated = w*x*y*z + w*x*y + w*z^2 + w*z + x + z + 1/z + 1/(w*x) + z/(w*x^2*y)
F_68_mutated = w*y*z + w + x + y*z + y + z + 1/z + z/(w*x) + 1/(w*x*y)
F_111_mutated = w*x*y + w + x + z + 1/z + 1/(w*x) + 1/(w*x^2*y)
F_115_mutated = w*y + w + x + y + z + 1/z + 1/(y*z) + y/(w*x) + 1/(w*x)
F_124_mutated = w*y + w*z + w*y/x + w*z/x + x + z + 1/z + 1/w + z/(w*y)
F_128_mutated = w*x*y + w*z + w + x + z + 1/z + 1/(w*x) + z/(w*x^2*y) + 1/(w*x^2*y)
F_210_mutated = w*x*y + w*z + w + x + y + z + 1/z + z/(w*x*y) + 1/(w*x*y)
F_216_mutated = w*z + w + x + y + z + 1/z + 1/y + z/(w*x) + 1/(w*x)
F_228_mutated = w*y*z^2 + w*y*z + w*z + w + x + 1/x + 1/(w*z) + 1/(w^2*y*z^2)
F_233_mutated = w*x*y*z^2 + w*x*y*z + w*x*z + w*x + w*z^2 + w*z + x + 1/(w*x*z) + 1/(w*x*y*z^2) + 1/(w*x^2*y*z)
F_256_mutated = w*z + w + x + y + z + 1/z + 1/y + 1/(w*x) + 1/(w*x*z)
F_261_mutated = w*y + w + x + y + z + 1/z + 1/y + 1/(y*z) + y/(w*x) + 1/(w*x)
F_279_mutated = w*y + w + x + y + z + y/(x*z) + 2/(x*z) + 1/(x*y*z) + 1/w + 1/(w*y)
F_284_mutated = w*y + w + x + y + y/z + z + 1/z + 1/y + 1/(w*x) + 1/(w*x*y)
F_302_mutated = w*x*y*z + w*x + w*y*z^2 + w*y*z + w*z + w + w*y*z^2/x + w*z/x + x + x/z + z + 1/(w*z) + 1/(w*y*z^2)
F_309_mutated = w*y + w*y/z + w + w/z + w*y/x + w/x + x + y + z + 1/w + 1/(w*y)
F_326_mutated = w*y*z^2 + w*y*z + w*z + w + x + 1/x + 1/(w*z) + 1/(w*y*z^2) + 1/(w^2*x*y*z^2)
F_341_mutated = w*y*z^2 + w*y*z + w*z + w + x + z + 1/x + 1/(w*z) + 1/(w*y*z^2)
F_343_mutated = w*x*z^2 + w*x*z + w*y*z^2 + w*y*z + w*z + w + x + 1/x + x/(w*y*z) + 1/(w*z) + 1/(w*y*z^2)
F_350_mutated = w*x*y*z^2 + w*x*z + w*y*z^2 + w*y*z + w*z + w + x*z + x + 1/x + 1/(w*z) + 1/(w*y*z^2)
F_368_mutated = w*y + w + w*y/x + w/x + x + y + z + 1/z + 1/w + 1/(w*y)
F_387_mutated = w*x*y*z^2 + w*x*z + w*y*z^2 + w*y*z + w*z + w + x + z + 1/x + 1/(w*z) + 1/(w*y*z^2)
F_19_mutated = w*y*z + w*y + w*z + w + x + y + 1/(w^2*z) + 1/(w^2*x*z) + 1/(w^2*x*y*z)
F_28_mutated = w*y^2*z + w*y*z + w + x + y + z + 1/(y*z) + 1/(w*x*z) + 1/(w*x*y*z) + 1/(w*x*y^2*z^2)
F_33_mutated = w*x*y^2*z + w*x*y*z + w*y + w + z + 1/(w*z) + 1/(w*y*z) + 1/(w*x*y*z^2) + 1/(w*x*y^2*z^2)
F_37_mutated = w*x*y^2*z + w*x*y*z + w + z + y/(w*z) + 2/(w*z) + 1/(w*y*z) + 2/(w*x*y*z^2) + 2/(w*x*y^2*z^2) + 1/(w*x^2*y^3*z^3)
F_71_mutated = w*y*z^2 + 2*w*y*z + w*y + w*z^2 + 2*w*z + w + x + z + 1/z + 1/(w^2*x*y*z) + 2/(w^2*x*y*z^2) + 1/(w^2*x*y*z^3)
F_106_mutated = w*x*z + w*x + w*z + w + y + z + 1/z + 1/(w^2*x) + 1/(w^2*x*y*z)
F_108_mutated = w*x*z + w*x + w*z + w + y + y/z + z + 1/z + 1/(w^2*x) + 1/(w^2*x*y)
F_139_mutated = w*x + w*y + w + x + y + z + 1/z + z/(w*y) + z/(w*x) + z/(w*x*y)
F_148_mutated = w*y + w + x + y + z + 1/w + 1/(w*z) + 1/(w*y*z) + 1/(w*x)
F_168_mutated = w*x*y*z + w*x*y + w*z + w + x + y + z + 1/z + 1/(w*y) + 1/(w*y*z) + 1/(w*x) + 1/(w*x*z)
F_223_mutated = w*y*z + w*y + w*z + x + y + z + 1/z + 1/y + 1/(w*x*y*z) + 1/(w*x*y*z^2) + 1/(w*x*y^2*z)
F_275_mutated = w*y*z^2 + w*y*z + w*z + w + x + 1/x + x/(w*z) + 1/(w*z) + 1/(w^2*y*z^2)
F_317_mutated = w*y*z + w*y + w*z + w + x + z + 1/z + 1/x + 1/(w^2*y) + 2/(w^2*y*z) + 1/(w^2*y*z^2)
F_347_mutated = w*y*z + w*y + w*z + w + x + y + 1/w + 1/(w*z) + 1/(w*y*z) + 1/(w*x)
F_352_mutated = w*y*z^2 + w*y*z + w*z + w + w*y*z/x + w/x + x + 1/x + 1/(w*z) + 1/(w*y*z^2) + 1/(w^2*y*z^2)
F_388_mutated = w*x + w*y + w + x + y + z + 1/z + 1/(w*y) + 1/(w*x) + 1/(w*x*y)
F_431_mutated = w*z + w + x + y + z + 1/z + z^2/(x*y) + 4*z/(x*y) + 6/(x*y) + 4/(x*y*z) + 1/(x*y*z^2) + 1/w + 1/(w*z)
F_444_mutated = w*y + w + x + y + y/z + z + 1/z + y/x + 1/x + 1/w + 1/(w*y)
F_547_mutated = w*x*z + w*x + w*z + w + x + y + 1/y + 1/x + 1/w + 1/(w*z) + 1/(w*x) + 1/(w*x*z)
F_551_mutated = w*x*y^2*z^2 + w*x*y^2*z + w*x*y*z^2 + w*x*y*z + w*z + w + z + y/(w*z) + 2/(w*z) + 1/(w*y*z) + 2/(w*x*y*z^2) + 2/(w*x*y^2*z^2) + 1/(w*x^2*y^3*z^3)
F_35_mutated = w*y + w + x + y + z + 1/(w*z) + 1/(w*y*z) + 1/(w*x) + 1/(w*x*z) + 1/(w*x*y) + 1/(w*x*y*z)
F_39_mutated = w*x*z + w*x + w*y*z + w*y + w*z + w + x + y + 1/(w^2*y*z) + 1/(w^2*x*z) + 1/(w^2*x*y*z)
F_140_mutated = w*x*y*z + w*x*y + w*z + w + x + z + 1/z + z/(w*x) + 1/(w*x) + z/(w*x^2*y) + 1/(w*x^2*y)
F_145_mutated = w*x*y*z + w*x*y + w + x + z + 1/z + 2/(x*y) + 2/(x*y*z) + 1/(x^2*y^2*z) + z/(w*x) + 1/(w*x) + 1/(w*x^2*y)
F_169_mutated = w*x*y^2 + w*x*y + w*z + z + 1/z + y/w + 2/w + 1/(w*y) + 2*z/(w*x*y) + 2*z/(w*x*y^2) + z^2/(w*x^2*y^3)
F_171_mutated = w*x*z + w*y*z + w + x*z + x + y*z + y + z + 1/z + 1/(w*y) + 1/(w*x) + 1/(w*x*y*z)
F_221_mutated = w*x*y*z + w*x*y + w*z + x + y*z + y + z + z/x + z/(x*y) + 1/(w*x*y) + 1/(w*x*y*z) + 1/(w*x^2*y^2)
F_291_mutated = w*y*z + w*z + w + x + y + z + 1/z + 1/(w*z) + 1/(w*x*z) + 1/(w*x*y*z) + 1/(w*x*y*z^2)
F_294_mutated = w*y*z + w*y + w*z + x + y + y/z + z + 1/z + 1/y + 1/(w*x*z) + 1/(w*x*z^2) + 1/(w*x*y*z)
F_300_mutated = w*x*y*z + w*x + w*y*z^2 + w*y*z + w*z + w + x + 1/(w*z) + 1/(w*y*z^2) + 1/(w*x*z) + 1/(w*x*y*z^2) + 1/(w^2*y*z^2)
F_357_mutated = w*x*y + w*x + w*y + w + x*y + x + y + z + 1/z + 1/(y*z) + 1/w + 1/(w*x)
F_358_mutated = w*y*z^2 + w*y*z + w*z + w + x + 1/x + 1/(w*z) + 1/(w*y*z^2) + 1/(w*x) + 1/(w*x*y*z) + 1/(w^2*x*y*z^2)
F_363_mutated = w*x*y*z + w*y*z + w + x*y + x + y + z + 1/z + 1/(y*z) + 1/(w*y*z) + 1/(w*x*y*z) + 1/(w*x*y^2*z^2)
F_401_mutated = w*y*z^2 + w*y*z + w*z + w + w*y*z/x + w/x + x + 1/x + 1/(w*z) + 1/(w*y*z^2) + x/(w^2*y*z^2) + 1/(w^2*y*z^2)
F_404_mutated = w*x*y*z + w*x*y + w*z + x + y + z + 1/z + 2*z/(x*y) + 2/(x*y) + z/(x^2*y^2) + 1/(w*x*y) + 1/(w*x*y*z) + 1/(w*x^2*y^2)
F_405_mutated = w*x + w*y + w + x + x/z + y + y/z + z + 1/z + 1/(w*y) + 1/(w*x) + 1/(w*x*y)
F_465_mutated = w*x*y + w*x + w + x + y + y/z + z + 1/z + 1/x + 1/(x*z) + 1/(w*x) + 1/(w*x*y) + 1/(w*x^2*y)
F_466_mutated = w*x + w + x + y + y/z + z + 1/z + 1/x + 1/w + 1/(w*y) + 1/(w*x) + 1/(w*x*y)
F_467_mutated = w*y + w + x + x/z + y + z + 1/z + y/x + 1/x + 1/w + 1/(w*y) + 1/(w*x) + 1/(w*x*y)
F_502_mutated = w*y*z^2 + w*y*z + w*z + w + x + 1/x + 1/w + 1/(w*z) + 1/(w*y*z) + 1/(w*y*z^2) + 1/(w^2*y*z^2)
F_38_mutated = w*x*y*z + w*x + w*y*z + w + z + 1/(w*z) + 1/(w*y*z^2) + 1/(w*x*z) + 1/(w*x*y*z^2) + 1/(w^2*x*y*z) + 1/(w^2*x*y*z^2)
F_161_mutated = w*x*y*z + w*x + w*y^2*z + w*y + y*z + y + z + 1/y + 1/(w*y*z) + 1/(w*x*z)
F_178_mutated = w*x*y*z + w*x*y + w*x*z + w*x + w*z + w + x*y^2 + z + 1/z + z/y + 2/y + 1/(y*z) + z/(x*y) + 2/(x*y) + 1/(x*y*z) + y/w + y/(w*z)
F_361_mutated = w*x*z + w*x + w*y*z + w*y + w*z + w + 1/w + 1/(w*y*z) + 1/(w*x*z) + 1/(w*x*y*z) + 1/(w^2*x*y*z)
F_402_mutated = w*x^2*z^2 + w*x*z + x*z + x + y + z + 1/z + 1/y + 1/(x*z) + 1/(w*x) + 1/(w*x^2*z)
F_419_mutated = w*y*z^2 + w*y*z + w*z + w + w*y*z/x + w/x + x + 1/x + x/(w*z) + x/(w*y*z^2) + 1/(w*z) + 1/(w*y*z^2) + x/(w^2*y*z^2) + 1/(w^2*y*z^2)
F_519_mutated = w*x*y + w*x + w*y*z + w*y + w*z + w + x + x/z + z + 2/z + z/x + 2/x + 1/(x*z) + 1/(w^2*y)
F_522_mutated = w*x*y + w*x + w*y*z + w*y + w*z + w + x + y + 1/x + 1/w + 1/(w*z) + 1/(w*y*z) + 1/(w*x) + 1/(w^2*y*z) + 1/(w^2*x*y*z)
F_549_mutated = w*x*z + w*x + w*z + x + 2*x/z + x/z^2 + y + z + 2/z + 1/y + 1/x + 1/(w*z) + 1/(w*z^2) + 1/(w*x*z)
F_572_mutated = w*x + w*y + w + x + y + z + 1/z + x/(w*y) + 2/w + 2/(w*y) + y/(w*x) + 2/(w*x) + 1/(w*x*y)
F_477_mutated = w*y*z + w*y + w*z + w + x + 1/w + 1/(w*z) + 1/(w*y) + 1/(w*y*z) + 1/(w*x) + 1/(w*x*z) + 1/(w*x*y) + 1/(w*x*y*z) + 1/(w^2*y*z) + 1/(w^2*x*y*z)
F_480_mutated = w*x*y*z + w*x*y + w*x*z + w*x + w*z + w + x + y + 1/x + 1/(w*x) + 1/(w*x*z) + 1/(w*x*y*z) + 1/(w*x^2*y*z) + 1/(w^2*x^2*y*z)
F_486_mutated = w*x*y*z + w*x*y + w*x*z + w*x + w*z + w + x + y + 1/x + 1/w + 1/(w*x*z) + 1/(w*x*y*z) + 1/(w*x^2*y*z) + 1/(w^2*x*y*z)
F_525_mutated = w*x*z + w*x + w*z + x + 2*x/z + x/z^2 + y + z + 2/z + 1/x + 1/(w*z) + 1/(w*z^2) + 1/(w*y*z) + 1/(w*y*z^2) + 1/(w*x*z) + 1/(w*x*y*z)
F_528_mutated = w*x*z + w*x + w*z + x + 2*x/z + x/z^2 + y + z + 2/z + 1/x + 1/(w*z) + 1/(w*z^2) + 1/(w*y) + 1/(w*y*z) + 1/(w*x*z) + 1/(w*x*y)
F_533_mutated = w*x*z + w*x + w*z + w + x + x/z + y + y/z + z + 2/z + 1/x + 1/(x*z) + 1/(w*z) + 1/(w*y) + 1/(w*x*z) + 1/(w*x*y)
F_545_mutated = w*x*y^2 + w*x*y*z + w*x*y + w*x*z + w*y*z + w*z + x*y^2/z^2 + 2*x*y/z + x*y/z^2 + x + 2*x/z + x/y + 2*y/z + z + 2/z + 2/y + 1/x + 1/(x*y) + 1/(w*z) + 1/(w*y) + 1/(w*x*y)
F_634_mutated = w*x^2*y + w*x*y + w*x + w + x^2*y + 2*x*y + 2*x + y + z + 1/z + 1/y + 2/x + 2/(x*y) + 1/(x^2*y) + 1/(w*z) + 1/(w*x*z) + 1/(w*x*y*z) + 1/(w*x^2*y*z)
F_663_mutated = w*x^2*z + w*x*z + w*x + w + x^2*z + 2*x*z + 2*x + y + z + 1/z + 1/y + 2/x + 2/(x*z) + 1/(x^2*z) + 1/w + 1/(w*x) + 1/(w*x*z) + 1/(w*x^2*z)
F_672_mutated = w*x*z + w*x + w*z^2 + w*z + x + 2*x/z + x/z^2 + y + 2*z + 2/z + 1/y + z^2/x + 2*z/x + 1/x + 2/(w*z) + 2/(w*z^2) + 2/(w*x) + 2/(w*x*z) + 1/(w^2*x*z^2)
F_688_mutated = w*x*y^2*z + w*x*y^2 + w*x*y*z + w*x*y + w*z + w + x*y^3 + 3*x*y^2 + 3*x*y + x + 3*y + 3/y + 3/(x*y) + 3/(x*y^2) + 1/(x^2*y^3) + 1/(w^2*x*y*z) + 1/(w^2*x*y^2*z) + 1/(w^2*x^2*y^3*z)
F_495_mutated = w*x*z + w*x + w*y*z + w*y + w*z + w + x + y + 1/(w*y) + 1/(w*y*z) + 1/(w*x) + 1/(w*x*z) + 1/(w*x*y) + 1/(w*x*y*z) + 1/(w^2*x*y*z)
F_531_mutated = w*x*y + w*x + w*y + w + x*y/z + x + x/z + y + y/z + z + 1/z + y/x + z/x + 1/x + 1/(w*z) + 1/(w*y*z) + 1/(w*x) + 1/(w*x*y)
F_537_mutated = w*x^2*y + w*x^2 + 2*w*x*y*z + 2*w*x*y + 2*w*x*z + 2*w*x + w*y*z^2 + 2*w*y*z + w*y + w*z^2 + 2*w*z + w + x + x/z + z + 2/z + z/x + 2/x + 1/(x*z) + 1/(w^2*x^2*y*z^2)
F_577_mutated = w*x*y*z + w*x*y + w*x*z + w*y*z + x + y + 2*y/z + y/z^2 + z + 2/z + 1/y + 2*y/x + 2*y/(x*z) + 2/x + y/x^2 + 1/(w*x*z) + 1/(w*x*z^2) + 1/(w*x*y*z) + 1/(w*x^2*z)
F_579_mutated = w*x*y*z + w*x*z + w*x + w*y*z + x + y + z + 2/z + 1/y + 2/(y*z) + 1/(y*z^2) + 2*y/x + 2/x + 2/(x*z) + y/x^2 + 1/(w*x*z) + 1/(w*x*y*z) + 1/(w*x*y*z^2) + 1/(w*x^2*z)
F_604_mutated = w*x*y*z + w*x*y + w*x*z + w*y*z + w*y + w*z + x + y + 2*y/z + y/z^2 + z + 2/z + 1/y + 1/x + 1/(w*z) + 1/(w*z^2) + 1/(w*y*z) + 1/(w*x*z) + 1/(w*x*z^2) + 1/(w*x*y*z)
F_612_mutated = w*x*y + w*x + 2*w*y + 2*w + w*y/x + w/x + x + x/z + y + z + 1/z + 1/y + y/x + 2/x + 1/(x*y) + 1/w + 1/(w*z) + 1/(w*y) + 1/(w*y*z)
F_624_mutated = w*x*y + w*x + w*y*z + w*y + w*z + w + x^2/z + x + x/z + y + 1/y + y*z/x + y/x + 2*z/x + 2/x + z/(x*y) + 1/(x*y) + x/(w*z) + x/(w*y*z) + 1/w + 1/(w*z) + 1/(w*y) + 1/(w*y*z)
F_629_mutated = w*x*z + w*x + w*y*z + w*y + w*z + w + x + x/y + y + z + 1/z + 2/y + y/x + 2/x + 1/(x*y) + 1/(w*y) + 1/(w*y*z) + 1/(w*x) + 1/(w*x*z) + 1/(w*x*y) + 1/(w*x*y*z)
F_631_mutated = w*y*z + w*y + w*z + w + w*y*z/x + w*z/x + x + x/z + y + 1/z + 1/y + y*z/x + 2*z/x + 1/x + z/(x*y) + x/(w*z) + x/(w*y*z) + 1/w + 1/(w*z) + 1/(w*y) + 1/(w*y*z) + 1/(w*x) + 1/(w*x*y)
F_674_mutated = w*x + w*x/z + w*x/y + w*y + w*y/z + 3*w + 3*w*y/x + w*y^2/x^2 + x + x/z + 2*x/y + z + 1/z + 2/y + 2*y/x + 4/x + 2*y/x^2 + x/(w*y) + 1/w + 2/(w*y) + 2/(w*x) + 1/(w*x*y) + 1/(w*x^2)
F_693_mutated = w*x^2*y^4 + 2*w*x^2*y^3 + w*x^2*y^2 + 2*w*x*y^2 + 2*w*x*y + w + x*y^3 + 3*x*y^2 + 3*x*y + x + 3*y + z + 1/z + 3/y + 3/(x*y) + 3/(x*y^2) + 1/(x^2*y^3) + 1/(w*x*y*z) + 1/(w*x*y^2*z) + 1/(w*x^2*y^3*z)
F_703_mutated = w*x^2*y^4 + 2*w*x^2*y^3 + w*x^2*y^2 + 2*w*x*y^2 + 2*w*x*y + w + x*y^3 + 3*x*y^2 + 3*x*y + x + 3*y + z + 1/z + 3/y + 3/(x*y) + 3/(x*y^2) + 1/(x^2*y^3) + 1/(w*x*y) + 1/(w*x*y^2) + 1/(w*x^2*y^3)
F_149_mutated = w*x*y + w*x + w*y + w + y + z + 1/z + 1/(w^2*x) + 1/(w^2*x*y)
F_245_mutated = w*y + w + x + z + 1/z + 1/(y*z) + 1/w + 1/(w*x) + 1/(w*x*y)
F_443_mutated = w*y + w + x + y + y/z + z + 1/z + 1/x + 1/w + 1/(w*y)
F_327_mutated = w*x*z + w*y*z + w + x + y + z + 1/z + 1/(w*y*z) + 1/(w*x*z) + 1/(w*x*y*z^2)
F_311_mutated = w*x*y + w*x*y/z + w*x + w*x/z + w*y + w*y/z + x + y + z + y/x + 1/(w*x) + 1/(w*x*y) + 1/(w*x^2)
F_158_mutated = w*x*y^2 + w*x*y + w*y*z + w*z + z + 1/z + 1/w + 1/(w*y) + z/(w*x*y) + z/(w*x*y^2)
F_378_mutated = w*y + w + x + x/z + y + z + 1/z + 1/w + 1/(w*y) + 1/(w*x) + 1/(w*x*y)
F_400_mutated = w*z + w + x + x/y + y + z + 1/z + 1/w + 1/(w*z) + 1/(w*x) + 1/(w*x*z)
F_448_mutated = w*x + w + x + y + z + 1/z + 1/x + 1/w + 1/(w*y*z) + 1/(w*x) + 1/(w*x*y*z)
F_461_mutated = w*z + w + x + y + z + 1/z + 1/x + 1/w + 1/(w*z) + 1/(w*y) + 1/(w*y*z)
F_342_mutated = w*y^3*z^2 + w*y^2*z^2 + w*y^2*z + w*y*z + w*y^4*z^2/x + w*y^3*z/x + x + 2*y*z + z + 1/(y*z) + y^2*z/x + 1/(w*y) + 1/(w*y^2*z)
F_214_mutated = w*x*z + w*z + w + x + y + 1/(x*y) + x/w + 1/w + 1/(w*z)
F_248_mutated = w*x*z + w + x + x/y + y + 1/w + 1/(w*x) + 1/(w*x*z) + 1/(w*x^2*z)
F_107_mutated = w*x^2*y*z^2 + w*x*z + x*y*z + x + y + z + 1/(x*z) + 1/(x*y*z) + 1/(w*x) + 1/(w*x^2*y*z)
F_437_mutated = w*x*y^2*z + w*x*y^2 + w*x*y*z + w*x*y + w*y*z + w*z + y + z + 1/z + 2*z/(x*y) + 2/(x*y) + z/(x^2*y^2) + 1/(w*y) + 1/(w*y*z) + 1/(w*x*y^2)
F_445_mutated = w*z + w + x + y + 1/y + 1/x + x/w + x/(w*z) + 1/w + 1/(w*z)
F_500_mutated = w*x^2*z + w*x*z + w + y + 1/y + x/w + 2/w + 1/(w*x) + 1/(w*x*z) + 1/(w*x^2*z)
F_299_mutated = w*x*y^2*z^2 + w*x*y^2*z + w*x*y*z + w*y*z + w*y + w + y*z + z + 1/(w*y*z) + 1/(w*x*y^2*z^2)
F_188_mutated = w*x*y^3*z^2 + w*x*y^2*z^2 + w*x*y*z + w*y^2*z + w*y*z + w + y + z + 1/(w*y*z) + 1/(w*x*y^2*z^2)
F_262_mutated = w*x*y^2*z^2 + w*x*y*z + w*y*z^2 + w*y*z + w*z + w + y + z + 1/(w*y*z) + 1/(w*x*y^2*z) + 1/(w*x*y^2*z^2)
F_301_mutated = w*x*y^2*z^2 + w*x*y*z^2 + w*x*y*z + w*y*z + w*z + w + x*y^2*z + x*y*z + y*z + z + 1/(w*y*z) + 1/(w*x*y^2*z^2)

F_list_mutated = [F_730_mutated, F_729_mutated, F_725_mutated, F_717_mutated, F_712_mutated, F_709_mutated, F_699_mutated, F_679_mutated, F_710_mutated, F_678_mutated, \
                  F_708_mutated, F_689_mutated, F_671_mutated, F_707_mutated, F_666_mutated, F_676_mutated, F_704_mutated, F_660_mutated, F_673_mutated, F_659_mutated, \
                  F_701_mutated, F_653_mutated, F_670_mutated, F_34_mutated, F_159_mutated, F_164_mutated, F_308_mutated, F_524_mutated, F_176_mutated, F_538_mutated, \
                  F_618_mutated, F_570_mutated, F_651_mutated, F_3_mutated, F_17_mutated, F_22_mutated, F_26_mutated, F_30_mutated, F_32_mutated, F_46_mutated, \
                  F_58_mutated, F_59_mutated, F_60_mutated, F_63_mutated, F_68_mutated, F_111_mutated, F_115_mutated, F_124_mutated, F_128_mutated, F_210_mutated, \
                  F_216_mutated, F_228_mutated, F_233_mutated, F_256_mutated, F_261_mutated, F_279_mutated, F_284_mutated, F_302_mutated, F_309_mutated, F_326_mutated, \
                  F_341_mutated, F_343_mutated, F_350_mutated, F_368_mutated, F_387_mutated, F_19_mutated, F_28_mutated, F_33_mutated, F_37_mutated, F_71_mutated, \
                  F_106_mutated, F_108_mutated, F_139_mutated, F_148_mutated, F_168_mutated, F_223_mutated, F_275_mutated, F_317_mutated, F_347_mutated, F_352_mutated, \
                  F_388_mutated, F_431_mutated, F_444_mutated, F_547_mutated, F_551_mutated, F_35_mutated, F_39_mutated, F_140_mutated, F_145_mutated, F_169_mutated, \
                  F_171_mutated, F_221_mutated, F_291_mutated, F_294_mutated, F_300_mutated, F_357_mutated, F_358_mutated, F_363_mutated, F_401_mutated, F_404_mutated, \
                  F_405_mutated, F_465_mutated, F_466_mutated, F_467_mutated, F_502_mutated, F_38_mutated, F_161_mutated, F_178_mutated, F_361_mutated, F_402_mutated, \
                  F_419_mutated, F_519_mutated, F_522_mutated, F_549_mutated, F_572_mutated, F_477_mutated, F_480_mutated, F_486_mutated, F_525_mutated, F_528_mutated, \
                  F_533_mutated, F_545_mutated, F_634_mutated, F_663_mutated, F_672_mutated, F_688_mutated, F_495_mutated, F_531_mutated, F_537_mutated, F_577_mutated, \
                  F_579_mutated, F_604_mutated, F_612_mutated, F_624_mutated, F_629_mutated, F_631_mutated, F_674_mutated, F_693_mutated, F_703_mutated, F_149_mutated, \
                  F_245_mutated, F_443_mutated, F_327_mutated, F_311_mutated, F_158_mutated, F_378_mutated, F_400_mutated, F_448_mutated, F_461_mutated, F_342_mutated, \
                  F_214_mutated, F_248_mutated, F_107_mutated, F_437_mutated, F_445_mutated, F_500_mutated, F_299_mutated, F_188_mutated, F_262_mutated, F_301_mutated];


# Repaired only to nested non-degeneracy, not transverse

F_719_mutated = x^4*y^12*z*w^2 + 4*x^4*y^11*z*w^2 + 6*x^4*y^10*z*w^2 + 4*x^4*y^9*z*w^2 + x^4*y^8*z*w^2 + 4*x^3*y^9*z*w^2 + 12*x^3*y^8*z*w^2 + 12*x^3*y^7*z*w^2 + 4*x^3*y^6*z*w^2 + 6*x^2*y^6*z*w^2 + 2*x^2*y^6*z*w + 12*x^2*y^5*z*w^2 + 2*x^2*y^6*w + 4*x^2*y^5*z*w + 6*x^2*y^4*z*w^2 + 4*x^2*y^5*w + 2*x^2*y^4*z*w + 2*x^2*y^4*w + 4*x*y^3*z*w^2 + 4*x*y^3*z*w + 4*x*y^2*z*w^2 + x*y^4 + 4*x*y^3*w + 4*x*y^2*z*w + 3*x*y^3 + 4*x*y^2*w + 3*x*y^2 + z*w^2 + x*y + 2*z*w + 3*y + z + 2*w + z^-1 + 3*y^-1 + 3*x^-1*y^-2 + x^-1*y^-2*w^-1 + 3*x^-1*y^-3 + x^-1*y^-2*z^-1*w^-1 + x^-1*y^-3*w^-1 + x^-1*y^-3*z^-1*w^-1 + x^-2*y^-5 + x^-2*y^-5*w^-1 + x^-2*y^-5*z^-1*w^-1
F_690_mutated = x*y^2*z^5*w^2 + 2*x*y^2*z^4*w^2 + x*y^2*z^3*w^2 + 2*x*y*z^4*w^2 + 4*x*y*z^3*w^2 + 2*x*y*z^3*w + 2*x*y*z^2*w^2 + x*z^3*w^2 + 2*x*y*z^2*w + 2*x*z^2*w^2 + 2*x*z^2*w + 2*y*z^2*w + x*z*w^2 + y*z^2 + 2*x*z*w + 2*y*z*w + x*z + 2*y*z + 2*z*w + y + 2*z + 2*w + w^-1 + 2*z^-1 + y^-1 + z^-1*w^-1 + 2*y^-1*z^-1 + x^-1*z^-1 + y^-1*z^-1*w^-1 + x^-1*z^-1*w^-1 + y^-1*z^-2 + y^-1*z^-2*w^-1 + x^-1*z^-2*w^-1 + x^-1*y^-1*z^-2*w^-1 + x^-1*y^-1*z^-3*w^-1
F_677_mutated = w^2*(x + y)*(x*z + y*z + 1)^2/(x*y) + 2*w*(x + y)*(y + 1)*(x*z + y*z + 1)/(x*y) + (x^3*z^2 + 2*x^2*y*z^2 + 2*x^2*z + x*y^2*z^2 + x*y^2*z + x*z + x + y^3*z + 2*y^2*z + y*z)/(x*y*z) + (y + 1)*(x*z + y*z + 1)/(w*y*z)
F_686_mutated = w^2*y*z + w^2*y/x + w^2/x + w*x*y^2*z^2 + w*x*y*z^2 + w*y^2*z + 5*w*y*z + w*z + 2*w + 2*w*y/x + 2*w/x + 2*w/(x*z) + w/(x*y*z) + 2*x*y^2*z^2 + 2*x*y*z^2 + 2*x*y*z + x*z + y^2*z + 5*y*z + y + z + 1/(y*z) + y/x + 1/x + 2/(x*z) + 1/(x*y*z) + 1/(x*y*z^2) + x*y^2*z^2/w + x*y*z^2/w + 2*x*y*z/w + x*z/w + x/w + y*z/w + 2/w + 1/(w*y*z)
F_655_mutated = x*y^3*z^-1*w^2 + 2*y^4*z^-1*w^2 + x^-1*y^5*z^-1*w^2 + 2*x*y^2*z^-1*w^2 + 4*y^3*z^-1*w^2 + 2*x^-1*y^4*z^-1*w^2 + 2*y^2*w + 2*x^-1*y^3*w + x*y*z^-1*w^2 + 2*y^2*z^-1*w^2 + x^-1*y^3*z^-1*w^2 + 2*y^2*z^-1*w + 2*x^-1*y^3*z^-1*w + 2*y*w + 2*x^-1*y^2*w + x + y + x^-1*y*z + 2*y*z^-1*w + 2*x^-1*y^2*z^-1*w + 2*x*y^-1 + 2*x^-1*y + y^-1*z*w^-1 + x^-1*y*z^-1 + x*y^-2 + y^-1 + y^-1*w^-1 + y^-2*z*w^-1 + y^-2*w^-1
F_614_mutated = w*x*z + w*x + w + w/(x^2*y*z^2) + x^2*y*z^3 + 2*x^2*y*z^2 + x^2*y*z + 2*x*z + 2*x + 1/x + 1/(x*z) + 1/(x^2*y*z^2) + x^2*y*z^3/w + 2*x^2*y*z^2/w + x^2*y*z/w + x*z/w + x/w + 1/w
F_585_mutated = w + x + x/y + y + z + 2*z/y + 2/y + y/x + 2*z/x + 2/x + z^2/(x*y) + 2*z/(x*y) + 1/(x*y) + x/w + x/(w*z) + y/w + y/(w*z) + z/w + 2/w + 1/(w*z)
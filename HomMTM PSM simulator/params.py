import math
import numpy as np

# the error here is (inten - theo)/theo

#U2OS(LabelFree) (total theo peak, miss peak):(27432, 3906)
#UCEC  (total theo peak, miss peak):(108628, 8575)

peak_miss_prob = (3906 + 8575) / ((27432 + 108628) * 2)
# peak_miss_prob = 0
mu_not_shared = -0.351019046350688
sigma_not_shared = math.sqrt(0.12801875919212194)


# mu_shared_error = 1.7539267814265496
# sigma_shared_error = math.sqrt(21.40564231561222)

mu_shared = 0.9927565875626165
sigma_shared = math.sqrt(6.561070222079408)

mu_max_shared = 9.906421076409224
sigma_max_shared = math.sqrt(107.72946099201894)



#small error:
# mu_not_shared = 0
# sigma_not_shared = 0.1
#
# mu_shared = 0
# sigma_shared = 0.1
#
# mu_max_shared = 0
# sigma_max_shared = 0.1


max_err = 20

random_abundance_upperbound = 10
charge = 2
WaterMass = 18.0105
ProtonMass = 1.007276
residue_mass = {"A":71.037114,
                "R":156.101111,
                "N":114.042927,
                "D":115.026943,
                "C":103.009185,
                "E":129.042593,
                "Q":128.058578,
                "G":57.021464,
                "H":137.058912,
                "I":113.084064,
                "L":113.084064,
                "K":128.094963,
                "M":131.040485,
                "F":147.068414,
                "P":97.052764,
                "S":87.032028,
                "T":101.047679,
                "U":150.95363,
                "W":186.079313,
                "Y":163.06332,
                "V":99.068414}

mods_mass = {"phosphorylation": 79.9663,
             "Oxidation":15.994915}


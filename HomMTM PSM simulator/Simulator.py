import math

import ReadRealCases
import numpy as np
import params
import random
import re

random.seed(66)
np.random.seed(66)

def getErr(mu,sigma):
    err = np.random.normal(mu, sigma)
    while err < -1 or err > params.max_err:
        err = np.random.normal(mu, sigma)
    return err

def getMaxSharedErr(mu,sigma):
    err = np.random.normal(mu, sigma)
    while err < 0 or err > params.max_err:
        err = np.random.normal(mu, sigma)
    return err

def mods_mass_pos_list(mods):
    mods_list = re.findall(r'<(.*?)>', mods)
    mass_pos_list = []
    for mod in mods_list:
        name = mod.split(",")[0]
        pos = int(mod.split(",")[1])
        mass = params.mods_mass[name]
        mass_pos_list.append((mass,pos))
    return mass_pos_list

def get_masses_mod(masses_no_mod, shift_pos_list):
    masses_mod = masses_no_mod
    for mod in shift_pos_list:
        shift = mod[0]
        pos = mod[1]
        masses_mod[pos] = masses_mod[pos] + shift
    return masses_mod

def getResMasses(pep_backbone, mods):
    masses_no_mod = [0]
    for res in pep_backbone:
        masses_no_mod.append(params.residue_mass[res])

    shift_pos_list = mods_mass_pos_list(mods)

    masses_mod = get_masses_mod(masses_no_mod, shift_pos_list)


    return masses_mod



def getTheoPeakMass(res_masses):
    peak_mass_list = []
    peak_mass = 0
    for i in range(len(res_masses) - 1):
        peak_mass = round(peak_mass + res_masses[i], 4)
        if i != 0:
            peak_mass_list.append(round(peak_mass + params.ProtonMass, 4))
    return peak_mass_list

# print(getTheoPeakMass(getResMasses("HAEPSSSPSK", "")))
# print(getTheoPeakMass(getResMasses("HAEPSSSPSK", "<phosphorylation,5>")))
# print(getTheoPeakMass(getResMasses("HAEPSSSPSK", "<phosphorylation,6>")))
# print(getTheoPeakMass(getResMasses("HAEPSSSPSK", "<phosphorylation,7>")))

def Missing(pep_backbone):
    peaks_missing = []
    for i in range(len(pep_backbone) - 1):
        prob = random.random()
        if prob < params.peak_miss_prob:
            peaks_missing.append(i)
    return peaks_missing

def getSimSpectrumMasses(pep_backbone, mods):
    res_masses = getResMasses(pep_backbone, mods)
    prec_mass = sum(res_masses) + params.WaterMass
    prec_mz = prec_mass / params.charge + params.ProtonMass
    theo_peak_mass_list = getTheoPeakMass(res_masses)


    return theo_peak_mass_list, prec_mass, prec_mz

def getMixSpectrum(peak_mass_1, peak_mass_2, q_1, q_2):
    shared_peaks = []
    not_shared_peaks = []
    errors = []

    shared_idx = []
    for i in range(len(peak_mass_1)):
        if peak_mass_1[i] in peak_mass_2 and peak_mass_1[i] != -1:
            shared_idx.append(i)

    max_shared_idx = random.sample(shared_idx,1)[0]
    max_shared_inten = (q_1 + q_2) * (1 + getMaxSharedErr(params.mu_max_shared, params.sigma_max_shared))
    while(max_shared_inten <= q_1 + q_2):
        max_shared_inten = (q_1 + q_2) * (1 + getMaxSharedErr(params.mu_max_shared, params.sigma_max_shared))
    shared_peaks.append((peak_mass_1[max_shared_idx], max_shared_inten))
    errors.append(abs(max_shared_inten - q_1 - q_2))

    for i in range(len(peak_mass_1)):
        if peak_mass_1[i] == -1:
            errors.append(q_1)

    for i in range(len(peak_mass_2)):
        if peak_mass_2[i] == -1:
            errors.append(q_2)

    for i in range(len(peak_mass_1)):
        if peak_mass_1[i] != -1:
            mass = peak_mass_1[i]
            if peak_mass_1[i] in peak_mass_2:
                if i != max_shared_idx:
                    inten = (q_1 + q_2) * (1 + getErr(params.mu_shared, params.sigma_shared))
                    while inten >= max_shared_inten:
                        inten = (q_1 + q_2) * (1 + getErr(params.mu_shared, params.sigma_shared))
                    shared_peaks.append((mass,inten))
                    errors.append(abs(inten - q_1 -q_2))
            else:
                inten = q_1 * (1 + getErr(params.mu_not_shared, params.sigma_not_shared))
                while inten >= max_shared_inten:
                    inten = q_1 * (1 + getErr(params.mu_not_shared, params.sigma_not_shared))
                not_shared_peaks.append((mass,inten))
                errors.append(inten - q_1)
    for i in range(len(peak_mass_2)):
        if peak_mass_2[i] != -1:
            mass = peak_mass_2[i]
            if peak_mass_2[i] not in peak_mass_1:
                inten = q_2 * (1 + getErr(params.mu_not_shared, params.sigma_not_shared))
                while inten >= max_shared_inten:
                    inten = q_2 * (1 + getErr(params.mu_not_shared, params.sigma_not_shared))
                not_shared_peaks.append((mass,inten))
                errors.append(inten - q_2)

    mix_peaks = shared_peaks + not_shared_peaks
    mix_peaks.sort(key = lambda x: x[0])
    return mix_peaks, sum(errors)

def RandomAbundance():
    a = random.randint(1, params.random_abundance_upperbound - 1)
    b = params.random_abundance_upperbound - a

    abund_1 = float(a/(a + b))
    abund_2 = float(b/(a + b))

    return abund_1, abund_2

def noise_at_sound(noise_mass, sound_peaks, prec_mass):
    for peak in sound_peaks:
        if abs(noise_mass - peak) < 0.11 or abs((prec_mass - peak + params.ProtonMass) - (peak - params.ProtonMass)) < 0.11:
            return True
    return False

def genNoisePeak(inten_mean, inten_var, min_inten, max_inten, prec_mass):
    noise_mass = np.random.normal(prec_mass/2, prec_mass/4)
    while noise_mass < 10 or noise_mass > prec_mass - 10:
        noise_mass = np.random.normal(prec_mass/2, prec_mass/4)



    noise_inten = np.random.normal(inten_mean, math.sqrt(inten_var))
    while noise_inten < min_inten or noise_inten > max_inten:
        noise_inten = np.random.normal(inten_mean, math.sqrt(inten_var))

    return (round(noise_mass,5), round(noise_inten,5))

def genNoise(case, prec_mass, theo_mass_1, theo_mass_2):
    mass_mean = case["Noise Mass Mean"]
    mass_var = case["Noise Mass Var"]
    inten_mean = case["Noise Inten Mean"]
    inten_var = case["Noise Inten Var"]
    noise_num = int(case["Number of Noise"] / 1.2)
    min_noise_inten = case["min Noise Inten"]
    max_noise_inten = case["max Noise Inten"]

    # max_inten_peak = max(mix_peaks,key=lambda x:x[1])
    # max_noise_inten = max_inten_peak[1]

    theo_mass = theo_mass_1 + theo_mass_2

    noise = []
    for i in range(noise_num):
        noisepeak = genNoisePeak(inten_mean, inten_var, min_noise_inten, max_noise_inten, prec_mass)
        if not noise_at_sound(noisepeak[0], theo_mass, prec_mass):
            noise.append(noisepeak)

    return noise

def Simulate(case):
    pep_backbone = case["Pep Backbone"]
    # abund_1, abund_2 = RandomAbundance()

    abund_1, abund_2 = 0.2, 0.8

    mods_1 = case["Mod 1"]
    mods_2 = case["Mod 2"]
    theo_iten_sum = case["Theo Inten Sum"]

    q_1 = round(theo_iten_sum * abund_1,5)
    q_2 = round(theo_iten_sum * abund_2,5)

    theo_mass_1, prec_mass, prec_mz  = getSimSpectrumMasses(pep_backbone, mods_1)

    theo_mass_2, _, _ = getSimSpectrumMasses(pep_backbone, mods_2)

    miss_peak_idx = Missing(pep_backbone)

    peak_mass_1 = theo_mass_1
    peak_mass_2 = theo_mass_2
    for idx in miss_peak_idx:
        if abs(peak_mass_1[idx] - peak_mass_2[idx]) < 0.001:
            peak_mass_1[idx] = -1
            peak_mass_2[idx] = -1
        else:
            p = random.random()
            if p > 0.5:
                peak_mass_1[idx] = -1
            else:
                peak_mass_2[idx] = -1


    miss_1 = []
    miss_2 = []
    for i, mass in enumerate(peak_mass_1):
        if mass == -1:
            miss_1.append(i + 1)
    for i, mass in enumerate(peak_mass_2):
        if mass == -1:
            miss_2.append(i + 1)

    mix_peaks, error= getMixSpectrum(peak_mass_1, peak_mass_2, q_1, q_2)


    offset_mix_peaks = []
    for peak in mix_peaks:
        mass = peak[0]
        offset = random.uniform(-0.05, 0.05)

        new_mass = round(mass + offset, 5)
        offset_mix_peaks.append((new_mass, peak[1]))

    noise_peaks = genNoise(case, prec_mass, theo_mass_1, theo_mass_2)
    # noise_peaks = []
    noise_peaks.sort(key = lambda x: x[0])
    simulated_spectrum = offset_mix_peaks + noise_peaks
    simulated_spectrum.sort(key = lambda x: x[0])

    return simulated_spectrum, prec_mass, prec_mz, abund_1, abund_2, q_1, q_2, miss_1, miss_2, error






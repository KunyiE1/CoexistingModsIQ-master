import ReadRealCases
import random
import Simulator
import writeFiles
import params

if __name__ == "__main__":
    simulation_num = 1000
    cases1 = ReadRealCases.readCases("./UCECTMT_isoformPSM.pkl")
    cases2 = ReadRealCases.readCases("./LabelFree_isoformPSM.pkl")
    cases_selected_1 = random.sample(cases1, int(simulation_num/2))
    cases_selected_2 = random.sample(cases2, int(simulation_num/2))
    cases_selected = cases_selected_1 + cases_selected_2
    sp_list = []
    for case in cases_selected:

        sp = {}
        sim_peaks, prec_mass, prec_mz, abund_1, abund_2, q_1, q_2, miss_1, miss_2, error = Simulator.Simulate(case)
        sp["peaks"] = sim_peaks
        sp["prec_mass"] = prec_mass
        sp["prec_mz"] = prec_mz
        sp["raw_seq"] = case["Pep Backbone"]
        sp["proteoforms"] = [case["Mod 1"], case["Mod 2"]]
        sp["Abundance"] = (abund_1, abund_2)
        sp["theory intensities"] = (q_1, q_2)
        sp["charge"] = params.charge
        sp["miss peak"] = [miss_1, miss_2]
        sp["error"] = error
        sp_list.append(sp)
        print(len(sp_list))
    writeFiles.toFiles(sp_list, "./sim_sp_diff_frac/28seed66_2/")
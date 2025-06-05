def toFiles(spectra, out_name):
    sp_list = []
    for spectrum in spectra:
        mass_inten = spectrum["peaks"]
        peaks = []
        for i in range(len(mass_inten)):
            peaks.append(str(mass_inten[i][0]) + ' ' + str(round(mass_inten[i][1],5)))

        sp = {'ID': len(sp_list),
              'FRACTION_ID': -1,
              'FILE_NAME': '',
              'SCANS': len(sp_list),
              'RETENTION_TIME': -1,
              'LEVEL': 2,
              'ACTIVATION': "HCD",
              'MS_ONE_ID': -1,
              'MS_ONE_SCAN': -1,
              'PRECURSOR_MZ': spectrum["prec_mz"],
              'PRECURSOR_CHARGE': spectrum["charge"],
              'PRECURSOR_MASS': spectrum["prec_mass"],
              'PRECURSOR_INTENSITY': -1,
              'PEAKS': peaks,
              'Fraction_ID': -1,
              'Fraction_feature_ID': -1,
              'Fraction_feature_intensity': -1,
              'Fraction_feature_score': -1,
              'Sample_feature_ID': -1,
              'Sample_feature_intensity': -1,
              'raw_seq': spectrum['raw_seq'],
              'proteoforms': spectrum['proteoforms'],
              'abundance':spectrum['Abundance'],
              'theory intensities': spectrum['theory intensities'],
              'missing peaks': spectrum["miss peak"],
              'error': spectrum["error"]}

        sp_list.append(sp)

    f_msalign = open(out_name + "sim_ms2.msalign", 'w', encoding='utf-8')
    for i in range(len(sp_list)):
        f_msalign.write('BEGIN IONS\n')
        f_msalign.write("ID=" + str(sp_list[i]['ID']) + "\n")
        f_msalign.write("FRACTION_ID=" + str(sp_list[i]['FRACTION_ID']) + "\n")
        f_msalign.write("FILE_NAME=" + str(sp_list[i]['FILE_NAME']) + "\n")
        f_msalign.write("SCANS=" + str(sp_list[i]['SCANS']) + "\n")
        f_msalign.write("RETENTION_TIME=" + str(sp_list[i]['RETENTION_TIME']) + "\n")
        f_msalign.write("LEVEL=" + str(sp_list[i]['LEVEL']) + "\n")
        f_msalign.write("ACTIVATION=" + str(sp_list[i]['ACTIVATION']) + "\n")
        f_msalign.write("MS_ONE_ID=" + str(sp_list[i]['MS_ONE_ID']) + "\n")
        f_msalign.write("MS_ONE_SCAN=" + str(sp_list[i]['MS_ONE_SCAN']) + "\n")
        f_msalign.write("PRECURSOR_MZ=" + str(sp_list[i]['PRECURSOR_MZ']) + "\n")
        f_msalign.write("PRECURSOR_CHARGE=" + str(sp_list[i]['PRECURSOR_CHARGE']) + "\n")
        f_msalign.write("PRECURSOR_MASS=" + str(sp_list[i]['PRECURSOR_MASS']) + "\n")
        f_msalign.write("PRECURSOR_INTENSITY=" + str(sp_list[i]['PRECURSOR_INTENSITY']) + "\n")
        peak_list = sp_list[i]['PEAKS']
        for j in range(len(peak_list)):
            # print(peak_list[j].split(' ')[0] + '\t' + peak_list[j].split(' ')[1] + '\t')

            f_msalign.write(peak_list[j].split(' ')[0] + '\t')
            f_msalign.write(peak_list[j].split(' ')[1] + '\t')
            f_msalign.write('1' + '\n')
        f_msalign.write("END IONS\n")
        f_msalign.write('\n')
    f_msalign.close()

    f_feature = open(out_name + "sim_ms2.feature", 'w', encoding='utf-8')
    f_feature.write("Spec_ID\tFraction_ID\tFile_name\tScans\tMS_one_ID\tMS_one_scans\tPrecursor_mass\tPrecursor_intensity\tFraction_feature_ID\tFraction_feature_intensity\tFraction_feature_score\tSample_feature_ID\tSample_feature_intensity")
    f_feature.write("\n")
    for i in range(len(sp_list)):
        f_feature.write(str(sp_list[i]['ID']) + "\t"
                        + str(sp_list[i]['FRACTION_ID']) + "\t"
                        + str(sp_list[i]['FILE_NAME']) + "\t"
                        + str(sp_list[i]['SCANS']) + "\t"
                        + str(sp_list[i]['MS_ONE_ID']) + "\t"
                        + str(sp_list[i]['MS_ONE_SCAN']) + "\t"
                        + str(sp_list[i]['PRECURSOR_MASS']) + "\t"
                        + str(sp_list[i]['PRECURSOR_INTENSITY']) + "\t"
                        + str(sp_list[i]['Fraction_feature_ID']) + "\t"
                        + str(sp_list[i]['Fraction_feature_intensity']) + "\t"
                        + str(sp_list[i]['Fraction_feature_score']) + "\t"
                        + str(sp_list[i]['Sample_feature_ID']) + "\t"
                        + str(sp_list[i]['Sample_feature_intensity']) + "\n")
    f_feature.close()

    f_ref_pep = open(out_name + "ref_peptide.txt", 'w', encoding='utf-8')
    for i in range(len(sp_list)):
        f_ref_pep.write(str(sp_list[i]['ID']) + "\t"
                        + str(sp_list[i]['raw_seq']) + "\n")
    f_ref_pep.close()

    f_gt = open(out_name + "gt.txt", 'w', encoding='utf-8')
    f_gt.write("ID\tError\tMod1\tMod2\tAbund1\tAbund2\tq1\tq2\tMissPeaks1\tMissPeaks2\n")
    for i in range(len(sp_list)):
        f_gt.write(str(sp_list[i]["ID"]) + "\t"
                   + str(sp_list[i]["error"]) + "\t"
                   + sp_list[i]["proteoforms"][0] + "\t"
                   + sp_list[i]["proteoforms"][1] + "\t"
                   + str(sp_list[i]["abundance"][0]) + "\t"
                   + str(sp_list[i]["abundance"][1]) + "\t"
                   + str(sp_list[i]['theory intensities'][0]) + "\t"
                   + str(sp_list[i]['theory intensities'][1]))
        for miss_peaks in sp_list[i]['missing peaks']:
            f_gt.write("\t[")
            for q, p in enumerate(miss_peaks):
                if q != 0:
                    f_gt.write(',')
                f_gt.write(str(p))
            f_gt.write(']')
        f_gt.write("\n")
    f_gt.close()

from pyopenms import *

debug = 0

# read mzML
mzML_path = "/Users/alka/Documents/work/projects/OpenSWATH_Metabolomics_data/20181121_full_data_mzml/04_PestMixes_individually_Solvent_DDA_20-50/03_rt_truncated/PestMix1_1ngSolventDDA20-50.mzML"
exp = MSExperiment()
MzMLFile().load(mzML_path, exp)

# FeatureFinderMetabo
# used for feature detection

# masstrace detection
# https://github.com/OpenMS/OpenMS/blob/develop/src/topp/FeatureFinderMetabo.cpp#L234

# TODO: set paramters - see ffm_param
m_traces = list()
mtdet = MassTraceDetection()
mtdet.run(exp, m_traces, 10000000000000) # third agrument - max_traces

# elution peak detection
# https://github.com/OpenMS/OpenMS/blob/develop/src/topp/FeatureFinderMetabo.cpp#L245

#TODO: set parameters - see ffm_param
splitted_traces = list()
m_traces_final = list()
epdet = ElutionPeakDetection()
epdet.detectPeaks(m_traces, splitted_traces)
m_traces_final = splitted_traces

# feature finding
# https://github.com/OpenMS/OpenMS/blob/develop/src/topp/FeatureFinderMetabo.cpp#L284

ffm_param = Param()
ffm_param = FeatureFindingMetabo().getParameters()
ffm_param.setValue("report_convex_hulls", 'true')

# check parameters
#for elem in ffm_param.items():
#    print(elem)

fmap = FeatureMap()
fchrom = list()
ffmet = FeatureFindingMetabo()
ffmet.setParameters(ffm_param)

ffmet.run(m_traces_final, fmap, fchrom)

print("done - FFM")

# HighResPrecusorMassCorrector
# to correct ms2 in the feature space to the monoisotopic trace for easier mapping
# confex_hull has to be reported in FFM
# https://github.com/OpenMS/OpenMS/blob/develop/src/topp/HighResPrecursorMassCorrector.cpp#L194

# TODO: set parameters
corrected_precursors = list()
rt_tolerance = 0.0
mz_tolerance = 100.0
mz_unit_ppm = True # true = ppm
believe_charge = False
keep_original = False
assign_all_matching = False
max_trace = 3 # Maximum isotopic trace considered in matching a precursor to a feature
debug_level_ = 0

corrected_to_nearest_feature = PrecursorCorrection().correctToNearestFeature(fmap,
                                                                             exp,
                                                                             rt_tolerance,
                                                                             mz_tolerance,
                                                                             mz_unit_ppm,
                                                                             believe_charge,
                                                                             keep_original,
                                                                             assign_all_matching,
                                                                             max_trace,
                                                                             debug_level_)

print("done - HighResPrecursorMassCorrection")

# Feature Mapping
# mapping of MS2 Indices over precursor information (mz/rt) and feature information
# https://github.com/OpenMS/OpenMS/blob/develop/src/openms/source/ANALYSIS/MAPMATCHING/FeatureMapping.cpp#L43
# https://github.com/OpenMS/OpenMS/blob/develop/src/openms/source/ANALYSIS/ID/SiriusAdapterAlgorithm.cpp#L236
#https://github.com/jpfeuffer/OpenMS/blob/feature/MetFragMSPFile/src/utils/MetFragAdapter.cpp#L171

list_fmap = list()
fp_map_kd = KDTreeFeatureMaps()

list_fmap.append(fmap)
fp_map_kd.addMaps(list_fmap)

precursor_mz_tol = 25
precursor_rt_tol = 7 # seconds
ppm_prec = True # true = ppm

feature_mapping = FeatureMapping().assignMS2IndexToFeature(exp,
                                                           fp_map_kd,
                                                           precursor_mz_tol,
                                                           precursor_rt_tol,
                                                           ppm_prec)

print("done - FeatureMapping")

if debug >= 2:
    om_path = "test_mzML.mzML"
    MzMLFile().store(om_path, exp)
    of_path = "test_feaureXML.featureXML"
    FeatureXMLFile().store(of_path, fmap)


# MetaboTargetedAssay
# method to extract a potential transistions based on the ms/ms based of the highest intensity precursor or a consensus spectrum
# https://github.com/OpenMS/OpenMS/blob/develop/src/openms/source/ANALYSIS/TARGETED/MetaboTargetedAssay.cpp

precursor_rt_tol = 7.0
precursor_mz_distance = 0.0001 # SpectraMerging: Max m/z distance of the precursor entries of two spectra to be merged in [Da]
cosine_sim_threshold = 0.98 # SpectraMerging: Threshold for cosine similarity of MS2 spectra from the same precursor used in consensus spectrum creation
transition_threshold = 0.00000000001 # Intensity threshold for MS2 peak used in MetaboTargetedAssay (% of max intensity of the spectrum)
min_fragment_mz = 0.0
max_fragment_mz = 2000.0
method_consensus_spectrum = True
exclude_ms2_precursor = False
file_counter = 0

print("debug_to_here")

tmp_mta = list()
tmp_mta = MetaboTargetedAssay().extractMetaboTargetedAssay(exp,
                                                           feature_mapping,
                                                           precursor_rt_tol,
                                                           precursor_mz_distance,
                                                           cosine_sim_threshold,
                                                           transition_threshold,
                                                           min_fragment_mz,
                                                           max_fragment_mz,
                                                           method_consensus_spectrum,
                                                           exclude_ms2_precursor,
                                                           file_counter)


print("done - TargetedExtraction")



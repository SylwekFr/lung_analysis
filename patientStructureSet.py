import pydicom
import os
import numpy as np
import collections
import statistics
import operator
import SimpleITK as sitk
from viewer import *
import re

class PatientStructureSet:

    def __init__(self, patient):
        self.patient = patient
        self.before_RT_structure_set_data_files_paths = patient.before_RT_data_files_paths
        self.structure_set_file = self.get_structure_set_file()
        self.all_slices_contours = self.get_contours_for_all_slices()


    def get_structure_set_file(self):
        for file in self.before_RT_structure_set_data_files_paths:
            if re.search("RS", file):
                 return pydicom.read_file(file)


    def get_contours_for_all_slices(self):
        avaliable_rois = {}
        all_roi_contours = {}

        for roi in self.structure_set_file.StructureSetROISequence:
            avaliable_rois[roi.ROINumber] = roi.ROIName

        for key, slice in self.patient.before_slices.items():
            all_slice_roi_contours = []
            for roi in self.structure_set_file.ROIContourSequence:
                one_slice_roi_contour = {}
                one_slice_roi_contour['color'] = roi.ROIDisplayColor
                one_slice_roi_contour['referenceROIName'] = avaliable_rois[roi.ReferencedROINumber]
                one_slice_roi_contour['contourPoints'] = []
                for contour_sequence in roi.ContourSequence:
                    if contour_sequence.ContourImageSequence[0].ReferencedSOPInstanceUID == slice['ID']:
                        slice_contour = [[contour_sequence.ContourData[i], contour_sequence.ContourData[i + 1]] for i in
                                         range(0, len(contour_sequence.ContourData), 3)]
                        one_slice_roi_contour['contourPoints'].append(slice_contour)
                all_slice_roi_contours.append(one_slice_roi_contour)
            all_roi_contours[key] = all_slice_roi_contours

        return all_roi_contours










from operator import itemgetter

import numpy
import math
import statistics
from os import system, name
from patient import Patient
from patientRtDose import PatientRtDose
import numpy as np
import SimpleITK as sitk
import pandas as pd
import time
import csv
from SeriesType import SeriesType

from matplotlib import pyplot

from patientStructureSet import PatientStructureSet
from viewer import *
from utils import *


class PreprocessedImagesViewer:

    def __init__(self, segmented_images, registered_images, patient: Patient):


        self.segmented_before_images = segmented_images[0]
        self.segmented_after_images = segmented_images[1]
        self.registered_images = registered_images
        self.patient = patient
        self.patientRtDose = patient.patient_RT_dose
        self.patientStructureSet = patient.patient_structure_set
        # self.elastic_registered_image = np.load("patient\\elastic_bspline_registration_fixed\\moving_resampled_arrays\\moving_elastic_array_33.npy")
        self.slice_ROIs = []


        ###########################################################################

        # DISPLAY RAW DICOM IMAGE (slice)

        # self._display_raw_dicom_slice(self.elastic_registered_image)
        # self._display_raw_dicom_slice(self.segmented_before_images[33])


        ###########################################################################

        # DISPLAY DICOM IMAGE WITH CHOOSEN ROI's (choosen ROIs list, slice number)

        # chosenROIs = ['Pluco (L)']
        # self.display_dicom_with_ROI(chosenROIs, 62, self.segmented_before_images[62])
        # self.display_dicom_with_ROI(chosenROIs, 33, self.elastic_registered_image)
        # self.display_dicom_with_ROI(chosenROIs, 60, self.patient.before_slices_image[60])


        ###########################################################################

        # DISPLAY DICOM IMAGE WITH RADIATION IN RANGE FOR GIVEN IMAGE (slice, low limit, high limit)

        # self.display_dicom_with_radiation_in_range_for_given_image(33, 10, 30, self.segmented_before_images[33])
        # self.display_dicom_with_radiation_in_range_for_given_image(33, 10, 50, self.elastic_registered_image)
        # self.display_dicom_with_radiation_in_range_for_given_image(33, 10, 30, self.patient.before_slices_image[33])

        ###########################################################################

        # DISPLAY DICOM IMAGE WITH FULL RADIATION AS HEATMAP (slice number)

        # self.display_dicom_with_full_radiation_as_heat_map(33, self.elastic_registered_image)
        # self.display_dicom_with_full_radiation_as_heat_map(33, self.patient.before_slices_image[33])
        # self.display_dicom_with_full_radiation_as_heat_map(33, self.segmented_before_images[33])

        ###########################################################################

        # DISPLAY AND CALCULATE VALUES IN SEGMENTED LUNGS IN CHOSEN RANGE BUT WITHOUT ADDITIONAL ROI

        # print("Only radiation range 10-30 - First")
        # self.display_and_calculate_dicom_with_radiation_in_range_without_roi(self.segmented_before_images[33], 33, 2, 5)
        # print("Only radiation range 10-30 - Second")
        # self.display_and_calculate_dicom_with_radiation_in_range_without_roi(self.elastic_registered_image, 33, 10, 30)

        # DISPLAY AND CALCULATE VALUES IN SEGMENTED LUNGS IN CHOSEN RANGE AND INSIDE CHOSEN ROI

        # print("Only ROI's - First")

        radiation_range = [[1,5],[5,10],[10,20],[20,30],[30,40]]
        #chosenROIs = ['Pluco (P)', 'Dmin']
        chosenROIs = ['Pluca']
        for radrange in radiation_range:
            csvBeforeData = []
            csvAfterData = []
            slices_count = 0
            for i, slice in enumerate(self.segmented_before_images):
                slices_count += 1
                print("Slice before: " + str(i))
                csv_line, calculated_values = self.display_and_calculate_dicom_with_radiation_in_range_inside_common_part(slice, i, radrange[0], radrange[1],
                                                                                                                          chosenROIs, SeriesType.BEFORE.name)
                if None not in calculated_values:
                    csvBeforeData.append(csv_line)
            patient_before_data_df = pd.DataFrame(csvBeforeData, columns = ['Series a', 'Slice Index', 'Chosen ROIs', 'Radiation Range', 'Median a', 'Mean a'])
            slices_count = 0
            for i, slice in enumerate(self.segmented_after_images):
                slices_count += 1
                print("Slice after: " + str(i))
                csv_line, calculated_values = self.display_and_calculate_dicom_with_radiation_in_range_inside_common_part(slice, i, radrange[0], radrange[1],
                                                                                                                          chosenROIs, SeriesType.AFTER.name)
                if None not in calculated_values:
                    csvAfterData.append(csv_line)
            patient_after_data_df = pd.DataFrame(csvAfterData,
                                                  columns=['Series b', 'Slice Index', 'Chosen ROIs b',
                                                           'Radiation Range b', 'Median b', 'Mean b'])
            patient_full_data_df = pd.merge(patient_before_data_df, patient_after_data_df, on='Slice Index')
            patient_full_data_df = patient_full_data_df.drop('Chosen ROIs b', 1)
            patient_full_data_df = patient_full_data_df.drop('Radiation Range b', 1)
            patient_full_data_df['Median diff'] = abs(
                patient_full_data_df['Median a'] - patient_full_data_df['Median b'])
            patient_full_data_df['Mean diff'] = abs(
                patient_full_data_df['Mean a'] - patient_full_data_df['Mean b'])
            patient_full_data_df.to_csv('statistics_range'+'to'.join(str(e) for e in radrange)+'.csv', sep=';', encoding='utf-8')



        # print("Only ROI's - Second")
        # chosenROIs = ['Pluco (L)']

        # self.display_and_calculate_dicom_with_radiation_in_range_with_added_rois(self.segmented_before_images[50], 50, 0, 0, chosenROIs, SeriesType.BEFORE)
        # self.display_and_calculate_dicom_with_radiation_in_range_with_added_rois(self.segmented_before_images[48], 48, 0, 0, chosenROIs, SeriesType.BEFORE)
        # self.display_and_calculate_dicom_with_radiation_in_range_with_rois(self.elastic_registered_image, 33, 10, 30, chosenROIs)


        # DISPLAY AND CALCULATE VALUES IN SEGMENTED LUNGS IN CHOSEN RANGE AND INSIDE CHOSEN ROI's COMMON AREA (OBSZAR WSPÓLNY)
        # chosenCommonROIs = ['Pluco (P)', 'Dmin']
        # chosenCommonROIs = []
        # self.display_and_calculate_dicom_with_radiation_in_range_inside_common_part(self.patient.before_slices_image[33], 33, 0, 0, chosenCommonROIs)
        # self.display_and_calculate_dicom_with_radiation_in_range_inside_common_part(self.elastic_registered_image, 33, 10, 50, chosenCommonROIs)
        # self.display_and_calculate_dicom_with_radiation_in_range_inside_common_part(self.segmented_before_images[33], 33, 0, 0, chosenCommonROIs)

        # DISPLAY AND CALCULATE VALUES IN SEGMENTED LUNGS IN CHOSEN RANGE AND INSIDE CHOSEN ROI's SUBSTRACED AREA (PIERWSZY - DRUGI)
        # chosenCommonROIs = ['Pluco (P)', 'Dmin']
        # self.display_and_calculate_dicom_with_radiation_in_range_inside_subtracted_part(self.segmented_before_images[33], 33, 0, 0, chosenCommonROIs)



    def display_raw_dicom_slice(self, slice):
        pyplot.close('all')
        pyplot.axes().set_aspect('equal', 'datalim')

        pyplot.pcolormesh(self.patient.before_pixel_axes[0], self.patient.before_pixel_axes[1],
                          slice)

        pyplot.set_cmap('gray')
        pyplot.gca().invert_yaxis()

    def display_dicom_with_ROI(self, chosenROIs, slice_number, slice):

        self.display_raw_dicom_slice(slice)
        self.draw_chosen_rois(chosenROIs, slice_number)

        pyplot.show()

    def display_dicom_with_radiation_in_range_for_given_image(self, slice_number, low_level, high_level, slice):
        self.display_raw_dicom_slice(slice)

        pyplot.pcolormesh(self.patientRtDose.dose_pixel_axes[0], self.patientRtDose.dose_pixel_axes[1],
                          self.prepare_radiation_matrix(low_level, high_level, slice_number), alpha=0.2, vmin=0, vmax=1)
        pyplot.set_cmap('jet')

        pyplot.show()

    def display_dicom_with_full_radiation_as_heat_map(self, slice_number, slice):

        self.display_raw_dicom_slice(slice)

        pyplot.axes().set_aspect('equal', 'datalim')
        pyplot.pcolormesh(self.patientRtDose.dose_pixel_axes[0], self.patientRtDose.dose_pixel_axes[1],
                          self.patientRtDose.dose_array_gray[slice_number], alpha=0.2)

        pyplot.set_cmap('jet')
        # pyplot.set_cmap('gray')
        pyplot.show()

    def display_and_calculate_dicom_with_radiation_in_range_without_roi(self, slice, slice_number, min_rad, max_rad):
        self.display_raw_dicom_slice(slice)

        start_mask_calc = time.time()
        radiation_in_given_range_in_slice_matrix = self.prepare_radiation_matrix(min_rad, max_rad, slice_number)
        CT_radiation_mask_matrix = self.prepare_original_size_radiation_mask(radiation_in_given_range_in_slice_matrix, slice)

        calculated_values = self.calculate_statistic_values(CT_radiation_mask_matrix, slice)
        end_mask_calc = time.time()
        time_all_calc = start_mask_calc - end_mask_calc

        print("Mask + Calculation time: " + str(time_all_calc))

        CT_radiation_to_display = CT_radiation_mask_matrix.astype(float)
        CT_radiation_to_display[CT_radiation_to_display == 0] = np.NAN

        pyplot.pcolormesh(self.patient.before_pixel_axes[0], self.patient.before_pixel_axes[1],
                          # CT_radiation_mask_matrix, alpha=0.1, vmin=0, vmax=1, cmap='brg')
                          CT_radiation_to_display, alpha=0.1, cmap='brg')

        pyplot.show()

    def display_and_calculate_dicom_with_radiation_in_range_with_added_rois(self, slice, slice_number, min_rad, max_rad, chosenROIs, series_type):
        self.display_raw_dicom_slice(slice)
        self.draw_chosen_rois(chosenROIs, slice_number)

        start_mask_calc = time.time()

        original_size_radiation_mask, original_size_contour_mask = self.initialize_masks()

        if min_rad != 0 or max_rad != 0:
            radiation_in_given_range_in_slice_matrix = self.prepare_radiation_matrix(min_rad, max_rad, slice_number)
            original_size_radiation_mask = self.prepare_original_size_radiation_mask(radiation_in_given_range_in_slice_matrix, slice)

        if len(chosenROIs) > 0:
            self.slice_ROIs = self.patientStructureSet.all_slices_contours[slice_number]
            original_size_contour_mask = self.prepare_original_size_contour_mask(slice, chosenROIs)

        original_final_mask = numpy.logical_and(original_size_radiation_mask.astype(bool),
                                                original_size_contour_mask.astype(bool))

        calculated_values = self.calculate_statistic_values(original_final_mask, slice)
        calculated_values_raw = calculated_values.copy()
        values_to_csv = self.prepare_csv_values(series_type, slice_number, chosenROIs, min_rad, max_rad, calculated_values)

        end_mask_calc = time.time()
        time_all_calc = end_mask_calc - start_mask_calc

        print("Mask + Calculation time: " + str(time_all_calc))

        original_final_mask_to_display = original_final_mask.astype(float)
        original_final_mask_to_display[original_final_mask_to_display == 0] = np.NAN

        pyplot.pcolormesh(self.patient.before_pixel_axes[0], self.patient.before_pixel_axes[1],
                          original_final_mask_to_display, alpha=0.1, cmap='Wistia')

        pyplot.show()

        return values_to_csv, calculated_values_raw

    def display_and_calculate_dicom_with_radiation_in_range_inside_common_part(self, slice, slice_number, min_rad, max_rad, chosenROIs, series_type):
        self.display_raw_dicom_slice(slice)
        roi_exists = self.draw_chosen_rois(chosenROIs, slice_number)

        original_size_radiation_mask, original_size_contour_mask = self.initialize_masks()

        start_mask_calc = time.time()

        if roi_exists:
            if min_rad != 0 or max_rad != 0:
                radiation_in_given_range_in_slice_matrix = self.prepare_radiation_matrix(min_rad, max_rad, slice_number)
                original_size_radiation_mask = self.prepare_original_size_radiation_mask(radiation_in_given_range_in_slice_matrix, slice)

            if len(chosenROIs) > 0:
                if numpy.amax(original_size_radiation_mask) != 0:
                    original_size_contour_mask = self.prepare_original_size_common_contour_mask(slice, chosenROIs)
                else:
                    original_size_contour_mask.fill(0)
        else:
            original_size_radiation_mask.fill(0)
            original_size_contour_mask.fill(0)
            print("No chosen ROI found in this slice")

        original_final_mask = numpy.logical_and(original_size_radiation_mask.astype(bool),
                                                original_size_contour_mask.astype(bool))

        calculated_values = self.calculate_statistic_values(original_final_mask, slice)
        calculated_values_raw = calculated_values.copy()
        values_to_csv = self.prepare_csv_values(series_type, slice_number, chosenROIs, min_rad, max_rad,
                                                calculated_values)

        end_mask_calc = time.time()
        time_all_calc = end_mask_calc - start_mask_calc

        print("Mask + Calculation time: " + str(time_all_calc))

        original_final_mask_to_display = original_final_mask.astype(float)
        original_final_mask_to_display[original_final_mask_to_display == 0] = np.NAN

        #pyplot.pcolormesh(self.patient.before_pixel_axes[0], self.patient.before_pixel_axes[1],
                         # original_final_mask_to_display, alpha=0.1, cmap='Wistia')

        # self.calculate_statistic_values(original_final_mask, slice)

        #pyplot.show()

        return values_to_csv, calculated_values_raw

    def display_and_calculate_dicom_with_radiation_in_range_inside_subtracted_part(self, slice, slice_number, min_rad, max_rad, chosenROIs):
        self.display_raw_dicom_slice(slice)
        self.draw_chosen_rois(chosenROIs, slice_number)

        original_size_radiation_mask, original_size_contour_mask = self.initialize_masks()

        if min_rad != 0 or max_rad != 0:
            radiation_in_given_range_in_slice_matrix = self.prepare_radiation_matrix(min_rad, max_rad, slice_number)
            original_size_radiation_mask = self.prepare_original_size_radiation_mask(radiation_in_given_range_in_slice_matrix, slice)

        if len(chosenROIs) != 0:
            original_size_contour_mask = self.prepare_original_size_subtracted_contour_mask(slice, chosenROIs)

        original_final_mask = numpy.logical_and(original_size_radiation_mask.astype(bool),
                                                original_size_contour_mask.astype(bool))

        original_final_mask_to_display = original_final_mask.astype(float)
        original_final_mask_to_display[original_final_mask_to_display == 0] = np.NAN

        pyplot.pcolormesh(self.patient.before_pixel_axes[0], self.patient.before_pixel_axes[1],
                          # original_final_mask.astype(int), alpha=0.1, vmin=0, vmax=1, cmap='jet')
                          original_final_mask_to_display, alpha=0.1, cmap='Wistia')

        self.calculate_statistic_values(original_final_mask, slice)

        pyplot.show()

    def draw_chosen_rois(self, chosenROIs, slice_number):
        roi_exists = False

        self.slice_ROIs = self.patientStructureSet.all_slices_contours[slice_number]
        index = 0
        colors = [(1,0,0),(0,1,1)]
        for slice_ROI in self.slice_ROIs:
            for chosenROI in chosenROIs:
                if slice_ROI['referenceROIName'] == chosenROI:
                    if slice_ROI['contourPoints']:
                        roi_exists = True
                        roi_color = (slice_ROI['color'][0]/255, slice_ROI['color'][1]/255, slice_ROI['color'][2]/255)
                        for point in slice_ROI['contourPoints'][0]:
                            circle = pyplot.Circle((point[0], point[1]), 0.3, color=roi_color)
                            # circle = pyplot.Circle((point[0], point[1]), 0.3, color=colors[index])
                            pyplot.gcf().gca().add_artist(circle)
                    index = index + 1
        return roi_exists

    def prepare_radiation_matrix(self, min, max, slice_number):

        radiation_in_range_in_slice = numpy.empty(
            [self.patientRtDose.slice_dimensions[1],
             self.patientRtDose.slice_dimensions[0]],
            dtype='uint16')
        radiation_in_range_in_slice.fill(0)
        radiation_in_range_in_slice[(self.patientRtDose.dose_array_gray[slice_number] >= min) & (
                self.patientRtDose.dose_array_gray[slice_number] <= max)] = 1

        return radiation_in_range_in_slice

    def prepare_original_size_radiation_mask(self, radiation_in_range_in_slice, slice):

        CT_radiation_mask_matrix = numpy.empty(
            [self.patient.slice_dimensions[1], self.patient.slice_dimensions[0]],
            dtype='uint16')
        CT_radiation_mask_matrix.fill(0)
        if numpy.amax(radiation_in_range_in_slice) != 0:
            non_zero_radiation_indexes = numpy.nonzero(radiation_in_range_in_slice)
            upper_left_pixel_index = [non_zero_radiation_indexes[1].min(), non_zero_radiation_indexes[0].min()]
            lower_right_pixel_index = [non_zero_radiation_indexes[1].max(), non_zero_radiation_indexes[0].max()]

            upper_left_pixel_center_coords = self.find_rad_matrix_centers(upper_left_pixel_index[0], upper_left_pixel_index[1])
            lower_right_pixel_center_coords = self.find_rad_matrix_centers(lower_right_pixel_index[0], lower_right_pixel_index[1])

            idx_xmin = find_nearest_point_index(self.patient.before_pixel_centers_axes[0], upper_left_pixel_center_coords[0] - 1)
            idx_ymin = find_nearest_point_index(self.patient.before_pixel_centers_axes[1], upper_left_pixel_center_coords[1] - 1)
            idx_xmax = find_nearest_point_index(self.patient.before_pixel_centers_axes[0], lower_right_pixel_center_coords[0] + 1)
            idx_ymax = find_nearest_point_index(self.patient.before_pixel_centers_axes[1], lower_right_pixel_center_coords[1] + 1)

            half_x_spacing = self.patientRtDose.slice_pixel_spacing[0]/2
            half_y_spacing = self.patientRtDose.slice_pixel_spacing[1]/2

            centers_of_rad = []
            for i, _ in enumerate(non_zero_radiation_indexes[0]):
                centers_of_rad.append(self.find_rad_matrix_centers(non_zero_radiation_indexes[1][i], non_zero_radiation_indexes[0][i]))

            sub_arr = slice[idx_ymin:idx_ymax, idx_xmin:idx_xmax]
            analized_pixels = 0
            for index, x in numpy.ndenumerate(sub_arr):
                hu_coord = self.find_CT_matrix_centers(index[1] + idx_xmin, index[0] + idx_ymin)
                for rad_pixel_center in centers_of_rad:
                    if (rad_pixel_center[0] - half_x_spacing) < hu_coord[0] < (rad_pixel_center[0] + half_x_spacing) and (rad_pixel_center[1] - half_y_spacing) < hu_coord[1] < (rad_pixel_center[1] + half_y_spacing) and slice[index[0] + idx_ymin, index[1] + idx_xmin] > - 980:
                        CT_radiation_mask_matrix[index[0] + idx_ymin][index[1] + idx_xmin] = 1
                        break
                analized_pixels = analized_pixels +1

        return CT_radiation_mask_matrix

    def find_rad_matrix_centers(self, x_axis_idx, y_axis_idx):
        x_axis_coord = self.patientRtDose.dose_pixel_centers[0][x_axis_idx]
        y_axis_coord = self.patientRtDose.dose_pixel_centers[1][y_axis_idx]
        return x_axis_coord, y_axis_coord

    def find_CT_matrix_centers(self, x_axis_idx, y_axis_idx):
        x_axis_coord = self.patient.before_pixel_centers_axes[0][x_axis_idx]
        y_axis_coord = self.patient.before_pixel_centers_axes[1][y_axis_idx]
        return x_axis_coord, y_axis_coord

    def prepare_original_size_contour_mask(self, slice, chosenROIs):

        original_size_contour_mask = numpy.empty(
            [self.patient.slice_dimensions[1], self.patient.slice_dimensions[0]],
            dtype='uint16')
        original_size_contour_mask.fill(0)

        for slice_ROI in self.slice_ROIs:
            for chosenROI in chosenROIs:
                if slice_ROI['referenceROIName'] == chosenROI:
                    if slice_ROI['contourPoints']:

                        upper_left_contour_pixel_cords = [min(slice_ROI['contourPoints'][0], key=itemgetter(0))[0],
                                                     min(slice_ROI['contourPoints'][0], key=itemgetter(1))[1]]
                        lower_right_contour_pixel_cords = [max(slice_ROI['contourPoints'][0], key=itemgetter(0))[0],
                                                      max(slice_ROI['contourPoints'][0], key=itemgetter(1))[1]]

                        original_array_x_min_index = find_nearest_point_index(self.patient.before_pixel_centers_axes[0], upper_left_contour_pixel_cords[0] - 1)
                        original_array_y_min_index = find_nearest_point_index(self.patient.before_pixel_centers_axes[1], upper_left_contour_pixel_cords[1] - 1)
                        original_array_x_max_index = find_nearest_point_index(self.patient.before_pixel_centers_axes[0], lower_right_contour_pixel_cords[0] + 1)
                        original_array_y_max_index = find_nearest_point_index(self.patient.before_pixel_centers_axes[1], lower_right_contour_pixel_cords[1] + 1)

                        orignal_sub_arr_containing_whole_roi = slice[original_array_y_min_index:original_array_y_max_index, original_array_x_min_index:original_array_x_max_index]
                        analized_pixels = 0
                        for sub_arr_pixel_index, x in numpy.ndenumerate(orignal_sub_arr_containing_whole_roi):
                            # print(index)
                            original_pixel_centers_coords = self.find_CT_matrix_centers(sub_arr_pixel_index[1] + original_array_x_min_index, sub_arr_pixel_index[0] + original_array_y_min_index)
                            if inside_polygon_ray_tracing(original_pixel_centers_coords[0], original_pixel_centers_coords[1], slice_ROI['contourPoints'][0]) and slice[sub_arr_pixel_index[0] + original_array_y_min_index, sub_arr_pixel_index[1] + original_array_x_min_index] > -980:
                                original_size_contour_mask[sub_arr_pixel_index[0] + original_array_y_min_index][sub_arr_pixel_index[1] + original_array_x_min_index] = 1
                            analized_pixels = analized_pixels + 1
                            # print(str(analized_pixels) + "/" + str(orignal_sub_arr_containing_whole_roi.size))
                    break

        return original_size_contour_mask

    def prepare_original_size_common_contour_mask(self, slice, chosenROIs):

        original_size_common_contour_mask = numpy.empty(
            [self.patient.slice_dimensions[1], self.patient.slice_dimensions[0]],
            dtype='uint16')
        original_size_common_contour_mask.fill(0)

        for slice_ROI in self.slice_ROIs:
            for i, chosenROI in enumerate(chosenROIs):
                if slice_ROI['referenceROIName'] == chosenROI:
                    if slice_ROI['contourPoints']:

                        upper_left_contour_pixel_cords = [min(slice_ROI['contourPoints'][0], key=itemgetter(0))[0],
                                                     min(slice_ROI['contourPoints'][0], key=itemgetter(1))[1]]
                        lower_right_contour_pixel_cords = [max(slice_ROI['contourPoints'][0], key=itemgetter(0))[0],
                                                      max(slice_ROI['contourPoints'][0], key=itemgetter(1))[1]]

                        original_array_x_min_index = find_nearest_point_index(self.patient.before_pixel_centers_axes[0], upper_left_contour_pixel_cords[0] - 1)
                        original_array_y_min_index = find_nearest_point_index(self.patient.before_pixel_centers_axes[1], upper_left_contour_pixel_cords[1] - 1)
                        original_array_x_max_index = find_nearest_point_index(self.patient.before_pixel_centers_axes[0], lower_right_contour_pixel_cords[0] + 1)
                        original_array_y_max_index = find_nearest_point_index(self.patient.before_pixel_centers_axes[1], lower_right_contour_pixel_cords[1] + 1)

                        orignal_sub_arr_containing_whole_roi = slice[original_array_y_min_index:original_array_y_max_index, original_array_x_min_index:original_array_x_max_index]
                        analized_pixels = 0
                        for sub_arr_pixel_index, x in numpy.ndenumerate(orignal_sub_arr_containing_whole_roi):
                            original_pixel_centers_coords = self.find_CT_matrix_centers(sub_arr_pixel_index[1] + original_array_x_min_index, sub_arr_pixel_index[0] + original_array_y_min_index)
                            if inside_polygon_ray_tracing(original_pixel_centers_coords[0], original_pixel_centers_coords[1], slice_ROI['contourPoints'][0]) and slice[sub_arr_pixel_index[0] + original_array_y_min_index, sub_arr_pixel_index[1] + original_array_x_min_index] > -980:
                                    original_size_common_contour_mask[sub_arr_pixel_index[0] + original_array_y_min_index][sub_arr_pixel_index[1] + original_array_x_min_index] += 1

                            analized_pixels = analized_pixels + 1
                            print(str(analized_pixels) + "/" + str(orignal_sub_arr_containing_whole_roi.size))
                    break

        tmp_common_contour_mask = original_size_common_contour_mask == np.amax(original_size_common_contour_mask)
        return tmp_common_contour_mask.astype(np.uint16)

    def prepare_original_size_subtracted_contour_mask(self, slice, chosenROIs):

        original_size_subtracted_contour_mask = numpy.empty(
            [self.patient.slice_dimensions[1], self.patient.slice_dimensions[0]],
            dtype='uint16')
        original_size_subtracted_contour_mask.fill(0)

        for slice_ROI in self.slice_ROIs:
            for i, chosenROI in enumerate(chosenROIs):
                if slice_ROI['referenceROIName'] == chosenROI:
                    if slice_ROI['contourPoints']:

                        upper_left_contour_pixel_cords = [min(slice_ROI['contourPoints'][0], key=itemgetter(0))[0],
                                                     min(slice_ROI['contourPoints'][0], key=itemgetter(1))[1]]
                        lower_right_contour_pixel_cords = [max(slice_ROI['contourPoints'][0], key=itemgetter(0))[0],
                                                      max(slice_ROI['contourPoints'][0], key=itemgetter(1))[1]]

                        original_array_x_min_index = find_nearest_point_index(self.patient.before_pixel_centers_axes[0], upper_left_contour_pixel_cords[0] - 1)
                        original_array_y_min_index = find_nearest_point_index(self.patient.before_pixel_centers_axes[1], upper_left_contour_pixel_cords[1] - 1)
                        original_array_x_max_index = find_nearest_point_index(self.patient.before_pixel_centers_axes[0], lower_right_contour_pixel_cords[0] + 1)
                        original_array_y_max_index = find_nearest_point_index(self.patient.before_pixel_centers_axes[1], lower_right_contour_pixel_cords[1] + 1)

                        orignal_sub_arr_containing_whole_roi = slice[original_array_y_min_index:original_array_y_max_index, original_array_x_min_index:original_array_x_max_index]
                        analized_pixels = 0
                        for sub_arr_pixel_index, x in numpy.ndenumerate(orignal_sub_arr_containing_whole_roi):
                            original_pixel_centers_coords = self.find_CT_matrix_centers(sub_arr_pixel_index[1] + original_array_x_min_index, sub_arr_pixel_index[0] + original_array_y_min_index)
                            if inside_polygon_ray_tracing(original_pixel_centers_coords[0], original_pixel_centers_coords[1], slice_ROI['contourPoints'][0]) and i == 0 and slice[sub_arr_pixel_index[0] + original_array_y_min_index, sub_arr_pixel_index[1] + original_array_x_min_index] > -980:
                                    original_size_subtracted_contour_mask[sub_arr_pixel_index[0] + original_array_y_min_index][sub_arr_pixel_index[1] + original_array_x_min_index] = 1
                            elif inside_polygon_ray_tracing(original_pixel_centers_coords[0], original_pixel_centers_coords[1], slice_ROI['contourPoints'][0]) and i != 0 and original_size_subtracted_contour_mask[sub_arr_pixel_index[0] + original_array_y_min_index][sub_arr_pixel_index[1] + original_array_x_min_index] == 1:
                                original_size_subtracted_contour_mask[sub_arr_pixel_index[0] + original_array_y_min_index][
                                    sub_arr_pixel_index[1] + original_array_x_min_index] = 0

                            analized_pixels = analized_pixels + 1
                            print(str(analized_pixels) + "/" + str(orignal_sub_arr_containing_whole_roi.size))
                    break

        # tmp_common_contour_mask = original_size_subtracted_contour_mask == np.amax(original_size_subtracted_contour_mask)
        return original_size_subtracted_contour_mask.astype(np.uint16)

    def calculate_statistic_values(self, original_final_mask, slice):

        non_zero_arrays = numpy.nonzero(original_final_mask)

        try:
            non_zero_mask = slice[non_zero_arrays]

            #mininum = non_zero_mask.min()
            #maximum = non_zero_mask.max()
            non_zero_array_dicom_in_mask_result_list = non_zero_mask.tolist()
            median_value_in_HU = statistics.median(non_zero_array_dicom_in_mask_result_list)
            mean_value_in_HU = statistics.mean(non_zero_array_dicom_in_mask_result_list)
        except ValueError:
            #mininum = None
            #maximum = None
            median_value_in_HU = None
            mean_value_in_HU = None

        return [median_value_in_HU, mean_value_in_HU]

    def initialize_masks(self):
        original_size_radiation_mask = numpy.empty(
            [self.patient.slice_dimensions[1], self.patient.slice_dimensions[0]],
            dtype='uint16')
        original_size_common_contour_mask = numpy.empty(
            [self.patient.slice_dimensions[1], self.patient.slice_dimensions[0]],
            dtype='uint16')
        original_size_radiation_mask.fill(1)
        original_size_common_contour_mask.fill(1)

        return original_size_radiation_mask, original_size_common_contour_mask


    def prepare_csv_values(self, series_type, slice_number, chosenROIs, min_rad, max_rad, calculated_values):
        calculated_values.insert(0, series_type)
        calculated_values.insert(1, slice_number)
        if len(chosenROIs) > 0:
            calculated_values.insert(2, chosenROIs)
        else:
            calculated_values.insert(2, '-')
        calculated_values.insert(3, str(min_rad) + '-' + str(max_rad))
        return calculated_values

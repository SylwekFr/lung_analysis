import pydicom
import os
import numpy as np
import collections
import statistics
import operator
import SimpleITK as sitk
from viewer import *
import re

class PatientRtDose:

    def __init__(self, patient):
        self.begin_index = patient.before_start_index-1
        self.end_index = patient.before_end_index-1

        self.before_RT_Dose_data_files_paths = patient.before_RT_data_files_paths
        self.total_number_of_files = patient.total_before_files_count
        self.before_first_RT_Dose_file = self.get_rt_dose_file()

        self.dose_array = self.get_significant_slices()
        self.dose_array = self.dose_array[::-1]
        self.dose_array_gray = self.dose_array * self.before_first_RT_Dose_file.DoseGridScaling

        self.slice_dimensions = self.get_slice_dimensions()
        self.slice_pixel_spacing = self.get_slice_pixel_spacing()
        self.image_position_patient = self.get_image_position_patient()
        self.slice_thickness = 3
        self.number_of_slices = int(self.end_index-self.begin_index)

        self.dose_pixel_centers = self.get_pixels_centers_axes()
        self.dose_pixel_axes = self.get_pixels_axes()


    def get_rt_dose_file(self):
        for file in self.before_RT_Dose_data_files_paths:
            if re.search("RD", file):
                 return pydicom.read_file(file)

    def get_all_slices(self):
        return self.before_first_RT_Dose_file.pixel_array

    def get_significant_slices(self):
         return self.before_first_RT_Dose_file.pixel_array[self.total_number_of_files- self.end_index:self.total_number_of_files-self.begin_index, :, :]

    def get_slice_dimensions(self):
        return [int(self.before_first_RT_Dose_file.Columns), int(self.before_first_RT_Dose_file.Rows)]

    def get_slice_pixel_spacing(self):
        return [float(self.before_first_RT_Dose_file.PixelSpacing[0]), float(self.before_first_RT_Dose_file.PixelSpacing[1])]

    def get_image_position_patient(self):
        return [float(self.before_first_RT_Dose_file.ImagePositionPatient[0]), float(self.before_first_RT_Dose_file.ImagePositionPatient[1]), float(self.before_first_RT_Dose_file.ImagePositionPatient[2])]

    def get_pixels_centers_axes(self):
        x_center_axis_begin = self.image_position_patient[0]
        x_center_axis_end = x_center_axis_begin + (self.slice_dimensions[0])*self.slice_pixel_spacing[1]
        x_center_axis_step = self.slice_pixel_spacing[1]
        x_center_axis = np.arange(x_center_axis_begin, x_center_axis_end, x_center_axis_step)

        y_center_axis_begin = self.image_position_patient[1]
        y_center_axis_end = y_center_axis_begin + (self.slice_dimensions[1]) * self.slice_pixel_spacing[0]
        y_center_axis_step = self.slice_pixel_spacing[0]
        y_center_axis = np.arange(y_center_axis_begin, y_center_axis_end, y_center_axis_step)

        z_center_axis_begin = self.image_position_patient[2] + (self.total_number_of_files-1)*self.slice_thickness
        z_center_axis_begin = z_center_axis_begin - self.begin_index*self.slice_thickness
        z_center_axis_end = z_center_axis_begin - (self.number_of_slices) * self.slice_thickness
        z_center_axis_step = -self.slice_thickness
        z_center_axis = np.arange(z_center_axis_begin, z_center_axis_end, z_center_axis_step)

        return [x_center_axis, y_center_axis, z_center_axis]

    def get_pixels_axes(self):
        half_row_pixel_spacing = 0.5 * self.slice_pixel_spacing[0]
        half_column_pixel_spacing = 0.5 * self.slice_pixel_spacing[1]

        x_axis_begin = self.image_position_patient[0] - half_column_pixel_spacing
        x_axis_end = x_axis_begin + (self.slice_dimensions[0]) * self.slice_pixel_spacing[1] + half_column_pixel_spacing
        x_axis_step = self.slice_pixel_spacing[1]
        x_axis = np.arange(x_axis_begin, x_axis_end, x_axis_step)

        y_axis_begin = self.image_position_patient[1] - half_row_pixel_spacing
        y_axis_end = y_axis_begin + (self.slice_dimensions[1]) * self.slice_pixel_spacing[0] + half_row_pixel_spacing
        y_axis_step = self.slice_pixel_spacing[0]
        y_axis = np.arange(y_axis_begin, y_axis_end, y_axis_step)

        z_axis_begin = self.image_position_patient[2] + (self.total_number_of_files - 1) * self.slice_thickness + 0.5 * self.slice_thickness
        z_axis_begin = z_axis_begin - self.begin_index * self.slice_thickness
        z_axis_end = z_axis_begin - (
            self.number_of_slices) * self.slice_thickness - 0.5 * self.slice_thickness
        z_axis_step = -self.slice_thickness
        z_axis = np.arange(z_axis_begin, z_axis_end, z_axis_step)

        return [x_axis, y_axis, z_axis]






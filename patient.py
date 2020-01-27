import dicom
import os
import numpy as np
import collections
import statistics
import operator
import SimpleITK as sitk
from patientRtDose import PatientRtDose
from patientStructureSet import PatientStructureSet
from viewer import sitk_show


class Patient:

    def __init__(self, data_path, before_indexes, after_indexes):
        self.patient_id = 0

        self.before_start_index = before_indexes[0]
        self.before_end_index = before_indexes[1]

        self.after_start_index = after_indexes[0]
        self.after_end_index = after_indexes[1]

        self.before_slices_files_paths = []
        self.before_RT_data_files_paths = []
        self.after_slices_files_paths = []

        self.generate_file_lists(data_path)

        self.before_patient_data = dicom.read_file(self.before_slices_files_paths[0])
        self.after_patient_data = dicom.read_file(self.after_slices_files_paths[0])

        self.before_slices_image, self.before_slices = self.get_uncompressed_slices()
        self.after_slices_image, self.after_slices = self.get_compressed_slices(data_path+'after_images')

        self.total_before_files_count = len(self.before_slices_files_paths)

        self.select_significant_slices()

        self.patient_RT_dose = PatientRtDose(self)
        self.patient_structure_set = PatientStructureSet(self)
        self.number_of_slices = int(self.before_end_index - self.before_start_index)

        self.image_position_patient = self.get_image_position_patient()
        self.slice_dimensions = self.get_slice_dimensions()
        self.slice_pixel_spacing = self.get_slice_pixel_spacing()
        self.slice_thickness = 3
        self.begin_index = self.before_start_index - 1

        self.before_pixel_axes = self.get_pixels_axes()
        self.before_pixel_centers_axes = self.get_pixels_centers_axes()

    def get_image_position_patient(self):
        return [float(self.before_patient_data.ImagePositionPatient[0]),
                float(self.before_patient_data.ImagePositionPatient[1]),
                float(self.before_patient_data.ImagePositionPatient[2])]

    def get_slice_dimensions(self):
        return [int(self.before_patient_data.Columns), int(self.before_patient_data.Rows)]

    def get_slice_pixel_spacing(self):
        return [float(self.before_patient_data.PixelSpacing[0]), float(self.before_patient_data.PixelSpacing[1])]


    def generate_file_lists(self, data_path):
        for dirName, subdirList, fileList in os.walk(data_path):
            for subdir in subdirList:
                    for afterDirName, afterSubdirList, afterFileList in os.walk(os.path.join(data_path, subdir)):
                        for file in afterFileList:
                            if subdir.lower() == "after_images":
                                self.after_slices_files_paths.append(os.path.join(afterDirName, file))
                            if subdir.lower() == "before_data":
                                self.before_RT_data_files_paths.append(os.path.join(afterDirName, file))
                            if subdir.lower() == "before_images":
                                self.before_slices_files_paths.append(os.path.join(afterDirName, file))


    def get_uncompressed_slices(self):
        slices_dict = {}
        slices_full = [dicom.read_file(path) for path in self.before_slices_files_paths]
        slices_full.sort(key=lambda x: int(x.InstanceNumber))
        rescale_intercept = slices_full[0].RescaleIntercept

        try:
            slice_thickness = np.abs(slices_full[0].ImagePositionPatient[2] - slices_full[1].ImagePositionPatient[2])
        except:
            slice_thickness = np.abs(slices_full[0].SliceLocation - slices_full[1].SliceLocation)

        slices_pixel_arrays = np.stack([s.pixel_array.astype('int16') for s in slices_full])
        temp = np.empty_like(slices_pixel_arrays)
        temp.fill(rescale_intercept)
        slices_pixel_arrays_housefield = slices_pixel_arrays + temp
        for i, slice in enumerate(slices_full):
            slices_dict[slice.InstanceNumber-1] = {'ID': slice.SOPInstanceUID,
                                                        'SliceThickness': slice_thickness}

        return slices_pixel_arrays_housefield, slices_dict

    def get_compressed_slices(self, compressed_dir_path):
        slices_dict = {}
        slices_full = [dicom.read_file(path) for path in self.after_slices_files_paths]
        slices_full.sort(key=lambda x: int(x.InstanceNumber))

        try:
            slice_thickness = np.abs(slices_full[0].ImagePositionPatient[2] - slices_full[1].ImagePositionPatient[2])
        except:
            slice_thickness = np.abs(slices_full[0].SliceLocation - slices_full[1].SliceLocation)

        reader = sitk.ImageSeriesReader()
        filenamesDICOMFixed = reader.GetGDCMSeriesFileNames(compressed_dir_path)
        reader.SetFileNames(filenamesDICOMFixed)
        imgOriginal = reader.Execute()
        slices_pixel_arrays_housefield = sitk.GetArrayFromImage(imgOriginal).astype('int16')
        slices_pixel_arrays_housefield[slices_pixel_arrays_housefield < -1000] = -1000
        slices_pixel_arrays_housefield = slices_pixel_arrays_housefield[::-1]

        for i, slice in enumerate(slices_full):
            slices_dict[slice.InstanceNumber-1] = {'ID': slice.SOPInstanceUID,
                                                        'SliceThickness': slice_thickness}

        return slices_pixel_arrays_housefield, slices_dict

    def select_significant_slices(self):
        before_significant_array = self.before_slices_image[self.before_start_index-1:self.before_end_index-1, :, :]
        before_significant_indexes = list(range(self.before_start_index-1, self.before_end_index-1))

        after_reduced_array = self.after_slices_image[self.after_start_index-1:self.after_end_index-1, :, :]
        after_reduced_array_indexes = list(range(self.after_start_index-1, self.after_end_index-1))
        after_significant_indexes = []
        after_significant_array = []


        for i, slice in enumerate(after_reduced_array):
            if i % 3 == 0 or i % 22 == 0 or i == 1 or i == 181:
            # if i % 3 == 0:
                after_significant_array.append(slice)
                after_significant_indexes.append(after_reduced_array_indexes[i])

        after_significant_array = np.array(after_significant_array)

        self.before_slices_image = before_significant_array
        self.after_slices_image = after_significant_array

        before_slice_thickness = next(iter(self.before_slices.values()))['SliceThickness']

        significant_before_slices = {k: v for k, v in self.before_slices.items() if k in before_significant_indexes}
        significant_before_slices = {k-self.before_start_index+1: v for k, v in significant_before_slices.items() }
        significant_after_slices = {k: v for k, v in self.after_slices.items() if k in after_significant_indexes}
        # significant_after_slices = {k - self.after_start_index: v for k, v in significant_after_slices.items()}
        for i, item in significant_after_slices.items():
            item['SliceThickness'] = before_slice_thickness

        self.before_slices = significant_before_slices
        self.after_slices = significant_after_slices

        print("Before slices count" + str(len(self.before_slices_image)))
        print("After slices count" + str(len(self.after_slices_image)))
        if len(self.before_slices_image) == len(self.after_slices_image):
            print("EQUAL LENGTHS")
        else:
            return


        np.save("patient2\\significant_slices\\significantBeforeArray", before_significant_array)
        np.save("patient2\\significant_slices\\significantAfterArray", after_significant_array)

        for i, before_slice in enumerate(before_significant_array):
            before_slice_image = sitk.GetImageFromArray(before_slice)
            sitk_show(before_slice_image, title="patient2\\significant_slices\\before\\slice" + str(i))
            print("Before: " + str(i))

        for i, after_slice in enumerate(after_significant_array):
            after_slice_image = sitk.GetImageFromArray(after_slice)
            sitk_show(after_slice_image, title="patient2\\significant_slices\\after\\slice" + str(i))
            print("After: " + str(i))

    def get_pixels_centers_axes(self):
        x_center_axis_begin = self.image_position_patient[0]
        x_center_axis_end = x_center_axis_begin + (self.slice_dimensions[0]) * self.slice_pixel_spacing[1]
        x_center_axis_step = self.slice_pixel_spacing[1]
        x_center_axis = np.arange(x_center_axis_begin, x_center_axis_end, x_center_axis_step)

        y_center_axis_begin = self.image_position_patient[1]
        y_center_axis_end = y_center_axis_begin + (self.slice_dimensions[1]) * self.slice_pixel_spacing[0]
        y_center_axis_step = self.slice_pixel_spacing[0]
        y_center_axis = np.arange(y_center_axis_begin, y_center_axis_end, y_center_axis_step)

        z_center_axis_begin = self.image_position_patient[2] - (self.begin_index) * self.slice_thickness
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

        z_axis_begin = self.image_position_patient[2] + 0.5 * self.slice_thickness
        z_axis_begin = z_axis_begin - (self.begin_index) * self.slice_thickness
        z_axis_end = z_axis_begin - (self.number_of_slices) * self.slice_thickness - 0.5 * self.slice_thickness
        z_axis_step = -self.slice_thickness
        z_axis = np.arange(z_axis_begin, z_axis_end, z_axis_step)

        return [x_axis, y_axis, z_axis]




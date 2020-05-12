import os
import numpy as np
import SimpleITK as sitk
import matplotlib.pyplot as plt
import operator
import statistics
import collections
import pickle
import time

from viewer import *

class PatientAffineRegistration:

    def __init__(self, segmented_images_arrays):
        self.segmented_images_arrays = segmented_images_arrays
        self.before_segmented_images_arrays = segmented_images_arrays[0]
        self.after_segmented_images_arrays = segmented_images_arrays[1]
        self.multires_iterations = []
        self.metric_values = []

    def save_registration_steps(self, folder_path, registration_method, slice_number):
        global metric_values, multires_iterations
        metric_values.append(registration_method.GetMetricValue())

        plt.plot(metric_values, 'r')
        plt.xlabel('Iteration Number', fontsize=12)
        plt.ylabel('Metric Value', fontsize=12)
        plt.savefig(folder_path + "metric_values\\" + "metric_val_" + str(slice_number) + ".png")
        plt.close()

        global iteration_number
        print(iteration_number)

        plt.close('all')
        iteration_number += 1


    def start_plot(self):
        global metric_values, multires_iterations

        metric_values = []
        multires_iterations = []

        self.metric_values = []
        self.multires_iterations = []


    def end_plot(self):
        global metric_values, multires_iterations

        del metric_values
        del multires_iterations

        self.metric_values = None
        self.multires_iterations = None

        plt.close()

    def update_multires_iterations(self):
        self.multires_iterations.append(len(self.metric_values))

    def save_one_step_images(self, fixed, moving, transform, folder, slice_number):
        moving_transformed = sitk.Resample(moving, fixed, transform,
                                           sitk.sitkLinear, 0.0,
                                           moving.GetPixelIDValue())

        show_overlayed_images(moving_transformed, fixed, "affine_final_overlayed_slice_" + str(slice_number),
                              folder)

    def affine_registration(self):
        save_information = False


        folder_path = 'patient\\affine_registration_fixed\\'
        final_transforms = []
        moving_resampled_images = []

        global iteration_number
        iteration_number = 0

        for index, slices in enumerate(zip(self.before_segmented_images_arrays, self.after_segmented_images_arrays)):
            print("Registrating affine slice: " + str(index))

            before_image = sitk.GetImageFromArray(slices[0])
            after_image = sitk.GetImageFromArray(slices[1])

            before_slice_image = sitk.Cast(before_image, sitk.sitkFloat32)
            after_slice_image = sitk.Cast(after_image, sitk.sitkFloat32)


            initial_transform = sitk.CenteredTransformInitializer(before_slice_image,
                                                                  after_slice_image,
                                                                  sitk.AffineTransform(before_slice_image.GetDimension()),
                                                                  sitk.CenteredTransformInitializerFilter.GEOMETRY)

            moving_resampled = sitk.Resample(after_slice_image, before_slice_image, initial_transform, sitk.sitkLinear,
                                             0.0, after_slice_image.GetPixelID())

            show_overlayed_images(moving_resampled, before_slice_image, "overlayed_before_affine_registration_" + str(index),
                                  folder_path+"initial_transforms\\")

            registration_method = sitk.ImageRegistrationMethod()

            registration_method.SetMetricAsMattesMutualInformation(numberOfHistogramBins=50)
            registration_method.SetMetricSamplingStrategy(registration_method.REGULAR)
            registration_method.SetMetricSamplingPercentage(1)

            registration_method.SetInterpolator(sitk.sitkLinear)

            # registration_method.SetOptimizerAsGradientDescent(learningRate=2, numberOfIterations=50, convergenceMinimumValue=1e-4, convergenceWindowSize=5)
            # registration_method.SetOptimizerAsGradientDescent(learningRate=2, numberOfIterations=100, convergenceMinimumValue=1e-5, convergenceWindowSize=10)
            registration_method.SetOptimizerAsGradientDescent(learningRate=1, numberOfIterations=255)

            registration_method.SetOptimizerScalesFromPhysicalShift()

            registration_method.SetShrinkFactorsPerLevel(shrinkFactors = [6,4,2,1])
            registration_method.SetSmoothingSigmasPerLevel(smoothingSigmas=[3,2,1,0])
            registration_method.SmoothingSigmasAreSpecifiedInPhysicalUnitsOn()

            registration_method.SetInitialTransform(initial_transform)


            if(save_information):
                registration_method.AddCommand(sitk.sitkStartEvent, self.start_plot)
                registration_method.AddCommand(sitk.sitkEndEvent, self.end_plot)
                registration_method.AddCommand(sitk.sitkMultiResolutionIterationEvent, self.update_multires_iterations())


                registration_method.AddCommand(sitk.sitkIterationEvent, lambda: self.save_registration_steps(folder_path, registration_method, index))
                registration_method.AddCommand(sitk.sitkIterationEvent, lambda: self.save_one_step_images(before_slice_image,
                                                                                                        after_slice_image,
                                                                                                        initial_transform,
                                                                                                        folder_path+"transformed_images\\",
                                                                                                        index))
            # start_patient_affine_time = time.time()

            final_transform = registration_method.Execute(sitk.Cast(before_slice_image, sitk.sitkFloat32),
                                                          sitk.Cast(after_slice_image, sitk.sitkFloat32))
            # end_patient_affine_time = time.time()
            # patient_affine_time = end_patient_affine_time - start_patient_affine_time


            # print('Total affine 2D registration time: ' + str(patient_affine_time))
            # print('End metric value: {0}'.format(registration_method.GetMetricValue()))
            # print('Optimizer\'s stopping condition, {0}'.format(
            #     registration_method.GetOptimizerStopConditionDescription()))


            moving_resampled = sitk.Resample(after_slice_image, before_slice_image, final_transform, sitk.sitkLinear, 0.0, after_slice_image.GetPixelID())

            # np.save(folder_path+"moving_resampled_arrays\\moving_resampled_array_" + str(index), sitk.GetArrayFromImage(moving_resampled))
            # sitk.WriteTransform(final_transform, folder_path+"transformations\\" + "affine_transform_"+str(index)+".tfm")

            final_transforms.append(final_transform)
            moving_resampled_images.append(sitk.GetArrayFromImage(moving_resampled))

        return [final_transforms, moving_resampled_images]

import os
import numpy as np
import SimpleITK as sitk
import matplotlib.pyplot as plt
import operator
import statistics
import collections
import dicom
import pickle
import time

from viewer import *


class PatientElasticBSplineRegistration:

    def __init__(self, segmented_images_arrays, affine_information):
        self.before_segmented_images = segmented_images_arrays[0]
        self.after_segmented_images = segmented_images_arrays[1]
        self.affine_transformations = affine_information[0]
        self.affine_image_arrays = affine_information[1]
        self.multires_iterations = []
        self.metric_values = []

    def save_registration_steps(self, folder_path, registration_method, slice_number):

        print("learnign rate: "+ str(registration_method.GetOptimizerLearningRate()))
        self.metric_values.append(registration_method.GetMetricValue())

        plt.plot(self.metric_values, 'r')
        plt.plot(self.multires_iterations, [self.metric_values[index] for index in self.multires_iterations], 'b*')
        plt.xlabel('Iteration Number', fontsize=12)
        plt.ylabel('Metric Value', fontsize=12)
        plt.savefig(folder_path+"metric_values\\" + "metric_val_" + str(slice_number) + ".png")
        plt.close()

        global iteration_number
        print("Iteration: " + str(iteration_number))

        plt.close('all')
        iteration_number += 1

    def start_plot(self):
        self.metric_values = []
        self.multires_iterations = []

    def end_plot(self):
        self.metric_values = None
        self.multires_iterations = None
        plt.close()


    def update_multires_iterations(self):
        self.multires_iterations.append(len(self.metric_values))

    def save_one_step_images(self, fixed, moving, transform, folder, slice_number):
        moving_transformed = sitk.Resample(moving, fixed, transform,
                                           sitk.sitkLinear, 0.0,
                                           moving.GetPixelIDValue())

        show_overlayed_images(moving_transformed, fixed, "elastic_final_overlayed_slice_" + str(slice_number),
                              folder)


    def elastic_registration(self):

        save_information = False

        folder_path = 'patient\\elastic_bspline_registration_fixed\\'
        final_transforms = []
        moving_resampled_images = []
        global iteration_number
        iteration_number = 0

        for index, slice in enumerate(self.before_segmented_images):
            print("Registering elastic slice: " + str(index))

            before_slice_array = slice
            before_slice_image = sitk.GetImageFromArray(before_slice_array)

            after_affine_slice_image_array = self.affine_image_arrays[index]
            after_affine_slice_image_array[after_affine_slice_image_array == 0.0] = -1000
            after_slice_image = sitk.GetImageFromArray(after_affine_slice_image_array)

            before_slice_image = sitk.Cast( before_slice_image, sitk.sitkFloat32 )
            after_slice_image = sitk.Cast( after_slice_image, sitk.sitkFloat32 )



            transformDomainMeshSize=[8]*after_slice_image.GetDimension()
            initial_transform = sitk.BSplineTransformInitializer(before_slice_image,
                                                                  transformDomainMeshSize)

            moving_resampled = sitk.Resample(after_slice_image, before_slice_image, initial_transform, sitk.sitkLinear,
                                             0.0, after_slice_image.GetPixelID())

            show_overlayed_images(moving_resampled, before_slice_image,
                                  "overlayed_before_elastic_registration_" + str(index),
                                  folder_path + "initial_transforms\\")


            registration_method = sitk.ImageRegistrationMethod()

            registration_method.SetMetricAsMattesMutualInformation(numberOfHistogramBins=50)
            registration_method.SetMetricSamplingStrategy(registration_method.REGULAR)
            registration_method.SetMetricSamplingPercentage(1)

            registration_method.SetInterpolator(sitk.sitkLinear)

            # registration_method.SetOptimizerAsGradientDescent(learningRate=1.0, numberOfIterations=100, convergenceMinimumValue=1e-6, convergenceWindowSize=5)
            # registration_method.SetOptimizerAsGradientDescent(learningRate=4.0, numberOfIterations=100, convergenceMinimumValue=1e-4, convergenceWindowSize=5)
            registration_method.SetOptimizerAsGradientDescent(learningRate=1.0, numberOfIterations=255)

            registration_method.SetOptimizerScalesFromPhysicalShift()


            registration_method.SetShrinkFactorsPerLevel(shrinkFactors = [2,1])
            registration_method.SetSmoothingSigmasPerLevel(smoothingSigmas=[1,0])
            registration_method.SmoothingSigmasAreSpecifiedInPhysicalUnitsOn()


            registration_method.SetInitialTransform(initial_transform)

            if(save_information):
                registration_method.AddCommand(sitk.sitkStartEvent, self.start_plot)
                registration_method.AddCommand(sitk.sitkEndEvent, self.end_plot)
                registration_method.AddCommand(sitk.sitkMultiResolutionIterationEvent, self.update_multires_iterations())
                registration_method.AddCommand(sitk.sitkIterationEvent, lambda: self.save_registration_steps(folder_path, registration_method, index))
                registration_method.AddCommand(sitk.sitkIterationEvent,
                                               lambda: self.save_one_step_images(before_slice_image,
                                                                                 after_slice_image,
                                                                                 initial_transform,
                                                                                 folder_path + "transformed_images\\",
                                                                                 index))

            # start_patient_elastic_time = time.time()

            final_transform = registration_method.Execute(sitk.Cast(before_slice_image, sitk.sitkFloat32),
                                                          sitk.Cast(after_slice_image, sitk.sitkFloat32))
            # end_patient_elastic_time = time.time()
            # patient_elastic_time = end_patient_elastic_time - start_patient_elastic_time

            # print('Total elastic 2D registration time: ' + str(patient_elastic_time))
            # print('End metric value: {0}'.format(registration_method.GetMetricValue()))
            # print('Optimizer\'s stopping condition, {0}'.format(
            #     registration_method.GetOptimizerStopConditionDescription()))


            moving_resampled = sitk.Resample(after_slice_image, before_slice_image, final_transform, sitk.sitkLinear, 0.0, after_slice_image.GetPixelID())

            # np.save(folder_path+"moving_resampled_arrays\\" + "moving_elastic_array_" + str(index), sitk.GetArrayFromImage(moving_resampled))
            # sitk.WriteTransform(initial_transform, folder_path+"transformations\\bspline_transform_"+str(index)+".tfm")

            final_transforms.append(final_transform)
            moving_resampled_images.append(sitk.GetArrayFromImage(moving_resampled))

        return [final_transforms, moving_resampled_images]

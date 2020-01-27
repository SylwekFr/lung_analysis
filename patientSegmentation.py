import matplotlib.pyplot as plt
import numpy as np
from patient import Patient
import imreg_dft as ird

import SimpleITK as sitk
from viewer import *


class PatientSegmentation:


    def __init__(self, patient: Patient):
        self.patient = patient
        self.images_to_be_segmented = [patient.before_slices_image, patient.after_slices_image]
        self.final_segmented_images_arrays = []


    def full_segmentation(self):

        openingRadius = (3, 3, 3)
        closingRadius = (7, 7, 7)
        holeFillingRadius = (16, 16, 16)
        erosionFinalRadius = (1, 1, 1)
        kernel = sitk.sitkBall

        for j, image_to_be_segmented in enumerate(self.images_to_be_segmented):
            segmentedImages = []
            segmentedImagesArray = []

            print("J is: " + str(j))

            for i, slice in enumerate(image_to_be_segmented):

                print("Segmenting: " + str(i))

                image = sitk.GetImageFromArray(slice)
                size = image.GetSize()

                imgSmooth = sitk.CurvatureFlow(image1=image,
                                               timeStep=0.125,
                                               numberOfIterations=5)

                imgThresholded = imgSmooth > -300
                # imgThresholded = imgSmooth > -280

                imgNeighborhoodConnected = sitk.NeighborhoodConnected(imgThresholded, [(3,3,3)], upper = 0, lower = 0, radius =(0,0,0), replaceValue = 1)

                imgSmoothInt = sitk.Cast(sitk.RescaleIntensity(imgSmooth), imgNeighborhoodConnected.GetPixelID())

                imgNeighborhoodConnectedArray = sitk.GetArrayFromImage(imgNeighborhoodConnected)
                imgThresholdedArray = sitk.GetArrayFromImage(imgThresholded)

                resultArray = imgNeighborhoodConnectedArray + imgThresholdedArray
                resultArray[resultArray > 0] = 1

                lungs = sitk.GetImageFromArray(resultArray)
                lungs = 1 - lungs

                imgConnectedComponents = sitk.ConnectedComponent(lungs)

                lungs_found = []
                ratio = 0.0025
                image_size_in_pixels = np.shape(slice)[0] * np.shape(slice)[1]
                imgConnectedComponentsArray = sitk.GetArrayFromImage(imgConnectedComponents)
                imgConnectedComponentsArrayReduced = np.empty_like(imgConnectedComponentsArray)
                imgConnectedComponentsArrayReduced.fill(0)
                uniqueValues, occurCount = np.unique(imgConnectedComponentsArray, return_counts=True)
                uniqueValuesOccures = zip(uniqueValues, occurCount)
                uniqueValuesOccures = sorted(uniqueValuesOccures, key=lambda x: x[1], reverse=True)

                for idx, element in enumerate(uniqueValuesOccures):
                    if element[0]!= 0 and idx < 3:
                        # print("Lungs found: " + str(lungs_found))
                        # print("Element number: " + str(element[0]) + ", percentage: " + str(element[1] / image_size_in_pixels))
                        if element[1]/image_size_in_pixels > ratio and len(lungs_found)<2:
                            imgConnectedComponentsArrayReduced[imgConnectedComponentsArray == element[0]] = 1
                            lungs_found.append(element[0])

                imgConnectedComponentsReduced = sitk.GetImageFromArray(imgConnectedComponentsArrayReduced)

                if(len(lungs_found)==0):
                    print("No lungs in given CT image found")
                else:
                    erosionDilation = sitk.BinaryMorphologicalOpening(imgConnectedComponentsReduced, openingRadius, kernel)
                    morphologicalClosing = sitk.BinaryMorphologicalClosing(erosionDilation, closingRadius, kernel)
                    finalResult = sitk.VotingBinaryHoleFilling(morphologicalClosing, holeFillingRadius)
                    finalResult = sitk.BinaryErode(finalResult, erosionFinalRadius, kernel)

                    # sitk_show(sitk.LabelOverlay(imgSmoothInt, finalResult), title="image")

                    finalResultArray = sitk.GetArrayFromImage(finalResult)
                    oneOriginalSliceArray = sitk.GetArrayFromImage(image)
                    oneOriginalSliceArray[finalResultArray == 0] = -1000
                    segmentedImagesArray.append(oneOriginalSliceArray)
                    finalImage = sitk.GetImageFromArray(oneOriginalSliceArray)
                    segmentedImages.append(finalImage)

                    # print("Result on original image")
                    # sitk_show(finalImage)

            # if j==0:
            #     path_to_save = "patient\\segmented_slices\\before\\"
            #     np.save("patient\\segmented_slices\\segmentedBeforeArray", segmentedImagesArray)
            # else:
            #     path_to_save = "patient\\segmented_slices\\after\\"
            #     np.save("patient\\segmented_slices\\segmentedAfterArray", segmentedImagesArray)

            self.final_segmented_images_arrays.append(segmentedImagesArray)

            # for i, segmentedImage in enumerate(segmentedImages):
            #     sitk_show(segmentedImage, title=path_to_save+"segmented"+str(i))
            #     print("Saving: " + str(i))

        return self.final_segmented_images_arrays







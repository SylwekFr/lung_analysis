import numpy as np  # linear algebra
import dicom
import os
import scipy.ndimage as ndimage
import matplotlib.pyplot as plt
import SimpleITK as sitk
from skimage import measure, morphology, segmentation
from patient import Patient

class PatientWatershedSegmentation:

    def __init__(self, patient: Patient):
        self.patient = patient
        self.images_to_be_segmented = [patient.before_slices_image, patient.after_slices_image]
        self.final_segmented_images_arrays = []

    def create_internal_markers(self, image, slice_idx, total_slices):
        marker_internal = image < -360
        marker_internal = segmentation.clear_border(marker_internal)
        marker_internal_labels = measure.label(marker_internal)

        if slice_idx < 0.08 * total_slices:

            areas = [r.area for r in measure.regionprops(marker_internal_labels)]
            areas.sort()
            # ratio = 0.0015
            ratio = 0.00013
            if len(areas) > 2:
                regionprops = measure.regionprops(marker_internal_labels)
                for region in regionprops:
                    region_ratio = region.area / (np.size(marker_internal_labels))
                    if region_ratio < ratio:
                        for coordinates in region.coords:
                            marker_internal_labels[coordinates[0], coordinates[1]] = 0

            if (len(measure.regionprops(marker_internal_labels)) == 2):
                # print("Only two region props")

                region_props = [r for r in measure.regionprops(marker_internal_labels)]
                regions_min_y = []
                for region in region_props:
                    min_y = 512
                    for coordinates in region.coords:
                        min_y = coordinates[0] if coordinates[0] < min_y else min_y
                    regions_min_y.append({'label': region.label,
                                          'area': region.area,
                                          'min_y': min_y})
                regions_to_delete_labels = []
                regions_to_delete = []
                if regions_min_y[0]['area'] < regions_min_y[1]['area'] and regions_min_y[0]['min_y'] < regions_min_y[1][
                    'min_y']:
                    regions_to_delete_labels.append(regions_min_y[0]['label'])
                elif regions_min_y[1]['area'] < regions_min_y[0]['area'] and regions_min_y[1]['min_y'] < \
                        regions_min_y[0][
                            'min_y']:
                    regions_to_delete_labels.append(regions_min_y[1]['label'])

                for region_prop in region_props:
                    if region_prop.label in regions_to_delete_labels:
                        regions_to_delete.append(region_prop)

                for region in regions_to_delete:
                    for coordinates in region.coords:
                        marker_internal_labels[coordinates[0], coordinates[1]] = 0

            if (len(measure.regionprops(marker_internal_labels)) == 3):
                # print("Measuring algorithm")

                region_props = [r for r in measure.regionprops(marker_internal_labels)]
                region_centroids = []
                for region_prop in region_props:
                    region_centroids.append({'label': region_prop.label,
                                             'centroid': region_prop.centroid})
                res = list(zip(region_centroids, region_centroids[1:] + region_centroids[:1]))
                centroid_distances = []
                for pair in res:
                    distance = np.sqrt((pair[0]['centroid'][1] - pair[1]['centroid'][1]) ** 2 + (
                            pair[0]['centroid'][0] - pair[1]['centroid'][0]) ** 2)
                    centroid_distance = {'pair_labels': [pair[0]['label'], pair[1]['label']],
                                         'distance': distance}
                    centroid_distances.append(centroid_distance)

                sorted_distances = sorted(centroid_distances, key=lambda k: k['distance'])
                lung_labels = sorted_distances[-1]['pair_labels']
                region_props_to_delete = []
                for region_prop in region_props:
                    if region_prop.label not in lung_labels:
                        region_props_to_delete.append(region_prop)

                for region in region_props_to_delete:
                    for coordinates in region.coords:
                        marker_internal_labels[coordinates[0], coordinates[1]] = 0

        else:

            # ---------------------------------------------------------------------

            areas = [r.area for r in measure.regionprops(marker_internal_labels)]
            areas.sort()
            ratio = 0.0025
            if len(areas) > 2:
                regionprops = measure.regionprops(marker_internal_labels)
                for region in regionprops:
                    region_ratio = region.area / (np.size(marker_internal_labels))
                    if (region.area < areas[-1] and region.area < areas[-2]) or region_ratio < ratio:
                        for coordinates in region.coords:
                            marker_internal_labels[coordinates[0], coordinates[1]] = 0

        # ---------------------------------------------------------------------

        marker_internal = marker_internal_labels > 0
        return  marker_internal


    def create_external_markers(self, marker_internal):

        external_a = ndimage.binary_dilation(marker_internal, iterations=8)
        external_b = ndimage.binary_dilation(marker_internal, iterations=55)
        marker_external = external_b ^ external_a

        marker_watershed = np.zeros((512, 512), dtype=np.int)
        marker_watershed += marker_internal * 255
        marker_watershed += marker_external * 128

        return marker_external, marker_watershed


    def create_sobel_gradient(self, image):
        sobel_filtered_dx = ndimage.sobel(image, 1)
        sobel_filtered_dy = ndimage.sobel(image, 0)
        sobel_gradient = np.hypot(sobel_filtered_dx, sobel_filtered_dy)
        sobel_gradient *= 255.0 / np.max(sobel_gradient)

        return sobel_gradient

    def fill_border_holes_with_black_hat(self, outline):

        blackhat_struct = [[0, 0, 1, 1, 1, 0, 0],
                           [0, 1, 1, 1, 1, 1, 0],
                           [1, 1, 1, 1, 1, 1, 1],
                           [1, 1, 1, 1, 1, 1, 1],
                           [1, 1, 1, 1, 1, 1, 1],
                           [0, 1, 1, 1, 1, 1, 0],
                           [0, 0, 1, 1, 1, 0, 0]]

        blackhat_struct = ndimage.iterate_structure(blackhat_struct, 2)
        outline += ndimage.black_tophat(outline, structure=blackhat_struct)

        return outline



    def segment_lungs(self, image, slice_idx, total_slices):

        marker_internal = self.create_internal_markers(image, slice_idx, total_slices)
        marker_external, marker_watershed = self.create_external_markers(marker_internal)

        sobel_gradient = self.create_sobel_gradient(image)

        watershed = morphology.watershed(sobel_gradient, marker_watershed)

        outline = ndimage.morphological_gradient(watershed, size=(3, 3))
        outline = outline.astype(bool)

        outline = self.fill_border_holes_with_black_hat(outline)
        lungfilter = np.bitwise_or(marker_internal, outline)

        lungfilter = ndimage.morphology.binary_opening(lungfilter, structure=np.ones((2, 2)), iterations=3)
        lungfilter = ndimage.morphology.binary_closing(lungfilter, structure=np.ones((4, 4)), iterations=5)
        lungfilter = ndimage.binary_erosion(lungfilter, structure=np.ones((3, 3)), iterations=2)

        segmented = np.where(lungfilter == 1, image, -1000 * np.ones((512, 512)))

        return segmented, lungfilter, outline, watershed, sobel_gradient, marker_internal, marker_external, marker_watershed

    def watershed_segmentation(self):
        segmented_images = []
        total_slices = len(self.images_to_be_segmented[0])

        for j, image_to_be_segmented in enumerate(self.images_to_be_segmented):

            segmented_slices = []

            # print("J is: " + str(j))

            for i, slice in enumerate(image_to_be_segmented):

                # print("Segmenting: " + str(i))

                segmented, lungfilter, outline, watershed, sobel_gradient, marker_internal, marker_external, marker_watershed = self.segment_lungs(slice, i, total_slices)

                segmented_slices.append(segmented)

                # print("Original Slice")
                # plt.imshow(slice, cmap='gray')
                # plt.show()
                # print("Internal Marker")
                # plt.imshow(marker_internal, cmap='gray')
                # plt.show()
                # print("External Marker")
                # plt.imshow(marker_external, cmap='gray')
                # plt.show()
                # print("Watershed Marker")
                # plt.imshow(marker_watershed, cmap='gray')
                # plt.show()
                # print ("Sobel Gradient")
                # plt.imshow(sobel_gradient, cmap='gray')
                # plt.show()
                # print ("Watershed Image")
                # plt.imshow(watershed, cmap='gray')
                # plt.show()
                # print ("Outline after reinclusion")
                # plt.imshow(outline, cmap='gray')
                # plt.show()
                # print ("Lungfilter after closing")
                # plt.imshow(lungfilter, cmap='gray')
                # plt.show()
                # print ("Segmented Lung")
                # plt.imshow(segmented, cmap='gray')
                # plt.show()

            # if j == 0:
            #     path_to_save = "patient\\watershed_segmented_slices_fixed\\before\\"
            #     np.save("patient\\watershed_segmented_slices_fixed\\segmentedBeforeArray", segmented_slices)
            # else:
            #     path_to_save = "patient\\watershed_segmented_slices_fixed\\after\\"
            #     np.save("patient\\watershed_segmented_slices_fixed\\segmentedAfterArray", segmented_slices)

            segmented_images.append(segmented_slices)

            # for i, segmented_slice in enumerate(segmented_slices):
            #     print("Saving image")
            #     plt.imsave(path_to_save+"segmented"+str(i), segmented_slice, cmap="gray")


        return segmented_images
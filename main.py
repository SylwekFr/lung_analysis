import matplotlib.pyplot as plt
import numpy as np
from patient import Patient
from patientSegmentation import PatientSegmentation
from patientAffineRegistration import PatientAffineRegistration
from patientRtDose import PatientRtDose
from patientStructureSet import PatientStructureSet
from patientWatershedSegmentation import PatientWatershedSegmentation
from preprocessedImagesViewer import PreprocessedImagesViewer
from patientElasticBSplineRegistration import PatientElasticBSplineRegistration
import imreg_dft as ird
import time
import pydicom

import SimpleITK as sitk
from viewer import *

data = 'Data\\'

example = pydicom.dcmread("example.dcm", force=True)

print("Start patient loading")
patient = Patient(data, [33, 102], [9, 191]) #od którego pliku z kolei licząc od 1, do którego z kolei
print("End patient loading")


print("Start SITK segmentation")
patient_segmentation = PatientSegmentation(patient)
segmentation_start = time.time()
segmented_images = patient_segmentation.full_segmentation()
segmentation_end = time.time()
segmentation_time = segmentation_end - segmentation_start
print("Segmentation time: " + str(segmentation_time))
print("End SITK segmentation")


#print("Start watershed segmentation")
#patient_watershed_segmentation = PatientWatershedSegmentation(patient)
#segmentation_start = time.time()
#segmented_images = patient_watershed_segmentation.watershed_segmentation()
#segmentation_end = time.time()
#segmentation_time = segmentation_end - segmentation_start
#print("Segmentation time: " + str(segmentation_time))
#print("End watershedsegmentation")


# segmented_images = [np.load("patient\\segmented_slices\\segmentedBeforeArray.npy"),
#                     np.load("patient\\segmented_slices\\segmentedAfterArray.npy")]

# segmented_images = [np.load("patient\\watershed_segmented_slices_fixed\\segmentedBeforeArray.npy"),
#                     np.load("patient\\watershed_segmented_slices_fixed\\segmentedAfterArray.npy")]
#
print("Start affine registration")
affine_start = time.time()
patient_affine_registration = PatientAffineRegistration(segmented_images)
affine_registration_output = patient_affine_registration.affine_registration()
affine_end = time.time()
print("Affine time: " + str(affine_end-affine_start))
print("End affine registration")


# print("Start elastic registration")
# elastic_start = time.time()
# patient_elastic_registration = PatientElasticBSplineRegistration(segmented_images, affine_registration_output)
# elastic_registration_output = patient_elastic_registration.elastic_registration()
# elastic_end = time.time()
# print("Elastic time: " + str(elastic_end-elastic_start))
# print("End elastic registration")

elastic_registration_output = [0,1]

preprocessed_images_viewer = PreprocessedImagesViewer(segmented_images, elastic_registration_output[1], patient)




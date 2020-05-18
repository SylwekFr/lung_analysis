#Lung analysis system  
## aims & scope  
This solution developed in Python have for objective to give an analysis of 
the effect of radiotherapy on the lungs. For completed this aim the software compare the Hounsfield Unit (HU) of the
medical images before and after the radiotherapy. The HU is " a relative quantitative measurement of radio density used by radiologists in the interpretation of computed
tomography (CT) images" [_Tami D. DenOtter; Johanna Schubert](https://www.ncbi.nlm.nih.gov/books/NBK547721/), to
simplify we can compare it to the grayscale measure.   
This analysis is based on CT scan of a patient before Radiation therapy and after radiation therapy in
[DICOM](https://www.dicomstandard.org/) format (standard format medical image and information).  
The first step is to process the images of the CT scan for the CT before and after the radiotherapy, a CT give several 
image as it virtually "cut" the lung into slices and take a picture    
Once the image are processed we need to make the slice of the after serial and before serial comparable, in another words
having the same number of slice in both serial, but also we need that for each comparison both slice have the same lungs
area, else the comparison would make no sense. So the software limit the number of slices of the bigger serial to be
the same than in the smaller one by taking every n slice.  
Then we can go to the calculation part, the software calculate the mean and median HU unit value for each slide of both
serial and compare those values for both before and after corresponding slices, all of this is done in data frame.
Finally the final dataframe in exported in .csv.
##requirements  
* Python V3.7 (if you don't yet have python, I recommend to install through [Anaconda](https://www.anaconda.com/) for 
automatic installation of some libraries)  
* Libraries:  
  * [time](https://docs.python.org/3/library/time.html)  
  * [pydicom](https://pydicom.github.io/)  
  * [os](https://docs.python.org/3/library/os.html)  
  * [math](https://docs.python.org/3/library/math.html)  
  * [numpy](https://numpy.org/)  
  * [SimpleITK](https://pypi.org/project/SimpleITK/)  
  * [scipy](https://www.scipy.org/)  
  * [skimage](https://scikit-image.org/)  
  * [aenum](https://pypi.org/project/aenum/)  
  * [operator](https://docs.python.org/3/library/operator.html)  
  * [statistics](https://docs.python.org/3/library/statistics.html)  
  * [pandas](https://pandas.pydata.org/)  
  * [matplotlib](https://matplotlib.org/)
  * [tqdm](https://github.com/tqdm/tqdm)  
##Recommended development environment  
We recommend the usage of PyCharm for this project, the free community edition is sufficient. However nothing forbid you
to use another environment such as for example Pyzo or Visual Studio Code.  
The solution was developed on Windows but Python is also running well on Apple and Linux, however some change might be 
required on those operating systems.  
###Computer resources  
* Minimal:
  * Processor I5
  * 8gb of Ram
  * SSD disk
* Recommended:
  * Processor i7
  * 32gb of Ram
  * SSD disk
##work  
- [x] image segmentation
- [x] image registration
- [x] HU calculation / Roi in slice
- [x] Before/After difference calculation
- [ ] Remove ROI from ROI
- [x] Automation of calculation for several range
- [ ] Automation of calculation for several patient
- [x] Data exportation to CSV 
- [x] Performance analysis 
- [ ] Refactoring  
##How to use :  
In the root in the folder Data make 3 folders: "before_data", "before-images" and "after_images". in the folder
before data put RS, RP and RD files, in before_images the iRT serial and in after_data the follow up serial
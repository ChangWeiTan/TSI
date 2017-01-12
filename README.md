# TSI
TSI: Indexing and classifying gigabytes of time series under time warping.

# Preamble
This is the source code for the work published in SDM 2017: Indexing and classifying gigabytes of time series under time warping

Authors: Chang Wei Tan, Geoffrey I. Webb, Francois Petitjean

When using this repository, please cite:
```
@INPROCEEDINGS{Tan2017-SDM,
  author = {Tan, Chang Wei and Webb, Geoffrey I. and Petitjean, Francois},
  title = {Indexing and classifying gigabytes of time series under time warping},
  booktitle = {SIAM International Conference on Data Mining},
  year = 2017,
  pages = {1--10}
}
```

If you want to use the code or found any bug in the code, please drop me an email at chang.tan@monash.edu. Thanks!

# Executing the code
Before using the code, ensure that these files and folders exist in your project directory which can be obtained from https://cloudstor.aarnet.edu.au/plus/index.php/s/pRLVtQyNhxDdCoM and https://drive.google.com/open?id=0B8Cg6Izm3IJybWxnWDJPeWZQWVk 
  1. Output folder: outputs/L experiment exists in your project directory
  2. SITS folder: dataset/SITS_2006_NDVI_C/SITS1M_fold*fold number* (e.g.dataset/SITS_2006_NDVI_C/SITS1M_fold1)
  3. UCR folder: dataset/UCR_Time_Series_Archive/*UCR datasets name*
  4. Also, make sure that in the dataset folders, you have a csv file for the properties of the dataset. The file should have the   following format: "nb of class,size of training set,size of testing set,time series length,warping window". Check the UCR Time Series website, http://www.cs.ucr.edu/~eamonn/time_series_data/ for the properties of each dataset. For SITS, the properties are "24,900000,100000,46,4". Feel free to change these parameters to suit your program. 
  5. Exact indices of 1NN-DTW using best warping window: index1NN/*SITS or UCR*/*SITS1M_fold#_1NN_LB_index1NN.csv or UCRDataset_1NN_LB_index1NN.csv*. The indices are sorted in the order downloaded from the UCR Time Series website. 

The main files to run the program are
  - SITS_NNDTW.java
  - SITS_NNED.java
  - SITS_TSI.java
  - UCR_NNDTW.java
  - UCR_TSI.java

## Example 1, Running TSI on UCR dataset 50words. 
Run UCR_TSI.java or SITS_TSI.java

Assuming using Eclipse, go to Run Configurations > Arguments 

The program arguments are 
  1. Project path                                             : *project directory*/src
  2. Dataset name                                             : 50words / change to SITS1M_fold# for SITS 
  3. Branching factor                                         : 3
  4. Max k-means iterations                                   : 10
  5. Number of NN                                             : 1
  6. Test number                                              : 10
  7. Candidate intervals to record results                    : 10
  8. Number of results to record before seeing 1st candidate  : 5
  9. Number of time series to examine per query               : 10
  10. Warping window size (in terms of length)                : 10

Alternatively, can change these parameters individually in the code

## Example 2, Running NNDTW on UCR dataset 50words. 
Run UCR_NNDTW.java, SITS_NNDTW.java or SITS_NNED.java (if NNED, there will be no warping window)

Assuming using Eclipse, go to Run Configurations > Arguments 

The program arguments are 
  1. Project path                                             : *project directory*/src
  2. Dataset name                                             : 50words / change to SITS1M_fold# for SITS 
  3. Test number                                              : 10
  4. Candidate intervals to record results                    : 10
  5. Warping window size (in terms of length)                 : 10

Alternatively, can change these parameters individually in the code

The results will be stored in outputs/L experiment/*dataset name*_*experiment*.csv

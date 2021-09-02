# LOCAL TUNING BIASES IN THE MOUSE dLGN

<img src="lgnASSETS/Three mechanisms of Tuning Biases.jpg" alt="Screen Shot 2021-09-02 at 2.56.34 PM" style="zoom: 50%;" />

## Overview

Local tuning biases is a hypothesized property of populations of visual neurons in the early visual system. Specifically, it is the notion that groups neurons/units process a given image region with a similar, or biased, set of tuned filters.  

The aim of this project is to assess the presence of these 'local tuning biases' among populations of LGN afferents innervating  mouse V1  (the thalamic axons that relay retinal signals to the cortex) . 

To answer this question, we estimated the spatial receptive fields and feature tuning of populations of thalamic boutons innervating mouse V1. Then, for each pair of boutons, we determined the extent to which their receptive fields overlapped and the similarity of their tuning profiles. If LGN boutons with overlapping receptive fields show similar tuning profiles, this would suggest the presence of local tuning biases.

<u>This Repository</u>

- The remainded of this README file provides a breif description of the methods that were carried out to assess this prediction. 
- The complete PDF manuscript for this project can be [ found here ](lgnSUBMIT/lgnManuscript_master.pdf) 
- The data used for this project can be found in [this folder](lgnDATA)
- The code for the analysis was written in MATLAB, and can be found in  [this folder](lgnANALYSIS)
  - [This MATLAB script ]( lgnSUBMIT/lgnManuscript_master.pdf) will guide you through the analysis 

## Methods

<u>Data collection</u> 

First we measure the visual properties of mouse V1 boutons using resonant, two-photon microscopy in the awake, behaving mouse. A total of three mice were used in this study. We imaged 3 unique cortical regions to obtain the data discussed here.

<u>Kernel Estimation</u> 

A bouton’s joint tuning to orientation and spatial frequency was measured by presenting a sequence of flashed, high-contrast sinusoidal gratings having pseudorandom orientations and spatial frequencies [FIGURE]. We estimated the tuning of each bouton by linearly regressing its response on the grating stimulus and denote the estimated tuning kernel of the ith bouton in an imaging field. The peak of the tuning kernel also yielded a spatial frequency preference and an orientation preference for each bouton.

We also measured the spatial receptive field maps of boutons by presenting a sequence of flickering, elongated bars at random orientations and positions across the visual field [FIGURE]. A boutons RF map was estimated by correlating its response and the location of the bars. 

<u>Image Processing</u>

After imaging a field, we manually segmented regions of interest  (ROIS) which corresponded micrometer-sized circular or elliptical boutons of dLGN axons. Then for each ROI we extracted signals by computing the mean of the calcium fluorescence within each region of interest and discounting the signals from the nearby neuropil. Spikes were then estimated via deconvolution.

<u>Analysis</u>

The primary goal is to determine if LGN boutons with overlapping receptive fields have similar tuning profiles. This would suggest  that neurons/units in the dLGN 'process' image regions with a similar, or biased, set of orientated filters. 

To do test this we take each possible pair of boutons (i,j)  in an imaging field, then compute the overlap between their receptive fields, and also compute the similarity of their tuning profiles. We can define the tuning similarity of bouton pairs as the as the correlation between their joint tuning kernels. The RF overlap of bouton pairs can also be defined as the correlation coefficient between their RF maps.

Once we compute the tuning similarty and receptive field overlap of bouton pairs, we can plot the relationship between the two variables.  The plot below shows that there exists a relationship between receptive field overlap (x-axis) and tuning similarity (yaxis). Neurons with overlapping receptive fields show similar tuning profiles. T

 





![Results_raw](lgnASSETS/Results_raw.png)


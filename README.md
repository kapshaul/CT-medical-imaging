# Medical Imaging Shot Noise Removal

**Authors:** Yong-Hwan Lee and Tony Storey

## Overview

This project simulates the Expectation Maximization (EM) algorithm for removing Poisson noise from medical images. The EM algorithm is particularly useful in this context because it allows for the estimation of the underlying image by iteratively maximizing the likelihood function, which accounts for the statistical nature of Poisson noise. Given the quantum nature of particles and their discrete arrival times, Poisson noise often manifests in medical imaging, especially in modalities like X-ray Computed Tomography (CT). This noise can lead to significant artifacts in the reconstructed images, which might result in diagnostic inaccuracies.

<br>

<div align="center">
    
<img src="https://github.com/kapshaul/ct-medical-imaging/blob/master/images/CT%20scan.jpg" width="400">&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;
<img src="https://github.com/kapshaul/ct-medical-imaging/blob/master/images/Computed_tomography_of_human_brain.png" width="400">

**Figure 1**: Left - Illustration of a CT scan procedure; Right - CT scan images

</div>

<br>

To mitigate these issues, the EM algorithm is applied as it effectively separates the noise from the actual signal in the image data. By modeling the image acquisition process and noise characteristics, the algorithm iteratively refines the image estimate, ultimately converging on a solution that minimizes the impact of noise.

Find more details from the report: [PDF](https://github.com/kapshaul/ct-medical-imaging/blob/master/Medical_Imaging.pdf)

## Problem Formulation

### 1. Observation

The distribution of $N_{ij}$ is given by,

$$
N_{ij} \sim Pois(a_{ij} \lambda_j)
$$

Here, $Pois$ denotes the Poisson distribution with parameter $\lambda$. Consider a model matrix $A$, where $A = (a_{ij})$ $i = 1, \dots, n$ $j = 1, \dots, m$.

Then, the observations $Y_{i=1,...n}$ can be written as below,

$$
Y_i = \sum_{j=1}^m N_{ij} \sim \sum_{j=1}^m Pois((a_{ij} \lambda_j)
$$

### 2. 3x3 cross section of voxel model

Each pixel represents the absorption coefficient. Below are examples of voxel models,

```
                                |        *--------------------*       |  
                                | y3---\ |      |      |      |       | 
                                |   ---/ |  p1  |  p2  |  p3  |       | 
                                -        *--------------------*       - 
                                | y2---\ |      |      |      |       | 
                                |   ---/ |  p4  |  p5  |  p6  |       | 
                                -        *--------------------*       - 
                                | y1---\ |      |      |      |       | 
                                |   ---/ |  p7  |  p8  |  p9  |       | 
                                |        *--------------------*       |  
```

```                            
                                         -------|------|-------
                                            y9     y10    y11
                                            ||     ||     ||
                                            \/     \/     \/
                                         *--------------------*          
                                         |      |      |      |        
                                         |  p1  |  p2  |  p3  |        
                                         *--------------------*        
                                         |      |      |      |        
                                         |  p4  |  p5  |  p6  |        
                                         *--------------------*        
                                         |      |      |      |        
                                         |  p7  |  p8  |  p9  |        
                                         *--------------------*         
                               
                                         -------|------|------- 
```

## EM Algorithm

### 1. Likelihood function

For the observations, the likelihood function can be written as,

$$
L(N_{ij})_ {ij}(\lambda) = \prod_i^n\prod_j^m e^{a_{ij} \lambda_j} \frac{(\lambda_j a_{ij})^{N_{ij}}}{N_{ij}!}
$$

Then, log-likelihood function can be below,

$$
l(N_{ij})_ {ij}(\lambda) = \sum_i^n \sum_j^m (-\lambda_j a_{ij} + N_{ij}\log{(\lambda_j a_{ij})} -\log{(N_{ij}!)})
$$

## Implementation

To implement the code, follow these steps:
1. Clone the repository, which includes the `Med5.m` file.
2. Run the `Med5.m` file to complete the EM algorithm and estimate the body model matrix coefficients.
   
The implementation will also plot Mean Squared Error (MSE) graphs, comparing the Cramer-Rao Lower Bound (CRLB) under the fixed body model matrix.

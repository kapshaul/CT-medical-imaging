# Medical Imaging Shot Noise Removal

**Authors:** Yong-Hwan Lee & Tony Storey

## Overview

This project simulates the Expectation Maximization (EM) algorithm for removing Poisson noise from medical images. The EM algorithm is particularly useful in this context because it allows for the estimation of the underlying image by iteratively maximizing the likelihood function, which accounts for the statistical nature of Poisson noise. Given the quantum nature of particles and their discrete arrival times, Poisson noise often manifests in medical imaging, especially in modalities like X-ray Computed Tomography (CT). This noise can lead to significant artifacts in the reconstructed images, which might result in diagnostic inaccuracies.

To mitigate these issues, the EM algorithm is applied as it effectively separates the noise from the actual signal in the image data. By modeling the image acquisition process and noise characteristics, the algorithm iteratively refines the image estimate, ultimately converging on a solution that minimizes the impact of noise.

Find more details from the report: [PDF](https://github.com/kapshaul/ct-medical-imaging/blob/master/Medical_Imaging.pdf)

### 1. The observation for each particle ray after poisson noise:

$$
y_n = \text{Poisson}((A \theta)_n)
$$

### 2. The observation for the poisson probability distribution:


$$
Poisson(Y_1, Y_2, Y_3, \ldots, Y_m | \theta_1, \theta_2, \theta_3, \ldots, \theta_n) = \frac{(A\theta)^Y e^{-(A\theta)}}{Y!}
$$

where:
- $Y$ represents the observation matrix.
- $\theta$ is the true coefficient that we aim to estimate.
- $A \in \mathbb{R}^{i \times j}$ is the body model matrix.

### 3. 3x3 cross section of voxel model example
Each pixel models the absorption coefficient.

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

## Implementation

To implement the code, follow these steps:
1. Clone the repository, which includes the `Med5.m` file.
2. Run the `Med5.m` file to complete the EM algorithm and estimate the body model matrix coefficients.
   
The implementation will also plot Mean Squared Error (MSE) graphs, comparing the Cramer-Rao Lower Bound (CRLB) under the fixed body model matrix.

# Medical Imaging Shot Noise Removal

**Authors:** YongHwan Lee & Tony Storey

## Overview

This project focuses on removing Poisson noise from medical images using the Expectation Maximization (EM) algorithm. The EM algorithm is applied for poisson noise removal.

Find more details from the report: [PDF](https://github.com/neurokimchi/t-medical-imaging/blob/master/Medical_Imaging.pdf)

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

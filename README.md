# Medical Imaging Shot Noise Removal

### Poisson Noise Removal via Expectation Maximization

**Authors:** YongHwan Lee & Tony Storey

---

### Overview

This project focuses on removing Poisson noise from medical images using the Expectation Maximization (EM) algorithm.

### Given

The following equations and concepts are used in this project:

- Probability Distribution:
  
$$
P(\Lambda) = \frac{\exp(A_{ij}\theta_i) \cdot (A_{ij}\theta_i)^{Y_{ij}}}{Y_{ij}!} \quad \quad where A \in \mathbb{R}^{i \times j}
$$

- Poisson Noise:

$$
y_k = \text{Poisson}((A \theta)_k) \quad \quad where k = 1,2,3,\ldots,n
$$

The EM algorithm is applied for noise removal.

### Voxel Model

#### 3x3 Cross Section of Voxel Model
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

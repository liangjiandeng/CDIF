# A New Context-Aware Details Injection Fidelity with Adaptive Coefficients Estimation for Variational Pansharpening [Matlab Code]

**Homepage:** https://liangjiandeng.github.io/
 
# Citation
```bibtex
@ARTICLE{xiao2022tgrs,
author={J.-L. Xiao, T.-Z. Huang, L.-J. Deng, Z.-C. Wu, and G. Vivone},
booktitle={IEEE Trans. Geosci. Remote Sens.},
title={A New Context-Aware Details Injection Fidelity with Adaptive Coefficients Estimation for Variational Pansharpening},
year={2022},
pages={DOI: 10.1109/TGRS.2022.3154480},
}
```

# Method

**Motivation:** Fig. 1. (a) Plot of the intensity for both the PAN image and the HRMS data
randomly choosing a row of the image; (b) Plot of the gradient intensity for
both the PAN image and the HRMS data considering blue (B), green (G), red
(R), and near-infrared (NIR) bands for the same row as in (a). It is worth to
be remarked that the behavior of PAN and HRMS images is more similar in
the gradient domain than in the intensity domain.

![motivation](fig-to-showw/fig1.jpg)




**Overall Framework:** The framework of our model. The details of our framework can be found in Sect. III.

![framework](fig-to-showw/fig2.jpg)



**Fidelity:** Fig. 4. (a) A close-up of the HRMS image of the Ple´iades dataset; (b)
Context-aware regions extracted from (a); (c) Scatter plot of ∇1X1 and ∇1P
for pixels belonging to Class 1 in (b); (d) A histogram of ∇1X1  ∇1P
for pixels belonging to Class 1 in (b). It is worth to be remarked that the
distribution of ∇1X1 and ∇1P is well described by a directly proportional
function, i.e., the black dotted line in (c). For other regions and bands, there is
a similar directly proportional relationship between HRMS and PAN images
in the gradient domain. Thus, the values of ∇X∇P are approximately the
same inside each region and for each band.

![framework](fig-to-showw/fig3.png)



**Handling the decimation operation:** Fig. 5. The graphic representation of (22) for a scale ratio equal to 4.
The white squares with a blank content indicate zero values. The first row
shows the processing from U to USST, which is equal to an element-wise
multiplication between DSST and U, i.e., the second row. It is worth noting
that DSST is produced from sparse matrices (i.e., the blue squares), whose
entries are 1 only in one position [44].

![framework](fig-to-showw/fig4.jpg)








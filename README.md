[![MIT License][license-shield]][license-url]

<!-- PROJECT LOGO -->
<br />
<p align="center">

  <h2 align="center">Quantifying COVID-19 importation risk in a dynamic network of domestic cities and international countries</h2>

  <p align="center">
Xiamen University; University of Illinois at Urbana Champaign; Pennsylvania State University; Nanjing Audit University; University of Wisconsin-Madison
  </p>
</p>

<!-- TABLE OF CONTENTS -->
## Table of Contents

* [Reference](#reference)
* [About the Project](#about-the-project)
* [File Descriptions](#file-descriptions)
* [License](#license)
* [Contact](#contact)
* [Acknowledgements](#acknowledgements)

<!-- Reference -->
## Reference
If you use this dataset and code in your research or applications, please refer to this source:

Han, X., Xu, Y., Fan, L., Huang, Y., Xu, M., and Gao, S. (2021). [Quantifying COVID-19 importation risk in a dynamic network of domestic cities and international countries.](https://doi.org/10.1073/pnas.2100201118) Proceedings of the National Academy of Sciences, Jun 2021, 118 (26), e2100201118, DOI: 10.1073/pnas.2100201118

```
@article{hou2020intra,
  title={Quantifying COVID-19 importation risk in a dynamic network of domestic cities and international countries},
  author={Han, Xiaoyi and Xu, Yilan and Fan, Linlin and Huang, Yi and Xu, Minhong and Gao, Song},
  journal={Proceedings of the National Academy of Sciences},
  year={2021},
  volume={118},
  number={26},
  pages={e2100201118},
  publisher={National Academy of Sciences}
}
```

<!-- ABOUT THE PROJECT -->
## About The Project

Since its outbreak in December 2019, the novel coronavirus 2019 (COVID-19) has spread to 191 countries and caused millions of deaths. Many countries have experienced multiple epidemic waves and faced containment pressures from both domestic and international transmission. In this study, we conduct a multi-scale geographic analysis of the spread of COVID-19 in a policy-influenced dynamic network to quantify COVID-19 importation risk under different policy scenarios using evidence from China. Our Spatial Dynamic Panel Data (SDPD) model explicitly distinguishes the effects of travel flows from the effects of transmissibility within cities, across cities, and across national borders. We find that within-city transmission was the dominant transmission mechanism in China at the beginning of the outbreak, and that all domestic transmission mechanisms were muted or significantly weakened before importation posed a threat. We identify effective containment policies by matching the change points of domestic and importation transmissibility parameters to the timing of various interventions. 
Our simulations suggest that importation risk is limited when domestic transmission is under control, but that cumulative cases would have been almost 13 times higher if domestic transmissibility had resurged to its pre-containment level after importation, and 32 times higher if domestic transmissibility had remained at its pre-containment level since the outbreak. Our findings provide practical insights into infectious disease containment and call for collaborative and coordinated global suppression efforts. 

## File Descriptions  
Written by: Han, Xiaoyi (han.293@buckeyemail.osu.edu.)

This Matlab program is to produce the counterfactual results in Figures 3 to 5 

Notes:
In this simulation code, we estimate the SDPD model with change points and simulated the counterfactual national cumulative cases, 
based upon the MCMC draws of parameter estimates. 

To implement the Bayesian 95% confidence interval, we include a function 'hpdi' provided by James P. LeSage in the file named 'jplv7', 
which can be downloaded  from [http://www.spatial-econometrics.com/](http://www.spatial-econometrics.com/)

We also include the data, as well as the time-varying spatial weights matrices used in the empirical study.

city_data_9_12_v14:        The excel file for data used in the empirical study
spatial_weights_9_29:     The M file for contemporaneous spatial weights matrix from population flows
spatial_weights_lag_9_29:  The M file for lagged spatial weights matrix from population flows

Main_empirical_figure3:  The code to generate figure 3. It reads the data and spatial weights, and implements mcmcbfnsbr1c3, which is a function to conduct MCMC estimation
for the SDPD model with change points and conduct counterfactual simulations for figure 3


Main_empirical_figure4:  The code to generate figure 4. It reads the data and spatial weights, and implements mcmcbfnsbr1c4, which is a function to conduct MCMC estimation
for the SDPD model with change points and conduct counterfactual simulations for figure 4


Main_empirical_figure5:  The code to generate figure 5. It reads the data and spatial weights, and implements mcmcbfnsbr1c5, which is a function to conduct MCMC estimation
for the SDPD model with change points and conduct counterfactual simulations for figure 5

jplv7: file for functions used in the econometric toolbox for matlab, which includes the function to generate Bayesian 95% confidence interval (written by James P. LeSage)

<!-- LICENSE -->
## License

Distributed under the MIT License. See `LICENSE` for more information.


<!-- CONTACT -->
## Contact
Xiaoyi Han - [Personal Website](https://hanxiaoyi.weebly.com/) - han.293 at buckeyemail.osu.edu </br>
Yilan Xu - [Personal Website](https://sites.google.com/site/xuyilan/) - yilanxu at llinois.edu </br>
Song Gao - [@gissong](https://twitter.com/gissong) - song.gao at wisc.edu  </br>


<!-- ACKNOWLEDGEMENTS -->
## Acknowledgements
The authors thank the China Data Lab (CDL) for providing the data [Resources for COVID-19 Study](https://projects.iq.harvard.edu/chinadatalab/resources-covid-19). Xiaoyi Han acknowledges the financial support of the National Natural Science Foundation of China (Awards No. 71973113 and No.71988101 ). Song Gao acknowledges the funding support provided by the U.S. National Science Foundation (Award No. BCS-2027375). Any opinions, findings, and conclusions, or recommendations expressed in this material are those of the author(s) and do not necessarily reflect the views of the funders. 

<!-- MARKDOWN LINKS & IMAGES -->
[license-shield]: https://img.shields.io/github/license/othneildrew/Best-README-Template.svg?style=flat-square
[license-url]: https://github.com/GeoDS/COVID19USFlows/blob/master/LICENSE.txt

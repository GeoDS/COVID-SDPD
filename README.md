[![MIT License][license-shield]][license-url]

<!-- PROJECT LOGO -->
<br />
<p align="center">
<a href="https://geods.geography.wisc.edu/">
    <img src="https://github.com/GeoDS/COVID19USFlows/raw/master/images/geods_safegraph_nsf_logo.jpg" alt="Logo" width="400">

  <h2 align="center">Quantifying COVID-19 importation risk in a dynamic network of domestic cities and international countries</h2>

  <p align="center">
Department of Mathematics, Department of Geography, Department of Life Sciences Communication, and School of Medicine and Public Health, 
University of Wisconsin-Madison
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

Hou, X., Gao, S., Li, Q., Kang, Y., Chen, N., Chen, K., Rao, J., Ellenberg, J.S. and Patz, J.A. (2021). [Intracounty modeling of COVID-19 infection with human mobility: assessing spatial heterogeneity with business traffic, age and race.](https://www.pnas.org/content/118/24/e2020524118) Proceedings of the National Academy of Sciences, Jun 2021, 118 (24), e2020524118, DOI: 10.1073/pnas.2020524118

```
@article{hou2020intra,
  title={Intra-county modeling of COVID-19 infection with human mobility: assessing spatial heterogeneity with business traffic, age and race},
  author={Hou, Xiao and Gao, Song and Li, Qin and Kang, Yuhao and Chen, Nan and Chen, Kaiping and Rao, Jinmeng and Ellenberg, Jordan S and Patz, Jonathan A},
  journal={Proceedings of the National Academy of Sciences},
  year={2021},
  volume={118},
  number={24},
  pages={e2020524118},
  publisher={National Academy of Sciences}
}
```

<!-- ABOUT THE PROJECT -->
## About The Project

Since its outbreak in December 2019, the novel coronavirus 2019 (COVID-19) has spread to 191 countries and caused millions of deaths. Many countries have experienced multiple epidemic waves and faced containment pressures from both domestic and international transmission. In this study, we conduct a multi-scale geographic analysis of the spread of COVID-19 in a policy-influenced dynamic network to quantify COVID-19 importation risk under different policy scenarios using evidence from China. Our Spatial Dynamic Panel Data (SDPD) model explicitly distinguishes the effects of travel flows from the effects of transmissibility within cities, across cities, and across national borders. We find that within-city transmission was the dominant transmission mechanism in China at the beginning of the outbreak, and that all domestic transmission mechanisms were muted or significantly weakened before importation posed a threat. We identify effective containment policies by matching the change points of domestic and importation transmissibility parameters to the timing of various interventions. 
Our simulations suggest that importation risk is limited when domestic transmission is under control, but that cumulative cases would have been almost 13 times higher if domestic transmissibility had resurged to its pre-containment level after importation, and 32 times higher if domestic transmissibility had remained at its pre-containment level since the outbreak. Our findings provide practical insights into infectious disease containment and call for collaborative and coordinated global suppression efforts. 

## File Descriptions  
A description of all files in the repository is shown below:  
data_0411_0814.mat: The daily COVID-19 confirmed cases data needed to run the model for Dane County, WI, provided by the [Wisconsin Department of Health Services](https://data.dhsgis.wi.gov/datasets/covid-19-historical-data-by-census-tract/data?orderBy=GEOID). (The codes of cleaning the data are provided in the folder ‘data_preparation’)

data_0311_0812.mat: The daily COVID-19 confirmed cases data needed to run the model for Milwaukee County, WI, provided by the [Wisconsin Department of Health Services](https://data.dhsgis.wi.gov/datasets/covid-19-historical-data-by-census-tract/data?orderBy=GEOID). (The codes of cleaning the data are provided in the folder ‘data_preparation’)

para_stochastic.m: A fixed-point iteration is used to determine the coefficients of the stochastic (Ornstein–Uhlenbeck) process.

fun_EKI_stochastic.m: Apply the Ensemble Kalman Filter.

SEIR_stochastic: Run the model forward

plots.m: Plot figures about the results of model.

case1.m, case2.m, case3.m: Three scenario studies of Dane County.


<!-- LICENSE -->
## License

Distributed under the MIT License. See `LICENSE` for more information.



<!-- CONTACT -->
## Contact
Xiaoyi Han - 
Yilan Xu - 
Song Gao - [@gissong](https://twitter.com/gissong) - song.gao at wisc.edu  </br>

Project Link: [https://github.com/GeoDS/IntraCounty-Mobility-SEIR](https://github.com/GeoDS/IntraCounty-Mobility-SEIR)  



<!-- ACKNOWLEDGEMENTS -->
## Acknowledgements
The authors thank the China Data Lab (CDL) for providing the data [Resources for COVID-19 Study](https://projects.iq.harvard.edu/chinadatalab/resources-covid-19). Xiaoyi Han acknowledges the financial support of the National Natural Science Foundation of China (Awards No. 71973113 and No.71988101 ). Song Gao acknowledges the funding support provided by the U.S. National Science Foundation (Award No. BCS-2027375). Any opinions, findings, and conclusions, or recommendations expressed in this material are those of the author(s) and do not necessarily reflect the views of the funders. 

<!-- MARKDOWN LINKS & IMAGES -->
[license-shield]: https://img.shields.io/github/license/othneildrew/Best-README-Template.svg?style=flat-square
[license-url]: https://github.com/GeoDS/COVID19USFlows/blob/master/LICENSE.txt
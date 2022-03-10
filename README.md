# wbc_sustained_obs

**Relevant code for Chandler et al. [in prep]. Western boundary current variability from sustained ocean observations.**
 
`XAA_ix21.m`, `XAA_px30.m`, and `XAA_px40.m` are the scripts to produce the velocity cross-sections for HR-XBT transects [IX21](http://www-hrx.ucsd.edu/ix15.html), [PX30](http://www-hrx.ucsd.edu/px31.html), and [PX40](http://www-hrx.ucsd.edu/px40.html).
 
`ix21_velocity.nc`, `px30_velocity.nc`, and `px40_velocity.nc` are the velocity cross-section data files output from the above scripts. *These files will be made permanently available through Zenodo upon publication of the manuscript.* In the meantime, preliminary versions for peer review can be found [in this Google Drive](https://drive.google.com/drive/folders/1UTvaPosz9Z--lX1g_9b-1OqlUm6AsMKI?usp=sharing).

`Fig1.m`, `Fig1_panels.m`, `Fig2.m`, `Fig3.m`, `Fig4.m`, and `Fig5.m` are the scripts used to produce the figures in the submitted manuscript.

`FigS1.m`, `FigS2.m`, and `FigS3.m` are the scripts used to produce the figures in the supporting information.

The `additional_functions` folder includes my functions I have used in the above scripts. In particular the function `XAA_processing.m` contains the bulk of the methodology (correcting for individual transect path differences, extending observations in depth, and producing regular monthly time series). Other functions I have used include those from [TEOS-10](http://www.teos-10.org/) and [MATLAB File Exchange](https://www.mathworks.com/matlabcentral/fileexchange/) - please contact me for specifics.

### Contact:
Mitchell Chandler  
Scripps Institution of Oceaongraphy, UC San Diego  
mlchandler@ucsd.edu  

---

https://mlchandler.github.io/

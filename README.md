# HCCA-aTO-Predictor
This repository provides the Shiny application source code for the dynamic prediction model of alternative Textbook Outcome (a-TO) in patients with perihilar cholangiocarcinoma (HCCA).
The model was developed using multi-center real-world clinical data and implemented using H2O Distributed Random Forest (DRF).
Online application:
https://hcca-to-predictor.shinyapps.io/shiny
Repository contents:
- app.R  
  Shiny application source code for interactive prediction.
- final_feature_names.rds  
  Final feature list used for model training and inference.
- final_factor_levels.rds  
  Factor level alignment metadata to ensure consistency between training and prediction.
Model availability:
Due to institutional data use agreements, the original dataset and trained model object are not publicly distributed. The trained H2O-DRF model may be provided for academic purposes upon reasonable request.
Reproducibility:
The modeling process, feature engineering procedures, and validation strategy are described in detail in the associated manuscript.
Disclaimer:
This tool is intended for research use only and does not replace clinical decision-making.

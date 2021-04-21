FROM inpred/tso500_base:v0.1

# directory "/inpred" is already present within the image

# copy InPreD scripts and PCGR/CPSR configuration files into the image
COPY internal_scripts /inpred/internal_scripts
COPY user_scripts /inpred/user_scripts
COPY configuration_files /inpred/configuration_files

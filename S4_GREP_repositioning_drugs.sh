
#--------------------
#
# GREP analysis 
# 
#--------------------

#Installation

#git clone https://github.com/saorisakaue/GREP
#cd ./GREP

module load anaconda/3.2019-10
#anaconda -h

cd ./GREP

python grep.py --genelist ./INSVD_analysis/mr_13_markers.txt --out mr_13_candidate_atc  --test ATC --output-drug-name
python grep.py --genelist ./INSVD_analysis/mr_13_markers.txt --out mr_13_candidate_icd  --test ICD --output-drug-name



 

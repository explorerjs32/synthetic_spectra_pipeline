#!/bin/csh -f

date
set mpath=/home/fmendez/Turbospectrum2019/COM-v19.1/models

foreach MODEL (model_5500g+4.6z-0.50vt01.TURBO)

#set Cabu = 8.656
#set Nabu = 7.78
#set Oabu = 8.66

set lam_min = '5000'
set lam_max = '5500'
set deltalam = '0.01'
set METALLIC = '-0.50'
set TURBVEL = '1.00'
set SUFFIX =
set result = ${SUFFIX}


# ABUNDANCES FROM THE MODEL ARE NOT USED !!!

/home/fmendez/Turbospectrum2019/exec-gf-v19.1/babsma_lu << EOF
'LAMBDA_MIN:'  '${lam_min}'
'LAMBDA_MAX:'  '${lam_max}'
'LAMBDA_STEP:' '${deltalam}'
'MODELINPUT:' '$mpath/${MODEL}'
'MARCS-FILE:' '.false.'
'MODELOPAC:' '/home/fmendez/Turbospectrum2019/COM-v19.1/contopac/${MODEL}opac'
'METALLICITY:'    '${METALLIC}'
'ALPHA/Fe   :'    '0.00'
'HELIUM     :'    '0.00'
'R-PROCESS  :'    '0.00'
'S-PROCESS  :'    '0.00'
'INDIVIDUAL ABUNDANCES:'   '0'
'XIFIX:' 'T'
$TURBVEL
EOF

########################################################################

/home/fmendez/Turbospectrum2019/exec-gf-v19.1/bsyn_lu <<EOF
'LAMBDA_MIN:'     '${lam_min}'
'LAMBDA_MAX:'     '${lam_max}'
'LAMBDA_STEP:'    '${deltalam}'
'INTENSITY/FLUX:' 'Flux'
'COS(THETA)    :' '1.00'
'ABFIND        :' '.false.'
'MODELOPAC:' '/home/fmendez/Turbospectrum2019/COM-v19.1/contopac/${MODEL}opac'
'METALLICITY:'    '${METALLIC}'
'ALPHA/Fe   :'    '0.00'
'HELIUM     :'    '0.00'
'R-PROCESS  :'    '0.00'
'S-PROCESS  :'    '0.00'
'RESULTFILE :' '/home/fmendez/Turbospectrum2019/COM-v19.1/syntspec/${result}'
'INDIVIDUAL ABUNDANCES:'   '2'
11 5.99
12 8.10
'ISOTOPES : ' '2'
3.006  0.075
3.007  0.925
'NFILES   :' '2'
DATA/Hlinedata
/home/fmendez/Documents/Research/Codes/linelist_handling/Fe-Mg-Na_5000-6000.lines
'SPHERICAL:'  'F'
  30
  300.00
  15
  1.30
EOF
########################################################################
date
end

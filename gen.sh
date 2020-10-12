# Pick a mechanism
#python -m pyjac -i Mechanisms/sanDiego_red.xml -l c -jf full -conv -b .
#python -m pyjac -i Mechanisms/sanDiego.xml -l c -jf full -conv -b .
python -m pyjac -i Mechanisms/SanDiegoN9.xml -l c -jf full -conv -b .
rm -r -f LeapIMEXMethod.H
rm -r -f LeapMRIMEXMethod.H
# Differential or algebraic treatment?
#python leap.pyz leap_imex_gen.py LeapIMEXMethod.H
python leap.pyz leap_imex_alg_gen.py LeapIMEXMethod.H
# IMEX-AM placeholder
python leap.pyz leap_mr_imex_gen.py LeapMRIMEXMethod.H

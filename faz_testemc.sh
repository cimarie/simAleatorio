#!/bin/bash

delta="0.001"
ngroup="5000"
nindiv="26"
mu="0.0001"
for b in "0.0" "1.0" "2.0" "5.0"; do
    for alpha in "0.2" "0.8" "1.0"; do 
        cd ~/Documents/simAleatorio
        python workflow_teste.py RodaTodosTestes --local-scheduler --benefit $b --alpha $alpha --nindiv 26 --ngroup 5000 --mu 0.0001 --delta $delta
        mv teste_mc_c\=*.txt figuras_mc
        mv teste_mc.log figuras_mc 
        cd figuras_mc
        python plotar_fig.py imprime_figuras --b $b --delta $delta --alpha $alpha --n $nindiv
        mkdir b\=${b}_delta\=${delta}_alpha\=${alpha}
        mv *.txt *.log *.png  b\=${b}_delta\=${delta}_alpha\=${alpha} 
    done 
done

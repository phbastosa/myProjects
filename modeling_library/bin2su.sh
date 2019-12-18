#!/bin/bash

entrada="Sismogramas_Trecho2.bin"    # nome do arquivp de entrada
nt=5000                              # número total de amostra no tempo
dt=120                               # discretização temporal em micro segundos
xf_0=188                             # posição do primeiro tiro em metros
df=8                                 # espaçamento entre fontes
nrec=72                              # quantidade de traços por tiro
rec_0=-188                           # distância do primeiro traço em relação à fonte  
inc_rec=8                            # espaçamento entre os receptores em metros
saida="Sismogramas_Trecho2.su"       # nome do arquivo de saída

suaddhead <$entrada ns=$nt | 
sushw key=dt a=$dt | 
sushw key=sx a=$xf_0 b=0 c=$df j=$nrec | 
sushw key=offset a=$rec_0 b=10 c=0 j=$nrec | 
suchw key1=gx key2=offset key3=sx b=1 c=1 d=1 | 
suchw key1=cdp key2=sx key3=gx b=1 c=1 d=2 | 
sushw key=fldr a=1001 b=0 c=1 >$saida

# Para visualizar o Header
sugethw <$saida sx gx offset | more



















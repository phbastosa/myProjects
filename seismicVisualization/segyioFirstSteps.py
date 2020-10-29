"""
Código: GGOP0003
Projeto: Desenvolvimento de habilidades computacionais e práticas do processamento sísmico

Conteúdo do algoritmo:
    
    - Primeiros passos usando a biblioteca segyio 
    - Leitura de arquivo SEGY contendo dados sísmicos pré-empilhados
    - Coleta de atributos do Header de arquivos SEGY
    - Plotando janelas de sismogramas no domínio do tiro

Construção: 28 / 09 / 2020
Professor: Marco Cetale
Disciplina: GGO00060 - Processamento Sísmico - A1

Monitor: Paulo Henrique Bastos 
"""
#%% Área do usuário - Opções binárias 0 [= False] ou 1 [= True]

inputFile = "seismic.segy"                               # Caminho do arquivo de entrada

showBinHeader = 0                                        # Mostrar o Header do conjunto de traços no terminal 
showTraceHeader = 0                                      # Mostrar o Header de cada traço no terminal
 
iSeismogram = 0                                          # Sismograma inicial levando em consideração o número de tiros 
fSeismogram = 0                                          # Sismograma final para ser mostrado na figura

colorRescal = 99                                         # Parâmetro para realçar as cores do sismograma, opções comuns 99,98,97...90 sendo 99 o realce mais intenso
showFigure = 0                                           # Para mostrar a figura num image plot

saveFigure = 0                                           # Para salvar a figura
filename = "shots.png"                                   # Nome do arquivo .PNG

#%% Importando as bibliotecas necessárias

import segyio                                             # Para leitura e edição de arquivos SEGY
import numpy as np                                        # Para gerar e editar arrays 
import matplotlib.pyplot as plt                           # Para gerar gráficos e imagens

from skimage import exposure                              # Usada para realçar as cores de imagens

#%% Leitura do dado sísmico

data = segyio.open(inputFile,ignore_geometry=True)   # Comando para a leitura

#%% Exportando as informações dos Headers (Adaptado de Felipe Timóteo - GISIS)

binheader = segyio.binfield.keys                          # Coletando as chaves do Binary Header

if showBinHeader:
    print("\n Check binary header \n")                    # Informando o que será exposto
    print("%25s %4s %5s \n" %("key","byte","value"))      # Associando posições e titulações
    for k,v in binheader.items():                         # Loop indexando chaves e valores
        if v in data.bin:                                 # Verificando se existem algum elemento da chave indicada
            print("%25s %d  %d" %(k,v,data.bin[v]))       # Mostrando informações

traceheaders = segyio.tracefield.keys                     # Coletando as chaves do Trace Header

if showTraceHeader:
    print("\n Check trace header \n")                     # Informando o que será exposto 
    print("%40s %5s %6s %6s \n" %("Trace header", "byte", # Associando posições e titulações
                                         "first", "last")) 
    for k,v in traceheaders.items():                      # Loop indexando chaves e velores    
        aux1 = data.attributes(v)[0]                      # Coletando o primeiro atributo associado
        aux2 = data.attributes(v)[data.tracecount-1]      # Coletando o último atributo associado
        print("%40s %5d  %6d %6d " %(k,v,aux1,aux2))      # Mostrando informações

#%% Extraindo a quantidade de fontes e receptores por tiro da aquisição

f = data.attributes(13)[data.tracecount-1][0]             # Coletando a informação final 
i = data.attributes(13)[0][0]                             # Coletando a informação inicial
nRec = f - i + 1                                          # Efetuando e registrando a informação do conteúdo de receptores por tiro

f = data.attributes(9)[data.tracecount-1][0]              # Coletando a informação final
i = data.attributes(9)[0][0]                              # Coletando a informação inicial
nSrc = f - i + 1                                          # Efetuando e registrando a informação do total de tiros da aquisição

#%% Preparando a janela de tiros 

shots = np.arange(nSrc)                                   # Array com o número de tiros

w = [iSeismogram , fSeismogram]                           # Janela a ser mostrada 
wind = slice(w[0]*nRec,nRec + w[1]*nRec)                  # Realizando o slice para ser parametro do objeto imagem

imgObj = data.trace.raw[wind].T                           # Coletando os dados de amplitude numa janela desejada 

perc = np.percentile(imgObj,[.5, colorRescal])            # Adicionando um ganho nas cores e as realçando
image = exposure.rescale_intensity(imgObj,                # Utilizando a biblioteca com as devidas sintaxes
               in_range=(perc[0],perc[1]),                # Adicionando o perc criado anteriormente 
                        out_range=(0,255))                # Range de cores padrão 

#%% Construindo atributos da figura 

nt = data.attributes(115)[0][0]                           # Número de amostras no tempo de cada traço
dt = data.attributes(117)[0][0] / 1e6                     # Tempo de uma amostra para outra - parâmetro de discretização temporal

locks_t = np.linspace(0,nt,10)                            # Para adicionar pontos de escala de tempo
label_t = np.linspace(0,nt,10) * dt                       # Calculando novos pontos de escala de tempo

ishot = data.attributes(9)[w[0]*nRec][0]                  # Número inicial do tiro analizado de acordo com o Header 
fshot = data.attributes(9)[w[1]*nRec][0]                  # Número final do tiro analizado de acordo com o Header 

#%% Gerando imagem com a janela de tiros atribuída

if showFigure:
    plt.figure(f"Shot {ishot} - {fshot}",figsize=(11,6))  # Título da figura com seu respectivo tamanho
    cbar = plt.colorbar(plt.imshow(imgObj,cmap="Greys"))  # Escala de cor sem o realce, conservando a escala do dado original 
    cbar.set_label("Amplitudes")                          # Atribuindo um título para a escala de cores

    plt.imshow(image, aspect="auto", cmap="Greys")        # Plot dos dados  

    plt.title("Shot gather window",fontsize=15)           # Título da figura
    plt.ylabel("Two way time [s]",fontsize=15)            # Título do eixo y 
    plt.xlabel("Trace number",fontsize=15)                # Título do eixo x

    plt.yticks(locks_t,np.around(label_t,decimals=1))     # Aplicando escala de tempo atribuida anteriormente

    if saveFigure: 
        plt.savefig(filename,dpi=200,bbox_inches="tight") # Salvando a figura
    
    plt.show()                                            # Mostrando a figura

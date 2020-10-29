"""
Código: GGOP0003

Projeto: Desenvolvimento de habilidades computacionais e 
práticas do processamento sísmico

Biblioteca de funções para ser usada nas aulas práticas de 
processamento sísmico.

Construção: 28/09 a --/12 de 2020.

Professor: Marco Cetale

Disciplina: GGO00060 - Processamento Sísmico - A1

Monitor: Paulo Henrique Bastos 
"""
#%% Importando as bibliotecas necessárias
import segyio                                                 # Para leitura e edição de arquivos SEGY
import numpy as np                                            # Para gerar e editar arrays 
import matplotlib.pyplot as plt                               # Para gerar gráficos e imagens

from skimage import exposure                                  # Usada para realçar as cores de imagens
#%% Função para ler arquivos segy usando segyio
def readSegyFile(sesimicFileName):
    '''
    Leitura de arquivos .segy ou .sgy

    input:
        seismicFileName - Caminho do arquivo a ser lido

    output:
        data - objeto de dado sísmico do segyio 
    '''
    data = segyio.open(sesimicFileName,ignore_geometry=True)  # Importando o dado na sintaxe do segyio
 
    return data                                               # Retornando objeto de dado sísmico do segyio 
#%% Função para exportar informações do binary header
def showBinaryHeader(data):
    '''
    Mostra na tela do terminal as informações do Binary Header
    
    input:
        data - objeto de dado sísmico do segyio
    '''
    binheader = segyio.binfield.keys                          # Coletando as chaves do Binary Header
    
    print("\n Check binary header \n")                        # Informando o que será exposto
    print("%25s %4s %5s \n" %("key","byte","value"))          # Associando posições e titulações
    for k,v in binheader.items():                             # Loop indexando chaves e valores
        if v in data.bin:                                     # Verificando se existem algum elemento da chave indicada
            print("%25s %d  %d" %(k,v,data.bin[v]))           # Mostrando informações
#%% Função para exportar informações do trace header
def showTraceHeader(data):
    '''
    Mostra na tela do terminal as informações do Trace Header
    
    input:
        data - objeto de dado sísmico do segyio
    '''
    traceheaders = segyio.tracefield.keys                     # Coletando as chaves do Trace Header

    print("\n Check trace header \n")                         # Informando o que será exposto 
    print("%40s %5s %6s %6s \n" %("Trace header", "byte",     # Associando posições e titulações
                                         "first", "last")) 
    for k,v in traceheaders.items():                          # Loop indexando chaves e velores    
        aux1 = data.attributes(v)[0]                          # Coletando o primeiro atributo associado
        aux2 = data.attributes(v)[data.tracecount-1]          # Coletando o último atributo associado
        print("%40s %5d  %6d %6d " %(k,v,aux1,aux2))          # Mostrando informações
#%% Função para realçar as cores dos sismogramas
def perc(matrix,value):
    '''
    Realça a escala de cores de uma matriz.
    
    Valores comuns: 99, 98, 97, ... 90. 
    
    Em uma escala de menor (99) para maior (90) realce. 

    input:
        matrix - matriz 2D para receber um ganho nos valores    
        value  - grandeza do ganho 

    output:
        image - matriz 2D com escala de cores realçada.
    '''
    p = np.percentile(matrix,[.5, value])                     # Adicionando um ganho nas cores e as realçando
    image = exposure.rescale_intensity(matrix,                # Utilizando a biblioteca com as devidas sintaxes
                         in_range=(p[0],p[1]),                # Adicionando o perc criado anteriormente 
                         out_range=(0,255))                   # Range de cores padrão 

    return image                                              # Retorna a imagem rescalada 
#%% Função para construir janela de shots em sequência
def buildShotWindow(data, ishot, fshot):
    '''
    Monta uma matriz com os tiros configurados em sequência 
    (ishot < fshot) de acordo com os identificadores do trace header
    
    input:
        data  - objeto de dado sísmico do segyio
        ishot - identificação do tiro inicial 
        fshot - identificação do tiro final    

    output:
        shotRange - matriz 2D com seguência de sismogramas 
    '''
    f = data.attributes(13)[data.tracecount-1][0]             # Coletando a informação final 
    i = data.attributes(13)[0][0]                             # Coletando a informação inicial
    nRec = f - i + 1                                          # Efetuando e registrando a informação do conteúdo de receptores por tiro

    f = data.attributes(9)[data.tracecount-1][0]              # Coletando a informação final
    i = data.attributes(9)[0][0]                              # Coletando a informação inicial
    nSrc = f - i + 1                                          # Efetuando e registrando a informação do total de tiros da aquisição

    w = [ishot , fshot]                                       # Janela a ser mostrada 
    wind = slice(w[0]*nRec,nRec + w[1]*nRec)                  # Realizando o slice para ser parametro do objeto imagem

    shotRange = data.trace.raw[wind].T                        # Coletando os dados de amplitude numa janela desejada 

    return shotRange                                          # retorna a janela de shots configurados
#%% Função para plotar janela de shots em sequência 
def plotShotWindow(data, shotRange, rescalePerc, figName):
    '''
    Gera figura de sequência de shots, ou shot único. Tem a opção 
    de reescalar as cores e salvar a figura.     
    
    input:
        data        - objeto de dado sísmico do segyio
        shotRange   - matriz 2D com seguência de sismogramas
        rescalePerc - parametro para ganho de cores na imagem
        figName     - caminho para armazenar a figura gerada
    
    output:         
        imagem .png com a sequência de shots configurada
    '''    
    nt = data.attributes(115)[0][0]                           # Número de amostras no tempo de cada traço
    dt = data.attributes(117)[0][0] / 1e6                     # Tempo de uma amostra para outra - parâmetro de discretização temporal

    locks_t = np.linspace(0,nt,10)                            # Para adicionar pontos de escala de tempo
    label_t = np.linspace(0,nt,10) * dt                       # Calculando novos pontos de escala de tempo

    image = perc(shotRange,rescalePerc)                       # Adicionando o ganho nas cores da figura

    plt.figure(1,figsize=(11,6))                              # Título da figura com seu respectivo tamanho
    cbar = plt.colorbar(plt.imshow(shotRange,cmap="Greys"))   # Escala de cor sem o realce, conservando a escala do dado original 
    cbar.set_label("Amplitudes")                              # Atribuindo um título para a escala de cores

    plt.imshow(image, aspect="auto", cmap="Greys")            # Plot dos dados  
    plt.title("Shot gather window",fontsize=15)               # Título da figura
    plt.ylabel("Two way time [s]",fontsize=15)                # Título do eixo y 
    plt.xlabel("Traces",fontsize=15)                          # Título do eixo x
    plt.yticks(locks_t,np.around(label_t,decimals=1))         # Aplicando escala de tempo atribuida anteriormente
    plt.savefig(figName,dpi=200,bbox_inches="tight")          # Salvando a figura
    plt.show()                                                # Mostrando a figura
#%% Plotando a organização dos tiros 
def plotAcquisitionGeometry(data,ishot,fshot,figName):
    '''
    Scatter plot da sequência de tiros (de acordo com o 
    trace header) em coordenadas relativas da geometria de 
    aquisição do dado    

    input:
        data    - objeto de dado sísmico do segyio 
        ishot   - indicador de tiro inicial
        fshot   - indicador de tiro final
        figName - caminho para armazenar a figura gerada

    output:
        imagem .png com a sequência de shots configurada
    '''   
    initShotNumber = data.attributes(9)[0][0]                 # Coletando o ponto de tiro inicial de acordo com o trace header 
    sx = data.attributes(73)[:]                               # Coletando todos os pontos de tiro em coordenadas relativas
    gx = data.attributes(81)[:]                               # Coletando todos os pontos de receptores em coordenadas relativas
    spread = data.attributes(13)[data.tracecount-1][0]        # Coletando o numero de receptores ativos por tiro
    psx = sx[::spread]                                        # Coletando cada ponto de tiro em X descartando os desnecessários
    pst = np.arange(len(psx)) + initShotNumber                # Fazendo um array para entrar como eixo Y na imagem, serão os identificadores de tiro

    psg = np.zeros(len(gx))                                   # Inicializando um array de posições de receptores 
    for p in range(len(pst)):                                 # Loop para coletar cada receptor por tiro
        psg[p*spread:p*spread + spread] = pst[p]              # Armazenando organizado para plotar 

    flagti = int(ishot - initShotNumber)                      # Ajustando parametros iniciais de identificador de tiro para encaixar na figura
    flagtf = int(fshot - initShotNumber + 1)                  # Ajustando parametros finais de identificador de tiro para encaixar na figura

    gslice = slice(flagti*spread,flagtf*spread)               # Realizando o range de receptores escolhido pelo usuário
    sslice = slice(flagti,flagtf)                             # Realizando o range de fontes escolhido pelo usuário

    plt.figure(1,figsize=(10,5))                              # Construindo figura
    plt.scatter(gx[gslice],psg[gslice],marker=".")            # Plotando receptores 
    plt.scatter(psx[sslice],pst[sslice],marker=".")           # Plotando shots
    plt.title("Geometria de aquisição",fontsize=15)           # Titulo da figura
    plt.xlabel("Distância [m]",fontsize=15)                   # Titulo do eixo x da figura 
    plt.ylabel("Indicador de sismograma",fontsize=15)         # Titulo do eixo y da figura

    plt.legend(["Recs","Shots"],loc="lower left")             # Legenda para indicar quais pontos são shots e quais são recs
    plt.gca().invert_yaxis()                                  # Invertendo o eixo para o menor valor no y começar de cima
    plt.savefig(figName,dpi=200,bbox_inches="tight")          # Salvando a figura
    plt.show()                                                # Mostrando a figura na tela
#%% Mostra as possibilidades de offset
def showOffsetPossibilities(data): 
    '''
    Mostra a quantidade de offsets e seus respectivos valores 
    em metro. Função usada para contruir seção de offset comum

    input:
        data - objeto de dado sísmico do segyio

    output:
        sections - array 1D com os valores de offset em metro
    '''
    ntrace = data.tracecount                                  # Coletando o número de traços na aquisição
    sx = data.attributes(73)[:]                               # Coletando as posições de tiros por traço
    gx = data.attributes(81)[:]                               # Coletando as posições de receptores por traço
    offset = sx - gx                                          # Calculando o offset para todos os traços

    sections = np.array([])                                   # Inicializando array para a coleta de offsets
    for i in range(ntrace):                                   # Loop varrendo todos os traços 
        if offset[i] not in sections:                         # Se offset atual ainda não foi coletado
            sections = np.append(sections,offset[i])          # Coletando offsets não coletados 

    return sections[::-1]                                     # Retornando array 1D com os offsets em ordem crescente (offset[0] = Near Offset)
#%% Montando seção de offset comum
def buildOffsetSection(data,setOffset):
    '''
    Constrói uma seção de offset comum dado um valor de offset 
    que contempla o header do dado - usar showOffsetPossibilities(data)   

    input:
        data      - objeto de dado sísmico do segyio
        setOffset - indicador de offset de valor único     

    output:
        section - matriz 2D sendo a seção de offset comum 
    '''
    nt = data.attributes(115)[0][0]                           # Coletando a quantidade de pontos no tempo do traço
    dt = data.attributes(117)[0][0] / 1e6                     # Tempo de uma amostra para outra - parâmetro de discretização temporal

    sx = data.attributes(73)[:]                               # Coletando as posições de tiros por traço
    gx = data.attributes(81)[:]                               # Coletando as posições de receptores por traço
    offset = sx - gx                                          # Calculando a distância tiro-receptor por traço

    search = np.where(setOffset == offset)[0]                 # Coletando todos os offsets compatíveis com o configurado
    section = np.zeros((nt,len(search)))                      # Inicializando o array de traços de offset comum
    for i in range(len(search)):                              # Loop varrendo todos os traços coletados 
        section[:,i] = data.trace.raw[search[i]]              # Encaixando cada traço na seção

    return section                                            # Retornando a seção completa
#%% Plot de seção do offset commum
def plotOffsetSection(data,section,setOffset,rescalePerc,figName):
    '''
    Plotando a seção de offset comum construida

    input:
        data        - objeto de dado sísmico do segyio  
        setion      - matriz 2D contendo a seção de offset comum
        setOffset   - offset configurado da seção
        rescalePerc - parametro para ganho de cores na imagem 
        figName     - caminho para armazenar a figura gerada

    output:
        imagem .png com a seção de offset comum configurada
    '''
    nt = data.attributes(115)[0][0]                           # Coletando a quantidade de pontos no tempo do traço
    dt = data.attributes(117)[0][0] / 1e6                     # Tempo de uma amostra para outra - parâmetro de discretização temporal

    img = perc(section,rescalePerc)                           # Adicionando o ganho nas cores da figura

    locks_t = np.linspace(0,nt,10)                            # Para adicionar pontos de escala de tempo
    label_t = np.linspace(0,nt,10) * dt                       # Calculando novos pontos de escala de tempo
    locks_x = np.linspace(0,len(section[0]),11,dtype=int)     # Delimitando posições para a escala de distância
    isn = data.attributes(9)[0][0]                            # Indicador inicial de tiro
    fsn = data.attributes(9)[data.tracecount-1][0]            # Indicador final de tiro
    label_x = np.linspace(isn,fsn,11,dtype=int)               # Montando array para encaixar na escala de distância 

    plt.figure(1,figsize=(10,5))
    cbar = plt.colorbar(plt.imshow(section,cmap="Greys"))     # Escala de cor sem o realce, conservando a escala do dado original 
    cbar.set_label("Amplitudes")                              # Atribuindo um título para a escala de cores    
    plt.imshow(img,cmap="Greys",aspect="auto")                # Plotando a imagem em escala de cinza
    plt.yticks(locks_t,np.around(label_t,decimals=1))         # Aplicando escala de tempo atribuida anteriormente
    plt.xticks(locks_x,label_x)                               # Aplicando escala de distancia construida anteriormente
    plt.title(f"Offset comum de {setOffset} m",fontsize=15)   # Titulo da figura
    plt.xlabel("Número do tiro",fontsize=15)                  # Titulo do eixo x
    plt.ylabel("Tempo [s]",fontsize=15)                       # Titulo do eixo y
    plt.savefig(figName,dpi=200,bbox_inches="tight")          # Salvando a figura
    plt.show()                                                # Mostrando a figura na tela
#%% Mostrando as possibilidades de posições de CMP  
def showCMPsPossibilities(data):
    '''
    Mostra as possibilidades de posições de ponto médio comum 
    na aquisição de acordo com a geometria e coordenadas relativas 

    input:
        data - objeto de dado sísmico do segyio

    output:
        line - array 1D com as posições dos CMPs em metros
    '''
    ntrace = data.tracecount                                  # Coletando a quantidade de traços na aquisição      
    sx = data.attributes(73)[:]                               # Coletando as posições de tiros por traço
    gx = data.attributes(81)[:]                               # Coletando as posições de receptores por traço

    cmps = sx - (sx - gx) / 2                                 # Calculando todos os CMPs por traço

    line = np.array([])                                       # Inicializando o array das posições dos CMPs
    for i in range(ntrace):                                   # Loop varrendo todos os traços
        if cmps[i] not in line:                               # Se o CMP ainda não foi coletado
            line = np.append(line,cmps[i])                    # Adicione o CMP na linha     

    return line                                               # Retorna a linha com as posições dos CMPs
#%% Plotando a organização dos CMPs
def plotCMPsPosition(data,line,figName):
    '''
    Mostra a configuração dos CMPs na aquisição, revelando os traços
    por ponto médio comum que a geometria conseguiu realizar

    input:
        data    - objeto de dado sísmico do segyio 
        line    - array com as posições dos CMPs 
        figName - caminho para armazenar a figura gerada

    output:
        imagem .png mostrando o conteúdo de traços por CMP
    '''
    ntrace = data.tracecount                                  # Coletando a quantidade de traços na aquisição      
    sx = data.attributes(73)[:]                               # Coletando as posições de tiros por traço
    gx = data.attributes(81)[:]                               # Coletando as posições de receptores por traço

    cmps = sx - (sx - gx) / 2                                 # Calculando os pontos médios comuns por traço

    tracSum = np.array([])                                    # Inicializando array 1D para contar os traços por CMP
    for trace in line:                                        # Loop varrendo todos os pontos possíveis de posições CMP
        x = np.where(trace == cmps)                           # Procura a quantidade de traços pertencentes àquele CMP
        tracSum = np.append(tracSum,len(x[0]))                # Adiciona ao array organizando cada CMP com seus traços 

    plt.figure(1,figsize=(10,5))                              # Inicializando uma figura
    plt.plot(line,tracSum)                                    # Plotando conteúdo de traços por posição de CMP
    plt.title("Conteúdo de traços por CMP",fontsize=15)       # Titulo da figura
    plt.xlabel("Distância [m]",fontsize=15)                   # Titulo do eixo x da figura
    plt.ylabel("Quantidade de traços",fontsize=15)            # Titulo do eixo y da figura
    plt.savefig(figName,dpi=200,bbox_inches="tight")          # Salvando a figura
    plt.show()                                                # Mostrando a figura na tela
#%% Constrói um grupo de traços por CMP 
def buildCMPSection(data,setCMP):
    '''
    Constrói uma seção de ponto médio comum usando um CMP existente
    na aquisição em metros - usar showCMPsPossibilities(data)

    input:
        data   - objeto de dado sísmico do segyio
        setCMP - CMP existente na aquisição         

    output:
        gather - seção de ponto médio comum configurada
    '''
    nt = data.attributes(115)[0][0]                           # Coletando as amostras temporais dos traços
    sx = data.attributes(73)[:]                               # Coletando as posições de tiros por traço
    gx = data.attributes(81)[:]                               # Coletando as posições de receptores por traço
    cmps = sx - (sx - gx) / 2                                 # Coletando as posições de CMP por traço

    x = np.where(setCMP == cmps)                              # Procurando e achando os traços pertencentes àquela posição CMP
    gather = np.zeros((nt,len(x[0])))                         # Inicializando o array 2D com o conteúdo de amplitude dos traços

    for i in range(len(x[0])):                                # Loop varrendo os traços coletados
        gather[:,i] = data.trace.raw[x[0][i]]                 # Adicionando o conteúdo de amplitude por traço no gather

    return gather                                             # Retornando gather CMP 
#%% Plotando familia de pontos médios comum
def plotCMPGather(data,gather,setCMP,rescalePerc,figName):
    '''
    Plota o gather CMP criado

    input:
        data        - objeto de dado sísmico do segyio
        gather      - matriz 2D contendo os traços por CMP 
        setCMP      - posição em metros do CMP escolhido
        rescalePerc - parâmetro para ganho de cores na imagem 
        figName     - caminho para armazenar a figura gerada

    output:
        imagem .png mostrando as amplitudes dos traços no CMP

    '''
    nt = data.attributes(115)[0][0]                           # Coletando as amostras temporais dos traços
    dt = data.attributes(117)[0][0] / 1e6                     # Tempo de uma amostra para outra - parâmetro de discretização temporal

    img = perc(gather,rescalePerc)                            # Adicionando o ganho nas cores da figura

    locks_t = np.linspace(0,nt,10)                            # Para adicionar pontos de escala de tempo
    label_t = np.linspace(0,nt,10) * dt                       # Calculando novos pontos de escala de tempo

    plt.figure(1,figsize=(8,6))                               # Inicializando uma figura
    cbar = plt.colorbar(plt.imshow(gather,cmap="Greys"))      # Escala de cor sem o realce, conservando a escala do dado original 
    cbar.set_label("Amplitudes")                              # Atribuindo um título para a escala de cores    
    plt.imshow(img,cmap="Greys",aspect="auto")                # Plotando a imagem propriamente dita
    plt.title(f"CMP gather na posição {setCMP} m",fontsize=15)# Titulo da figura
    plt.xlabel("Traços",fontsize=15)                          # titulo do eixo x da figura
    plt.ylabel("Tempo [s]",fontsize=15)                       # Titulo do eixo y da figura
    plt.yticks(locks_t,np.around(label_t,decimals=1))         # Aplicando escala de tempo atribuida anteriormente
    plt.savefig(figName,dpi=200,bbox_inches="tight")          # Salvando a figura
    plt.show()                                                # Mostrando a figura na tela
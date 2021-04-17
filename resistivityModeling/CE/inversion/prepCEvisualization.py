#%% Cabeçalho do arquivo 
'''
    Algoritmo para visualizar geometria de arranjos 
    dipolo-dipolo em aquisições de eletrorresisitividade.
    
    Funcionalidades principais: 
    
    - Leitura da tabela de dados de campo feita pela Horizonte 
      para levantamentos manuais de caminhamento elétrico. 
    
    - Gera um arquivo de dados para ser usado no programa RES2DINV 
      com a sintaxe correta. 
      
    - Gera uma imagem dos pontos em profundidade com seus respectivos 
      valores de resistividade aparente. Caso não haja dados, o
      algoritmo pode ser usado para gerenciamento de aquisições 
      calculando a quantidade de pontos na linha de acordo com o 
      espaçamento entre os eletrodos e os níveis de investigação.  

    Um legado que eu deixo para a Horizonte Soluções Geofísicas Jr.
    
    3 de outubro de 2020.

    Paulo Henrique Bastos Alves - Diretor de projetos na gestão de 2018. 
'''
#%% Importando as bibliotecas necessárias

import numpy as np                        # Realização de operações numéricas com arrays 
import matplotlib.pyplot as plt           # Plotar figuras
import xlrd                               # Ler arquivos .xlsx
#%% Apresentação do arquivo no terminal

print(38*"-=")
print("Dados adquiridos com arranjo dipolo-dipolo para tomografia de resistividade")
print(38*"-=")
#%% Área do usuário

lineLength = 100.0                    # Tamanho completo da linha de caminhamento elétrico
electSpace = 5.0                      # Espaçamento unitário dos eletrodos 
invesLevel = 5                        # Níveis de investigação desejados

justSee = False                       # Flag para gerenciamento de aquisições, quando verdadeira, será mostrada somente o arranjo com pontos sem valores

table = 2                             # Tabela do excel que se deseja gerar o arquivo para inversão
sheetPath = "dataBase.xlsx"     # Caminho da planilha completa com os dados de campo

outputPath = "linha3.dat"             # Nome do arquivo de saida com os dados formatados para a inversão
title = "H37 - Linha3 01/10/2020"     # Título que será colocado no arquivo de saída, importante colocar a data de aquisição
#%% Construindo geometria de aquisição

nElectrodes = int(lineLength / electSpace + 1)                              # Cálculo do número de eletrodos na linha 

diags = (nElectrodes - (invesLevel + 2)) * invesLevel                       # Quantidade de diagonais completas
others = sum(range(1,invesLevel))                                           # Calculando outros pontos que não formam diagonais
nPoints = diags + others                                                    # Calculando pontos completos na aquisição

electrodesPosition = np.linspace(0,lineLength,nElectrodes)                  # Calculando a posição dos eletrodos de acordo com o tamanho da linha

displPointsPosition = np.array([])                                          # Inicializando array para calcular a projeção de cada ponto em X
depthPointsPosition = np.array([])                                          # Inicializando array para calcular a profundidade em cada ponto

for i in range(invesLevel):
    displInc = 3*electSpace/2 + i*electSpace/2                              # Incremento de cada posição em X
    levelPoints = np.arange(displInc,lineLength - displInc + 1,electSpace)  # Cálculo de acordo com o nivel de investigação, começando no primeiro
    displPointsPosition = np.append(displPointsPosition,levelPoints)        # Montando o array completo de posições em X  

    depthInc = i*electSpace/2                                               # Incremento de cada posição em profundidade
    depth = np.ones(len(levelPoints)) * electSpace                          # Cálculo de acordo com o ivel de investigação
    depthPointsPosition = np.append(depthPointsPosition,depth+depthInc)     # Montando o array completo de profundidades 
#%% Lendo dados de campo - se houver 

if (not justSee):                                                           # Lendo somente se a flag justSee for falsa, ou seja, existe uma tabela direcionada para a leitura       
    wb = xlrd.open_workbook(sheetPath)                                      # Abrindo a tabela de dados de campo
    sheet = wb.sheet_by_index(table)                                        # Selecionando a aba da planilha, setar na area do usuário  
    
    index = sheet.col_values(1)[5:]                                         # Coletando os níveis de investigação da planilha
    apparentResistivity = sheet.col_values(19)[5:]                          # Coletando a resistividade aparente da planilha
    orderedData = np.array([],dtype=float)                                  # Inicializando um array para ordenar os dados por nível de investigação 

    o = open(outputPath, "w")                                               # Abrindo um arquivo texto
    o.write(title+"\n")                                                     # Escrevendo o título do arquivo configurado na area do usuário
    o.write(str(electSpace)+"\n")                                           # Escrevendo o espaçamento dos eletrodos  
    o.write("3\n")                                                          # Escrevendo a flag para arranjo Dipolo-dipolo 
    o.write(str(nPoints)+"\n")
    o.write("1\n")                                                          # Escrevendo a flag para indicar que as posições X estão dadas a partir do primeiro eletrodo  
    o.write("0\n")                                                          # Escrevendo a flag para indicar a ausência de prolaridade induzida no dado 

    cont = 0                                                                # Contador para indexar os os arrays     
    for levels in range(invesLevel):                                        # Loop varrendo os níveis de investigação
        for k,n in enumerate(index):                                        # Loop varrendo os indices na tabela de dados de campo
            if (n == "n" + str(levels+1)):                                  # Fixando níveis de investigação para construir arquivo texto
                x = displPointsPosition[cont]                               # Posição X de cada ponto
                r = apparentResistivity[k]                                  # Posição Z de cada ponto
                o.write(f"{x:.1f}   {electSpace}    {levels+1}  {r:.2f}\n") # Escrevendo em arquivo no formato necessário para a inversão
                orderedData = np.append(orderedData,float(r))               # Coletando o arquivo ordenado para o plot a seguir
                cont += 1                                                   # Adicionando valores ao contador para passar pro próximo ponto

    o.write("0,0,0,0,0,0,0,0\n")                                            # Escrevendo alguns zeros necessários para o arquivo de inversão
    o.close()                                                               # Fechando o arquivo texto     
#%% Imagem da geometria e, quando houver dados, os valores de rhoa 

plt.figure(1,figsize=(15,5))                                                # Criando uma figura com um tamanho pré-definido          
 
plt.scatter(electrodesPosition,np.zeros(nElectrodes))                       # Plotando somente os eletrodos 

X = displPointsPosition                                                     # Coletando as posições X dos dados
Z = depthPointsPosition                                                     # Coletando as posições Z dos dados
plt.scatter(X,Z)                                                            # Plotando os dados da aquisição configurada

if (not justSee):                                                           # Se existir tabela de campo 
    for i, values in enumerate(orderedData):                                # Coleto dos valores dos dados ordenados 
        plt.annotate(np.round(values,decimals=1),(X[i],Z[i]))               # Anoto no mapa os valores das resistividades aparentes coletadas em campo 

plt.gca().invert_yaxis()                                                    # Invertendo eixo X para a profundidade ser escrita para baixo no mapa

plt.title(f"Linha com {nPoints} pontos e {nElectrodes} eletrodos")          # Título da imagem 
 
plt.xlabel("Distância [m]",fontsize=15)                                     # Eixo X da imagem  
plt.ylabel("Profundidade [m]",fontsize=15)                                  # Eixo Z da imagem
plt.legend(["Eletrodos","Dados"],fontsize=12,loc="lower left")              # Legenda da imagem
plt.savefig("geometria.png")                                                # Salvando a figura em uma imagem .png
plt.show()                                                                  # Mostrando a imagem na tela   

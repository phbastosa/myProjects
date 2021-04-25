#%% Biblioteca utilizadas

import os
import plotly
import requests
import openpyxl

import numpy as np
import pandas as pd
import yfinance as yf

from bs4 import BeautifulSoup

#%% Valores a serem filtrados

agilMin = 0.8       # Preço por valor patrimonial mínimo
agilMax = 1.0       # Preço por valor patrimonial máximo

precMin = 80.0      # Preço mínimo da cota
precMax = 100.0     # Preço máximo da cota

dividendMin = 0.8   # Último dividendo mínimo pago

#%% Coletando os dados do site fundsExplorer

url = "https://www.fundsexplorer.com.br/ranking"

response = requests.get(url)

soup = BeautifulSoup(response.text, "lxml")

table = soup.find_all("table")[0]

dataFrame = pd.read_html(str(table), decimal=".", thousands="")[0]

#%% Mudando os títulos das colunas do dataFrame

parameters = dataFrame.columns
columns = [0,1,2,3,4,5,18,23,24,25]
names = ["Fundo", "Setor", "Preço", "Liquidez", "Dividendo", "DY", "P/VP", "Vac Física", "Vac Financeira", "Qtd Ativos"]

funds = dataFrame[parameters[columns]]

#%% Organizando os dados para serem analizados 

for i in range(len(columns)):
    funds = funds.rename({parameters[columns[i]]:names[i]}, axis=1)

for i in funds.index:
    if funds["Preço"].isnull()[i]:
        funds = funds.drop(i)

for i in funds.index:
    funds["Preço"][i] = funds["Preço"][i][2:].replace(".","")
    funds["Preço"][i] = float(funds["Preço"][i].replace(",","."))
    funds["P/VP"][i] = float(funds["P/VP"][i].replace(",","."))
    funds["Dividendo"][i] = float(funds["Dividendo"][i][2:].replace(",","."))

#%% Filtrando para coletar os fundos que se encaixam na filtragem

bestFunds = []

os.system("clear")
for i in funds.index:
    if agilMin < funds["P/VP"][i] < agilMax: 
        if precMin < funds["Preço"][i] < precMax:     
            if funds["Dividendo"][i] > dividendMin:
                bestFunds.append(i)

for i in funds.index:
    if i not in bestFunds:
        funds = funds.drop(i)

#%% Apresentação dos fundos na tela do terminal

print(funds)
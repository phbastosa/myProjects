import functions as fnc 
from time import sleep
from random import choice

posicoes = [1,2,3,4,5,6,7,8,9]
somatorio = [0,0,0,0,0,0,0,0,0]
jogadas = ['_','_','_','_','_','_','_','_','_']
print(("-=")*30)
print("Instruções:")
print(("-=")*30)
print("Os números mostrados abaixo são \nas posições possíveis para as jogadas.")
fnc.show_matrix(posicoes)
print("Jogador será marcado com X e computador com O.")
print(("-=")*20)
print("Situação atual:")

jog = fnc.player(jogadas,somatorio) # jogada do jogador
print(("-=")*30)
print("Vez do computador...")
print("Joguei aqui!")

if jog == 5:
    jogadas[0] = 'O'
    somatorio[0] = -1
    
    jog = fnc.player(jogadas,somatorio) # jogada do jogador
    print(("-=")*30)
    print("Vez do computador...")
    print("Joguei aqui!")

    b = fnc.block(somatorio)
    if b >= 0:
        b = fnc.block(somatorio)
        jogadas[b] = 'O'
        somatorio[b] = -1
    else:
        jogadas[6] = 'O'
        somatorio[6] = -1
    
    jog = fnc.player(jogadas,somatorio) # jogada do jogador
    print(("-=")*30)
    print("Vez do computador...")
    print("Joguei aqui!")

    a = fnc.attack(somatorio)
    b = fnc.block(somatorio)
    if a >= 0:
        jogadas[a] = 'O'
        somatorio[a] = -1
        fnc.show_matrix(jogadas)
        print(("-=")*30)
        fnc.winner(somatorio)
        exit()
    elif b >= 0:
        jogadas[b] = 'O'
        somatorio[b] = -1
    else:
        for i in range(len(somatorio)):
            possibilidade = []
            if somatorio[i] == 0:
                possibilidade.append(i)

        p = choice(possibilidade)
        jogadas[p] = 'O'
        somatorio[p] = -1
                
    jog = fnc.player(jogadas,somatorio) # jogada do jogador
    print(("-=")*30)
    print("Vez do computador...")
    print("Joguei aqui!")

    a = fnc.attack(somatorio)
    b = fnc.block(somatorio)
    if a >= 0:
        jogadas[a] = 'O'
        somatorio[a] = -1
        fnc.show_matrix(jogadas)
        print(("-=")*30)
        fnc.winner(somatorio)
        exit()
    elif b >= 0:
        jogadas[b] = 'O'
        somatorio[b] = -1
    else:
        possibilidade = []
        for i in range(len(somatorio)):
            if somatorio[i] == 0:
                possibilidade.append(i)

        p = choice(possibilidade)
        jogadas[p] = 'O'
        somatorio[p] = -1

    jog = fnc.player(jogadas,somatorio) # jogada do jogador
    print(("-=")*30)
    print("Velha!\n")    

else:
    jog_ant = jog
    if jog in [1,3,7,9]:
        jogadas[4] = 'O'
        somatorio[4] = -1
    else:
        if jog in [2,6,8]:
            jogadas[jog] = 'O'
            somatorio[jog] = -1
        else:
            jogadas[0] = 'O'
            somatorio[0] = -1
               
    jog = fnc.player(jogadas,somatorio) # jogada do jogador
    print(("-=")*30)
    print("Vez do computador...")
    print("Joguei aqui!")

    b = fnc.block(somatorio)
    if b >= 0:
        jogadas[b] = 'O'
        somatorio[b] = -1
    else:
        if jog_ant and jog in [1,3,7,9]:
            if somatorio[5] == 0:
                jogadas[5] = 'O'
                somatorio[5] = -1
            else:
                jogadas[2] = 'O'
                somatorio[2] = -1
                    
        if jog_ant in [1,3,7,9] and jog in [2,4,6,8]:
            if somatorio[3] == 0:
                jogadas[2] = 'O'
                somatorio[2] = -1
            elif somatorio[5] == 0:
                jogadas[0] = 'O'
                somatorio[0] = -1

        if jog_ant in [2,6,8] and jog in [1,3,7,9]:
            if somatorio[0] == 1:
                jogadas[6] = 'O'
                somatorio[6] = -1
            else:
                jogadas[8] = 'O'
                somatorio[8] = -1

    jog = fnc.player(jogadas,somatorio) # jogada do jogador
    print(("-=")*30)
    print("Vez do computador...")
    print("Joguei aqui!")

    a = fnc.attack(somatorio)
    b = fnc.block(somatorio)
    if a >= 0:
        jogadas[a] = 'O'
        somatorio[a] = -1
        fnc.show_matrix(jogadas)
        print(("-=")*30)
        fnc.winner(somatorio)
        exit()
    elif b >= 0:
        jogadas[b] = 'O'
        somatorio[b] = -1
    else:
        possibilidade = []
        for i in range(len(somatorio)):
            if somatorio[i] == 0:
                possibilidade.append(i)

        p = choice(possibilidade)
        jogadas[p] = 'O'
        somatorio[p] = -1

    jog = fnc.player(jogadas,somatorio) # jogada do jogador
    print(("-=")*30)
    print("Vez do computador...")
    print("Joguei aqui!")

    a = fnc.attack(somatorio)
    b = fnc.block(somatorio)
    if a >= 0:
        jogadas[a] = 'O'
        somatorio[a] = -1
        fnc.show_matrix(jogadas)
        print(("-=")*30)
        fnc.winner(somatorio)
        exit()
    elif b >= 0:
        jogadas[b] = 'O'
        somatorio[b] = -1
    else:
        possibilidade = []
        for i in range(len(somatorio)):
            if somatorio[i] == 0:
                possibilidade.append(i)

        p = choice(possibilidade)
        jogadas[p] = 'O'
        somatorio[p] = -1

    jog = fnc.player(jogadas,somatorio) # jogada do jogador
    print(("-=")*30)
    print("Velha!\n")    

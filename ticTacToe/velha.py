import functions as fnc 
from time import sleep
from random import choice

posicoes = [1,2,3,4,5,6,7,8,9]
print(("-=")*30)
print("Instruções:")
print(("-=")*30)
print("Os números mostrados abaixo são \nas posições possíveis para as jogadas.")
fnc.show_matrix(posicoes)
print("Jogador será marcado com X e computador com O.")
print(("-=")*30)

while True:
    somatorio = [0,0,0,0,0,0,0,0,0]
    jogadas = ['_','_','_','_','_','_','_','_','_']

    print("Situação atual:")
    jog = fnc.player(jogadas,somatorio) # jogada do jogador
    print(("-=")*30)
    print("Vez do computador...")
    sleep(1)
    print("Joguei aqui!")

    if jog == 5:
        jogadas[0] = 'O'
        somatorio[0] = -1
        
        jog = fnc.player(jogadas,somatorio) # jogada do jogador
        print(("-=")*30)
        print("Vez do computador...")
        sleep(1)
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
        sleep(1)
        print("Joguei aqui!")

        a = fnc.attack(somatorio)
        b = fnc.block(somatorio)
        if a >= 0:
            jogadas[a] = 'O'
            somatorio[a] = -1
            fnc.show_matrix(jogadas)
            print(("-=")*30)
            fnc.winner(somatorio)
            print(("-=")*30)

            print("Deseja jogar novamente? [s ou n]")
            r = str(input("Resposta: "))
            print(("-=")*30)

            if r == 'n':
                print("\nFoi muito bom jogar com você!!\n")
                exit()
            else:
                continue                    

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
        sleep(1)
        print("Joguei aqui!")

        a = fnc.attack(somatorio)
        b = fnc.block(somatorio)
        if a >= 0:
            jogadas[a] = 'O'
            somatorio[a] = -1
            fnc.show_matrix(jogadas)
            print(("-=")*30)
            fnc.winner(somatorio)
            print(("-=")*30)

            print("Deseja jogar novamente? [s ou n]")
            r = str(input("Resposta: "))
            print(("-=")*30)
            if r == 'n':
                print("\nFoi muito bom jogar com você!!\n")
                exit()
            else:
                continue                    

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
        print(("-=")*30)

        print("Deseja jogar novamente? [s ou n]")
        r = str(input("Resposta: "))
        print(("-=")*30)
        if r == 'n':
            print("\nFoi muito bom jogar com você!!\n")
            exit()
        else:
            continue                    

    else:  # Divisor de águas
    
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
        sleep(1)
        print("Joguei aqui!")

        b = fnc.block(somatorio)
        if b >= 0:
            jogadas[b] = 'O'
            somatorio[b] = -1
        else:                    
            if jog and jog_ant in [1,3,7,9]:
                if somatorio[1] == 0:
                    jogadas[1] = 'O'
                    somatorio[1] = -1
            
            if jog_ant in [2,4,6,8]:
                if jog in [1,3,7,9]:
                    jogadas[4] = 'O'
                    somatorio[4] = -1
                
            if jog_ant in [4,8]:
                if jog in [2,4,6,8]:
                    jogadas[4] = 'O'
                    somatorio[4] = -1

            if jog_ant == 2:
                if jog in [4,6]:
                    jogadas[4] = 'O'
                    somatorio[4] = -1

            if jog_ant == 6:
                if jog in [2,8]:
                    jogadas[0] = 'O'
                    somatorio[0] = -1


        jog = fnc.player(jogadas,somatorio) # jogada do jogador
        print(("-=")*30)
        print("Vez do computador...")
        sleep(1)
        print("Joguei aqui!")

        a = fnc.attack(somatorio)
        b = fnc.block(somatorio)
        if a >= 0:
            jogadas[a] = 'O'
            somatorio[a] = -1
            fnc.show_matrix(jogadas)
            print(("-=")*30)
            fnc.winner(somatorio)
            print(("-=")*30)

            print("Deseja jogar novamente? [s ou n]")
            r = str(input("Resposta: "))
            print(("-=")*30)
            if r == 'n':
                print("\nFoi muito bom jogar com você!!\n")
                exit()
            else:
                continue                    

        elif b >= 0:
            jogadas[b] = 'O'
            somatorio[b] = -1
        else:
            if jog == 6 and somatorio[6] == 0:
                jogadas[6] = 'O'
                somatorio[6] = -1    
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
        sleep(1)
        print("Joguei aqui!")

        a = fnc.attack(somatorio)
        b = fnc.block(somatorio)
        if a >= 0:
            jogadas[a] = 'O'
            somatorio[a] = -1
            fnc.show_matrix(jogadas)
            print(("-=")*30)
            fnc.winner(somatorio)
            print(("-=")*30)

            print("Deseja jogar novamente? [s ou n]")
            r = str(input("Resposta: "))
            print(("-=")*30)
            if r == 'n':
                print("\nFoi muito bom jogar com você!!\n")
                exit()
            else:
                continue                    

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

        print("Deseja jogar novamente? [s ou n]")
        r = str(input("Resposta: "))
        print(("-=")*30)
        if r == 'n':
            print("\nFoi muito bom jogar com você!!\n")
            exit()
        else:
            continue                    


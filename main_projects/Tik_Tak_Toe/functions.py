def show_matrix(matrix):
    print(f"\n           {matrix[0]} | {matrix[1]} | {matrix[2]}\n")
    print(f"           {matrix[3]} | {matrix[4]} | {matrix[5]}\n")
    print(f"           {matrix[6]} | {matrix[7]} | {matrix[8]}\n")

def player(jogada,somatorio):
    show_matrix(jogada)
    print(("-=")*30)        
    jog = int(input("Faça a sua jogada numa posição possível: "))
    while True:
        if jogada[jog-1] == '_':
            break
        else:
            print("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!")
            print("Jogada inválida, tente novamente em uma nova posição...")
            print("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n")
            show_matrix(jogada)        
            print(("-=")*30)
            jog = int(input("Faça a sua jogada numa posição possível: "))
    
    jogada[jog-1] = 'X'
    somatorio[jog-1] = 1
    show_matrix(jogada)

    return jog

def attack(s):
    lin1 = s[0] + s[1] + s[2]
    lin2 = s[3] + s[4] + s[5]
    lin3 = s[6] + s[7] + s[8]
    col1 = s[0] + s[3] + s[6]
    col2 = s[1] + s[4] + s[7]
    col3 = s[2] + s[5] + s[8] 
    diag1 = s[0] + s[4] + s[8]
    diag2 = s[2] + s[4] + s[6]

    pos = -1
    if lin1 == -2:
        if s[0] == 0:
            pos = 0
            return pos
        elif s[1] == 0:
            pos = 1
            return pos
        elif s[2] == 0:
            pos = 2
            return pos        

    if lin2 == -2:
        if s[3] == 0:
            pos = 3
            return pos
        elif s[4] == 0:
            pos = 4
            return pos
        elif s[5] == 0:
            pos = 5
            return pos        

    if lin3 == -2:
        if s[6] == 0:
            pos = 6
            return pos
        elif s[7] == 0:
            pos = 7
            return pos
        elif s[8] == 0:
            pos = 8
            return pos        

    if col1 == -2:
        if s[0] == 0:
            pos = 0
            return pos
        elif s[3] == 0:
            pos = 3
            return pos
        elif s[6] == 0:
            pos = 6
            return pos        

    if col2 == -2:
        if s[1] == 0:
            pos = 1
            return pos
        elif s[4] == 0:
            pos = 4
            return pos
        elif s[7] == 0:
            pos = 7
            return pos        

    if col3 == -2:
        if s[2] == 0:
            pos = 2 
            return pos
        elif s[5] == 0:
            pos = 5
            return pos
        elif s[8] == 0:
            pos = 8
            return pos        

    if diag1 == -2:
        if s[0] == 0:
            pos = 0
            return pos
        elif s[4] == 0:
            pos = 4
            return pos
        elif s[8] == 0:
            pos = 8
            return pos        

    if diag2 == -2:
        if s[2] == 0:
            pos = 2
            return pos
        elif s[4] == 0:
            pos = 4
            return pos 
        elif s[6] == 0:
            pos = 6
            return pos

    return pos                

def block(s):
    lin1 = s[0] + s[1] + s[2]
    lin2 = s[3] + s[4] + s[5]
    lin3 = s[6] + s[7] + s[8]
    col1 = s[0] + s[3] + s[6]
    col2 = s[1] + s[4] + s[7]
    col3 = s[2] + s[5] + s[8] 
    diag1 = s[0] + s[4] + s[8]
    diag2 = s[2] + s[4] + s[6]

    pos = -1
    if lin1 == 2:
        if s[0] == 0:
            pos = 0
            return pos
        elif s[1] == 0:
            pos = 1
            return pos
        elif s[2] == 0:
            pos = 2
            return pos        

    if lin2 == 2:
        if s[3] == 0:
            pos = 3
            return pos
        elif s[4] == 0:
            pos = 4
            return pos
        elif s[5] == 0:
            pos = 5
            return pos        

    if lin3 == 2:
        if s[6] == 0:
            pos = 6
            return pos
        elif s[7] == 0:
            pos = 7
            return pos
        elif s[8] == 0:
            pos = 8
            return pos        

    if col1 == 2:
        if s[0] == 0:
            pos = 0
            return pos
        elif s[3] == 0:
            pos = 3
            return pos
        elif s[6] == 0:
            pos = 6
            return pos        

    if col2 == 2:
        if s[1] == 0:
            pos = 1
            return pos
        elif s[4] == 0:
            pos = 4
            return pos
        elif s[7] == 0:
            pos = 7
            return pos        

    if col3 == 2:
        if s[2] == 0:
            pos = 2 
            return pos
        elif s[5] == 0:
            pos = 5
            return pos
        elif s[8] == 0:
            pos = 8
            return pos        

    if diag1 == 2:
        if s[0] == 0:
            pos = 0
            return pos
        elif s[4] == 0:
            pos = 4
            return pos
        elif s[8] == 0:
            pos = 8
            return pos        

    if diag2 == 2:
        if s[2] == 0:
            pos = 2
            return pos
        elif s[4] == 0:
            pos = 4
            return pos 
        elif s[6] == 0:
            pos = 6
            return pos 
    
    return pos               

def winner(s):
    lin1 = s[0] + s[1] + s[2]
    lin2 = s[3] + s[4] + s[5]
    lin3 = s[6] + s[7] + s[8]
    col1 = s[0] + s[3] + s[6]
    col2 = s[1] + s[4] + s[7]
    col3 = s[2] + s[5] + s[8] 
    diag1 = s[0] + s[4] + s[8]
    diag2 = s[2] + s[4] + s[6]

    if lin1 == 3 or lin2 == 3 or lin3 == 3 or col1 == 3 or col2 == 3 or col3 == 3 or diag1 == 3 or diag2 == 3:
        print("Jogador ganhou!\n")

    if lin1 == -3 or lin2 == -3 or lin3 == -3 or col1 == -3 or col2 == -3 or col3 == -3 or diag1 == -3 or diag2 == -3:
        print("Computador ganhou!\n")


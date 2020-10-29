# include <iostream>
# include <unistd.h>
# include <stdlib.h>
# include <stdio.h>
# include "SnakeUtilities.h"
# include "kbhitgetch.h"

using namespace std;

int main()
{
    char dir;
    bool start = true;
    int cont = 0;
    int xHead = 28;
    int yHead = 15;

    Grid * BG = new Grid();
    Point * SH = new Point(xHead,yHead);        
    Tail * TL = new Tail(SH);
    Fruit * FT = new Fruit(BG,TL);

    while (start)
    {
        if (cont == 0) {cout << "Press any key to start!" << endl; cont++;}

        dir = getch();   
        if (dir == 'x')
        {
            system("clear");
            cout << "Game Over!" << endl;
            cout << "Thank you for playing!!" << endl;   
            break;       
        }

        system("clear");
        
        BG-> refreshGrid();
        FT-> catchFruit(SH,TL,BG, SH-> movePoint(dir,SH));

        if (GameOver(BG,SH,TL)) 
        { 
            TL-> zeroingTail();
            FT-> setHighScore(FT-> getScore());
            FT-> zeroingScore();    
        } 

        Show(BG,SH,TL,FT);         
        sleep(0.01);
    }

    return 0;
}
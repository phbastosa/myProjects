# include <time.h>
# include <stdio.h>
# include <unistd.h>
# include <iostream>
# include "gameClass.h"

using namespace std; // comment

int main(int argc,char **argv)
{
    srand(time(NULL));

    Game * blackJack = new Game();

    blackJack->play();

    return 0;
}

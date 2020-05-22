# ifndef CLASSES_H_INCLUDED
# define CLASSES_H_INCLUDED

# include "essentialClasses.h"

class Game
{
    private:

        int pushs = 0;        /* Number of pushes */  
        int dealerWins = 0;   /* number of dealer wins */  
        int playerWins = 0;   /* Number of player wins */ 
    
        int playerPoints = 0; /* Player points at each round */
        int dealerPoints = 0; /* Dealer points at each round */ 
        int move = 99;        /* Parameter to choice a movement */

        bool doubleGame;      /* Flag to initiate a split game*/
        bool usingHand1;      /* Flag to indicate use of the first hand */
        int  flagDouble;      /* Flag to indicate how many times the second hand was used */
        bool blackJackf;      /* Flag to indicate if the player get a BlackJack */
        bool hand1bustf;      /* Flag to indicate uf hand 1 make a bust game */

        const char * message; /* String used to indicates a result */ 
        const char * hidden1; /* Message to hidden the type of the dealer's card hidden */ 
        const char * hidden2; /* Message to hidden the suit of the dealer's card hidden */ 
        const char * hidden3; /* Message to expose the BlackJack result */

        Deck * deck;          /* Deck object */ 
        Hand * dealer;        /* Dealer's hand object */ 
        Hand * playerN;       /* Player's first hand object */ 
        Hand * playerD;       /* Player's secundary hand object */ 

    public:

        /* Constructor method */
        Game()
        {
            this->restart(); /* This method initialize the atributes and the objects */
        }

        /* This method write in terminal the game interface */ 
        void display()
        {
            /* Condition to get a right points in screen */
            if (this->usingHand1) /* If the player using the first hand */
            {
                this->playerPoints = this->playerN->points(); /* Put on the screen the first hand points */
            }
            else
            {
                this->playerPoints = this->playerD->points(); /* Put on the screen the second hand points */   
            }
            
            /* The game screen */
            printf("----------------------|--------------------------------------------------------|\n");
            printf("      BlackJack       |                       Dealer's hand                    |\n");
            printf("----------------------| |----------|----------|----------|----------|----------|\n");
            printf("Statistics:           | |  Card 1  |  Card 2  |  Card 3  |  Card 4  |  Card 5  |\n");
            printf("                      | |          |          |          |          |          |\n");
            printf("Player wins: %2.0f       | |%s|%s|%s|%s|%s|\n",(float) this->playerWins,this->hidden1,this->dealer->getCardsType(1),this->dealer->getCardsType(2),this->dealer->getCardsType(3),this->dealer->getCardsType(4));
            printf("Dealer wins: %2.0f       | |%s|%s|%s|%s|%s|\n",(float) this->dealerWins,this->hidden2,this->dealer->getCardsSuit(1),this->dealer->getCardsSuit(2),this->dealer->getCardsSuit(3),this->dealer->getCardsSuit(4));
            printf("----------------------| |----------|----------|----------|----------|----------|\n");
            printf("Results:              |                                                        |\n");
            printf("                      |                       Player's hand                    |\n");
            printf("%s| |----------|----------|----------|----------|----------|\n",this->message);
            printf("%s| |  Card 1  |  Card 2  |  Card 3  |  Card 4  |  Card 5  |\n",this->hidden3);
            printf("----------------------| |          |          |          |          |          |\n");
            printf("Options:              | |%s|%s|%s|%s|%s|\n",this->playerN->getCardsType(0),this->playerN->getCardsType(1),this->playerN->getCardsType(2),this->playerN->getCardsType(3),this->playerN->getCardsType(4));
            printf("[1] Pull card         | |%s|%s|%s|%s|%s|\n",this->playerN->getCardsSuit(0),this->playerN->getCardsSuit(1),this->playerN->getCardsSuit(2),this->playerN->getCardsSuit(3),this->playerN->getCardsSuit(4));        
            printf("[2] Hold hand         | |----------|----------|----------|----------|----------|\n");        
            printf("[3] Double game       | |  Card 1  |  Card 2  |  Card 3  |  Card 4  |  Card 5  |\n");        
            printf("[4] Exit game         | |          |          |          |          |          |\n");
            printf("----------------------| |%s|%s|%s|%s|%s|\n",this->playerD->getCardsType(0),this->playerD->getCardsType(1),this->playerD->getCardsType(2),this->playerD->getCardsType(3),this->playerD->getCardsType(4));
            printf("Dealer's points: %2.0f   | |%s|%s|%s|%s|%s|\n",(float) this->dealerPoints,this->playerD->getCardsSuit(0),this->playerD->getCardsSuit(1),this->playerD->getCardsSuit(2),this->playerD->getCardsSuit(3),this->playerD->getCardsSuit(4));        
            printf("Player's points: %2.0f   | |----------|----------|----------|----------|----------|\n",(float) this->playerPoints);
            printf("-------------------------------------------------------------------------------|\n");               
        }

        /* The method to execute the game */
        void play()
        {
            while (true)  /* Infinite loop */
            {
                system("clear");             /* Clean the screen at each loop */

                this->display();             /* Show the game on terminal */   
    
                printf("Type your move: ");  /* Just a message to the user */
                scanf("%i",&this->move);     /* Collect the number typed */

                if (this->move == 4)         /* Flag to exit game */ 
                {
                    this->gameOverMessage(); /* Message to user when the game is closed */
                    break;                   /* Closing game */
                }

                this->action(this->move);    /* Method to carry out the collected instructions */
            }
        }

        /* This method makes the decision from the information entered by the user */
        void action(int arg)
        {   
            switch (arg) 
            {
                case 0: /* Only to debug the code making a restart to test possibilities */                     
                    this->restart();     
                    break;
                
                case 1:  /* Flag to pull a card of the deck */
                    if (this->doubleGame)  /* This condition separates two possibilities, a game with a single hand and a game with a two hands */ 
                    {    
                        if (this->usingHand1) /* This condition is to pull card for the first hand */
                        {
                            this->pullCard(this->playerN,this->playerD); /* Method to pull a single card */
                            this->isBust();                              /* Verifying if this hand is dead */

                            if(this->playerN->points() > 21)  /* This condition in case of Bust, moving the pull card method to a second hand */
                            {
                                this->usingHand1 = false; /* Finishing works with a first hand */
                            }    
                        }
                        else /* Condition to work with a second hand */
                        {
                            this->pullCard(this->playerD,this->playerN); /* Pulling card with a second hand */
                            this->isBust();                              /* Verifying if this hand is dead */

                            if (this->playerD->points() > 21)  /* Condition in case of Bust with the second hand */
                            {   
                                this->flagDouble++;                            /* Triggering flag to use in the next condition */
                                this->dealerTime(this->playerN,this->playerD); /* Making a dealer moves */ 
                            }
                        }
                    }
                    else  /* In case the game is happen with a single hand */
                    {
                        this->pullCard(this->playerN,this->playerD); /* Pulling card with the first hand (the only one)*/
                        this->isBust(); /* checking if the value has passed of 21 */
                    }            
                    break;    
                
                case 2:   /* Is the option to hold hand */      
                    if (this->doubleGame)   /* This condition determinates if the game is in single or two hands */
                    {
                        this->usingHand1 = false; /* In the first time when the user come here, he still with the first hand, so this flag activate the secoond one */
                        
                        if (this->flagDouble > 0) /* This condition decides make a dealer moves when the user hold the second hand */ 
                        {
                            this->dealerTime(this->playerN,this->playerD); /* Making a dealer moves */
                        }

                        this->flagDouble++; /* Adapting flag if the user holded the first hand */
                    }
                    else /* Simple hand case */
                    {
                        this->dealerTime(this->playerN,this->playerD); /* Just pulling and comparing the dealer's hand and the player's hand */   
                    }
                    break;    

                case 3:   /* It makes a split when the first type card of the first hand is equal to the second type card */
                    if (this->playerN->getCardsType(0) == this->playerN->getCardsType(1))  /* Comparing the types card and initializing parameters to do bouble hands game */
                    {
                        this->doubleGame = true;                                                 
                        this->playerD->setCards(0,this->playerN->getCardsSuit(1),this->playerN->getCardsType(1));
                        this->playerN->setCards(1,"    X     ","    X     ");
                        this->playerN->setTotalCards(1);
                        this->playerD->setTotalCards(1);
                    }
                    break;
            }
        }    

        /* Method to initializes all necessary atributes and objects each round */
        void restart()
        {
            // Necessary objects
            this->deck = new Deck();
            this->dealer = new Hand(deck);
            this->playerN = new Hand(deck);
            this->playerD = new Hand(deck);
            
            // Object seted to make a second hand 
            this->playerD->setXcards();
            this->playerD->setTotalCards(0);

            // Messages that change during the game 
            this->message = "   Game in progress   ";
            this->hidden1 = "  Hidden  ";
            this->hidden2 = "          ";
            this->hidden3 = "                      ";
            
            this->dealerPoints = (int) NULL;
            
            this->move = 0;
            
            this->doubleGame = false;
            this->usingHand1 = true;
            this->flagDouble = 0;
            this->blackJackf = false;
            this->hand1bustf = false;

            this->BlackJack(); /* Verifying if some hand has a blackjack */
        }

        /* Method to cacth a blackjack */
        void BlackJack()
        {
            if ((this->playerN->getTotalCards() == 2) && (this->dealer->getTotalCards() == 2)) /* It works only if dealer and player have a two cards in hand */
            {
                if ((this->dealer->points() == this->playerN->points()) && (this->dealer->points() == 21)) /* BlackJack Push */
                {
                    this->setVictory(3); /* flag to make a Push situation */ 
                }
                else if(this->dealer->points() == 21) /* Dealer BlackJack */
                {
                    this->hidden3 = "      BlackJack!!     ";  /* Hidden message exposed */
                    this->blackJackf = true;                   /* Activate flag to sum more points */
                    this->setVictory(0);                       /* Dealer wins */
                }
                else if (this->playerN->points() == 21) /* Player BlackJack */ 
                {
                    this->hidden3 = "      BlackJack!!     "; /* Hidden message exposed */
                    this->blackJackf = true;                  /* Activate flag to sum more points */
                    this->setVictory(1);                      /* Player wins */
                }
            }       
        }

        /* Method to verify Bust, when the points value goes out from 21*/
        void isBust()
        {
            if (this->doubleGame)
            {
                if(this->usingHand1)
                {
                    if(this->playerN->points() > 21)      /* Player's first hand Bust */
                    {
                        this->hand1bustf = true;
                    }
                }
                else
                {
                    if (this->playerD->points() > 21)    /* Player's second hand Bust */    
                    {
                        if (this->hand1bustf)            /* Two hand are bust*/
                        {
                            this->dealerPoints += 1;     /* Dealer wins, double points */
                            setVictory(0);              
                        }
                        else                            /* If second hand is bust but first isnt */
                        {
                            this->doubleGame = false;
                            this->dealerTime(this->playerN,this->playerD);  /* Dealer time with simple compare */   
                        }
                    }
                }
            }
            else   /* Checking the simple hand game */
            {
                if(this->playerN->points() > 21)      /* Player's first hand Bust */
                {
                    this->setVictory(0);
                }
            }            
        }

        /* Method to make dealer moves */
        void dealerTime(Hand * hand1, Hand * hand2)
        {            
            while (this->dealer->points() < 17)         /* Pulling cards until reach 17 points*/ 
            {   
                if(this->dealer->getTotalCards() < 5)   /* Until get the max cards in hand */
                {
                    this->dealer->pullCard(this->deck); /* Pulling the card */                
                }
                else
                {
                    break;                              /* Going out of the loop */
                }
            }

            if (this->dealer->points() >= 17)           /* Checking dealer points */ 
            {
                if (this->doubleGame)                   /* If the game is a split game */
                {
                    this->compareTwice(this->dealer,hand1,hand2); /* Comparing the points values of all of hands*/      
                }
                else                                    /* If the game is a single hand game */
                {
                    compare(this->dealer,hand1);        /* Comparing only the dealer and the player first hand */        
                }
            }
        }

        /* This method compares the first user's hand and the dealer's hand*/
        void compare(Hand * dealer, Hand * player)
        {
            if(dealer->points() <= 21)  /* If the dealer had a valid points */
            {    
                if(dealer->points() > player->points())        /* Dealer wins */
                {
                    this->setVictory(0);
                }
                else if (dealer->points() < player->points())  /* Player wins */
                {
                    this->setVictory(1);
                }
                else if (dealer->points() == player->points())
                {
                    this->setVictory(2);                       /* Push */
                }
            }
            else
            {
                this->setVictory(1);                  /* It dealer Bust */
            }
        }

        /* Method to compare dealer's hand with the two player's hand */
        void compareTwice(Hand * dealer, Hand * hand1, Hand * hand2)
        {
            if ((dealer->points() <= 21) && (hand1->points() <= 21) && (hand2->points() <= 21)) /* Almost the same to the previous method */
            {
                if ((dealer->points() == hand1->points()) && (dealer->points() == hand2->points()))
                {
                    this->setVictory(2); /* Push */ 
                }
                else if ((dealer->points() < hand1->points()) && (dealer->points() == hand2->points()))
                {
                    this->setVictory(1); /* Players wins, single point */
                }
                else if ((dealer->points() == hand1->points()) && (dealer->points() < hand2->points()))
                {
                    this->setVictory(1); /* Players wins, single point */
                }
                else if ((dealer->points() == hand1->points()) && (dealer->points() > hand2->points()))
                {
                    this->setVictory(0); /* Dealer wins, single point */
                }
                else if ((dealer->points() > hand1->points()) && (dealer->points() == hand2->points()))
                {
                    this->setVictory(0); /* Dealer wins, single point */
                }
                else if ((dealer->points() < hand1->points()) && (dealer->points() < hand2->points()))
                {
                    this->playerPoints += 1;
                    this->setVictory(1);     /* Player wins, double points */
                }
                else if ((dealer->points() > hand1->points()) && (dealer->points() > hand2->points()))
                {
                    this->dealerPoints += 1;
                    this->setVictory(0);     /* Dealer wins, double points */
                }
                else if ((dealer->points() < hand1->points()) && (dealer->points() > hand2->points()))
                {
                    this->setVictory(1);   /* Player wins, single points */
                }
                else if ((dealer->points() > hand1->points()) && (dealer->points() < hand2->points()))
                {
                    this->setVictory(1);   /* Player wins, single points */
                }
            }
            else /* If dealer bust */
            {
                this->playerPoints += 1;
                this->setVictory(1);     /* Player wins, double points */
            }
        }

        /* Method to change the messages on the screen when player stop playing */    
        void changeMessages()
        {
            this->dealerPoints = this->dealer->points();
            this->hidden1 = this->dealer->getCardsType(0);
            this->hidden2 = this->dealer->getCardsSuit(0);
        }

        void setVictory(int arg)
        {
            switch (arg)
            {
            case 0:
                this->message = "      Dealer wins     ";

                if (this->blackJackf) 
                {
                    this->dealerWins += 2;
                }
                else
                {
                    this->dealerWins += 1;
                }
                
                this->restartMessage();
                break;
            
            case 1:
                this->message = "      Player wins     ";
                
                if (this->blackJackf)
                {
                    this->playerWins += 2;
                }
                else
                {    
                    this->playerWins += 1;                
                }
                
                this->restartMessage();
                break;

            case 2:
                this->message = "         Push         ";
                this->pushs += 1;
                this->restartMessage();
                break;
            }
        }

        void pullCard(Hand * hand1, Hand * hand2)
        {
            if(hand1->getTotalCards() < 5)
            {
                hand1->pullCard(this->deck);
            }
            else
            {
                this->dealerTime(hand1,hand2);    
            }
        }

        void restartMessage()
        {
            this->changeMessages();
            this->display();
            printf("A new game will start automatically in a moment...\n");
            sleep(3);
            system("clear");
            this->restart();
        }

        void gameOverMessage()
        {
            printf("\n################################################################################\n");
            printf("                            Thank you for playing!!\n");
            printf("################################################################################\n");
        }
};

# endif
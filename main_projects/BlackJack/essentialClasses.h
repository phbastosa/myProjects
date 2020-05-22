# ifndef ESSENTIALCLASSES_H_INCLUDED
# define ESSENRIALCLASSES_H_INCLUDED

class Card
{
    private:
        const char * suit;
        const char * type;

    public:
        // Método contrutor: carta tem q receber naipe e tipo quando criada
        Card(const char * s,const char * t)
        {
            this->setSuit(s);
            this->setType(t);        
        }

        void setSuit(const char * s) 
        {
            this-> suit = s;
        }
        
        void setType(const char * t) 
        {
            this-> type = t;
        }

        const char * getSuit() 
        {
            return this-> suit;
        }
        
        const char * getType() 
        {
            return this-> type;
        }
};

class Deck
{
    private:
        static const int nsuits = 4;
        static const int ntypes = 12; 
        const char * suits[nsuits] = {"  Hearts  "," Diamonds ","  Spades  ","  Clubs   "};
        const char * types[ntypes] = {"    A     ","    2     ","    3     ","    4     ",
                                      "    5     ","    6     ","    7     ","    8     ",
                                      "    9     ","    J     ","    Q     ","    K     "};

    public:    
        // Método contrutor
        Deck(){}

        // Retira um tipo de carta aleatoriamente
        const char * randTypes() 
        {
            return this->types[rand()%this->getNT()];
        } 
        
        // Retira um naipe aleatoriamente 
        const char * randSuits() 
        {
            return this->suits[rand()%this->getNS()];
        }

        // Coleta o número máximo de naipes
        int getNS() 
        {
            return this-> nsuits;
        }    

        // Coleta o número máximo dos tipos de cartas
        int getNT() 
        {
            return this-> ntypes;
        }
};

class Hand
{
    private:
        static const int maxCards = 5; // Máximo de cartas na mão   
        Card * Cards[maxCards];        // Mão (conjunto de cartas)
        int totalCards;                // Cartas presentes na mão

    public:
    
    // Método construtor
    Hand(Deck * deck)
    {
        this->setXcards(); // Ocultando as cartas que não estão sendo utilizadas
        
        this->Cards[0] = new Card(deck->randSuits(),deck->randTypes()); // nova carta 1 aleatória
        this->Cards[1] = new Card(deck->randSuits(),deck->randTypes()); // nova carta 2 aleatória
        
        this->setTotalCards(2); // atualização da quantidade de cartas na mão
    }

    // Coleto o naipe da carta 
    const char * getCardsSuit(int n) 
    {
        return this->Cards[n]-> getSuit();
    }
    
    // Coleto o tipo da carta 
    const char * getCardsType(int n) 
    {
        return this->Cards[n]-> getType();
    }
    
    // Atribuindo qualquer valor para uma carta 
    void setCards(int position,const char * suit, const char * type)
    {
        this->Cards[position] = new Card(suit,type);
    }

    // Inicializando cartas não utilizadas
    void setXcards()
    {
        for (int i = 0; i < this->getMaxCards(); i++)
        {
            this->setCards(i,"    X     ","    X     ");
        }            
    }

    // Atribuindo carta à mão em formato de contador
    void setTotalCards(int num) 
    {
        this->totalCards = num;
    }
    
    // Coletando a quantidade de cartas na mão
    int getTotalCards() 
    {
        return this->totalCards;
    }

    // Coletando a quantidade máxima de cartas possiveis
    int getMaxCards() 
    {
        return this->maxCards;
    }

    // Puxando uma carta do baralho
    void pullCard(Deck * deck)
    {
        // se a quantidade de cartas na mão não ultrapassou o máximo, puxe a carta
        if (this->getTotalCards() < this->getMaxCards())
        {
            this->Cards[this->getTotalCards()] = new Card(deck->randSuits(),deck->randTypes());
            this->setTotalCards(this->getTotalCards()+1); // atribuindo mais uma carta no contador
        }
        else
        {
            std::cout<<"Já retirou o máximo de cartas!"<<std::endl;
        }
    }

    // Contando pontos numa mão qualquer
    int points()
    {
        int p = 0;
        static const int st = 8;  // quantidade de cartas simples
        static const int tt = 3;  // quantidade de cartas com valor 10 fixo
        static const int at = 1;  // quantidade de As (em um único naipe)

        // salvando os simbolos das cartas para o calculo
        const char * simpleTypes[st] = {"    2     ","    3     ","    4     ","    5     ","    6     ","    7     ","    8     ","    9     "};
        const char * tenTypes[tt] = {"    J     ","    Q     ","    K     "};
        const char * asType[at] = {"    A     "};

        for (int c = 0; c < this->getTotalCards(); c++)
        {
            // Atribuindo pontos de cartas comuns
            for (int simple = 0; simple < st; simple++)
            {
                if (this->getCardsType(c) == simpleTypes[simple])
                {
                    p += simple + 2;
                }
            }

            // Atribuindo pontos de cartas com valor 10 fixo
            for (int ten = 0; ten < tt; ten++)
            {
                if (this->getCardsType(c) == tenTypes[ten])
                {
                    p += 10;
                }
            }
        }

        for (int c = 0; c < this->getTotalCards(); c++)
        {        
            // Atribuindo valor de As independente da posição
            if (this->getCardsType(c) == asType[0])
            {
                if (p <= 10) p += 11;
        
                else if (p > 10) p += 1;
            }
        }

        return p;
    }
};

# endif
# ifndef SNAKEUTILITIES_H_INCLUDED
# define SNAKEUTILITIES_H_INCLUDED

class Grid
{
    private:
        int width = 60;
        int height = 20;

    public:
        char * B;

        Grid()
        {   
            this-> B = (char *) malloc(this-> getWidth()* this-> getHeight()*sizeof(char));
        }    

        void refreshGrid()
        {
            int index,i,j;

            for(index = 0; index < this-> getWidth()*this-> getHeight(); index++)
                this-> B[index] = (char) NULL;

            for(i = 0; i < this-> getHeight(); i++)
            {
                this-> B[i * this-> getWidth()] = '#';
                this-> B[i * this-> getWidth() + (this-> getWidth()-1)] = '#';   
                
                for(j = 1; j < this-> getWidth()-1; j++)
                {
                    this-> B[0*this-> getWidth() + j] = '#';
                    this-> B[(this-> getHeight()-1)*this-> getWidth() + j] = '#'; 
                }
            }
        }

        void setPoint(int x, int y, char parameter)
        {
            this-> B[y * this-> getWidth() + x] = parameter;
        }

        int getWidth()   {return this-> width; }
        int getHeight() {return this-> height;} 
        void setWidth(int w)  {this-> width = w; }
        void setHeight(int h) {this-> height = h;}
};

class Point
{
    private:
        int x;
        int y;

    public:
        Point(int x, int y)
        {
            this-> setX(x);
            this-> setY(y);        
        }

        Point * movePoint(char dir, Point * SH)
        {
            Point * p = new Point(SH-> getX(),SH->getY()); 
           
            if (dir == 'w') {this-> setY(this-> y -= 1);}       
            if (dir == 'a') {this-> setX(this-> x -= 1);}
            if (dir == 'd') {this-> setX(this-> x += 1);}
            if (dir == 's') {this-> setY(this-> y += 1);}
        
            return p;
        }

        void setX(int x) { this-> x = x;}
        void setY(int y) { this-> y = y;}
        int getX() {return this-> x;}
        int getY() {return this-> y;}
};

typedef struct
{
    int x;
    int y;

} Position;

class Tail
{
    private:
        int nTail = 1000;
        int SizeTail = 4;    
        int x_i; 
        int y_i; 

    public:
        Position * pos = (Position *) malloc(nTail*sizeof(Position));

        Tail(Point * SH)
        {     
            this-> setXi(SH-> getX());
            this-> setYi(SH-> getY() + 1);

            for (int i = 0; i < this-> SizeTail; i++)
            {
                this-> pos[i].x = this-> getXi();
                this-> pos[i].y = this-> getYi() + i;
            }
        }

        void followHead(Point * SH, Point * Neck)
        {         
            for (int i = this-> getSizeTail()-1; i > 0 ; i--)
            {
                this-> pos[i].x = this-> pos[i-1].x;
                this-> pos[i].y = this-> pos[i-1].y;
            }       

            this-> pos[0].x = Neck-> getX();
            this-> pos[0].y = Neck-> getY();
        }

        void addTail()
        {
            this-> setSizeTail(this-> getSizeTail() + 1);
        }

        void zeroingTail()
        {
            this-> setSizeTail(0);
            for (int i = 0; i < this-> nTail; i++)
            {
                this->pos[i].x = (int) NULL;
                this->pos[i].y = (int) NULL;
            }
        }

        int getSizeTail()  {return this-> SizeTail;}
        void setSizeTail(int s) {this-> SizeTail = s;}
        void setXi(int x) {this-> x_i = x;}
        void setYi(int y) {this-> y_i = y;}
        int getXi() {return this-> x_i;}
        int getYi() {return this-> y_i;}
};

class Fruit
{
    private:
        int x_fruit;
        int y_fruit;
        int score = 0;
        int highScore = 0;
    public:
        Fruit(Grid * BG, Tail * TL)
        {
            this-> setFruitPosition(BG,TL);        
        }

        void setFruitPosition(Grid * BG, Tail * TL)
        {   
            int x = (int) rand() % (BG-> getWidth()-4) + 2;
            int y = (int) rand() % (BG-> getHeight()-4) + 2;

            for (int i = 0; i < TL-> getSizeTail(); i++)
            {
                if((TL->pos[i].x == x) && (TL->pos[i].y == y))
                {
                    x = (int) rand() % (BG-> getWidth()-4) + 2;
                    y = (int) rand() % (BG-> getHeight()-4) + 2;
                }
            }
            

            this-> x_fruit = x; 
            this-> y_fruit = y;
        }

        void catchFruit(Point * SH, Tail * TL, Grid * BG, Point * Neck)
        {
            if ((SH-> getX() == this-> getFruitX()) && (SH-> getY() == this-> getFruitY()))
            {
                TL-> addTail();
                this-> setScore(100);

                for (int i = TL-> getSizeTail()-1; i > 0 ; i--)
                {
                    TL-> pos[i].x = TL-> pos[i-1].x;
                    TL-> pos[i].y = TL-> pos[i-1].y;
                }       

                TL-> pos[0].x = Neck->getX();
                TL-> pos[0].y = Neck->getY();
                
                this-> setFruitPosition(BG,TL);
            } 
            else 
            {
                TL-> followHead(SH,Neck);
            }
        }

        int getFruitX() {return this-> x_fruit;}
        int getFruitY() {return this-> y_fruit;}
        int getScore() {return this-> score;}
        void setScore(int s) {this-> score += s;}
        void zeroingScore() {this-> score = 0;}
        int getHighScore() {return this-> highScore;}
        void setHighScore(int hs) {this-> highScore = hs;}
};

bool GameOver(Grid * BG, Point * SH, Tail * TL)
{    
    if (BG-> B[SH-> getY() * BG->getWidth() + SH-> getX()] == '#')
    {
        return true;
    }
    
    for (int i = 0; i < TL-> getSizeTail(); i++)
    {
        if ((TL-> pos[i].x == SH->getX()) && (TL-> pos[i].y == SH->getY()))
        {
            return true;
        }
    }
}

void Show(Grid * BG, Point * SH, Tail * TL, Fruit * FT)
{
    int i, j;

    for (i = 0; i < TL-> getSizeTail(); i++)
    {
        BG-> setPoint(TL-> pos[i].x, TL-> pos[i].y,'o');
    }

    BG-> setPoint(SH-> getX(), SH-> getY(),'O');

    BG-> setPoint(FT-> getFruitX(), FT-> getFruitY(),'%');

    for(i = 0; i < BG-> getHeight(); i++)
    {
        for(j = 0; j < BG-> getWidth(); j++)
        {
            if (BG->B[i * BG-> getWidth() + j] == (char) NULL)
                std::cout << " ";
            else
                std::cout << BG->B[i * BG-> getWidth() + j];
        }

        std::cout << std::endl;
    }
    
    std::cout << "Score: "<< FT-> getScore() << std::endl;
    std::cout << "Highscore: "<<FT-> getHighScore() << std::endl;
    std::cout << "Press x to exit..."<< std::endl;
};

# endif
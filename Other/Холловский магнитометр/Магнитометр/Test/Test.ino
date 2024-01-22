#include <GyverButton.h>

int i = 0;
bool ch[16][4] = {{0, 0, 0, 0},
{0, 0, 0, 1},
{0, 0, 1, 0},
{0, 0, 1, 1},
{0, 1, 0, 0},
{0, 1, 0, 1},
{0, 1, 1, 0},
{0, 1, 1, 1},
{1, 0, 0, 0},
{1, 0, 0, 1},
{1, 0, 1, 0},
{1, 0, 1, 1},
{1, 1, 0, 0},
{1, 1, 0, 1},
{1, 1, 1, 0},
{1, 1, 1, 1}};


GButton butt1(9);
GButton butt2(8);

void setup(void)
{
  Serial.begin(9600);
  butt2.setTimeout(1000);
}

void loop(void)
{
  butt1.tick();
  butt2.tick(); 
  if (butt1.isClick()) 
  {
    i++;
    if (i == 16) i = 0;
    Serial.print("Chanel: ");Serial.print(i);Serial.print("; ");Serial.print(ch[i][0]);Serial.print("|");Serial.print(ch[i][1]);Serial.print("|");Serial.print(ch[i][2]);Serial.print("|");Serial.println(ch[i][3]);         // проверка на один клик
  }
  
  if (butt2.isHolded())
  {
    Serial.println("Set zero");         // проверка на один клик
  }
}

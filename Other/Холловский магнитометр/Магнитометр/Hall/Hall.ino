#include <GyverOLED.h>

#define PIN A1

GyverOLED<SSD1306_128x32, OLED_BUFFER> oled;

float U, B;
float last_t = 0;
float zero = 2.4917;
float sens = 1.4; // мВ/Гаусс 

void setup() {
  Serial.begin(9600);
  oled.init();        // инициализация
  oled.clear();       // очистка
}

void loop() {
if (millis() - last_t > 300){
    U = analogRead(PIN);
    U = U*(5.0/1023.0);
    Serial.println(U, 3);
    B = (U-zero)/(sens/1000);
    printOled();
    last_t = millis();
  }
}

void printOled(){        //---------------------------------------------------------------- Print
  oled.clear();       // очистка  
  oled.setScale(1);
  oled.setCursorXY(10,15);
  oled.print("B = "); oled.print(B,0); oled.print(" Гаусс");  
  oled.update();
}

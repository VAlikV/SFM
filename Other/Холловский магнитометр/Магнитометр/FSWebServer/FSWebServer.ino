#include <WiFi.h>
#include <WiFiClient.h>
#include <WebServer.h>
#include <FS.h> 
#include <SPIFFS.h>               

IPAddress apIP(192, 168, 4, 1);

// Web интерфейс для устройства
WebServer HTTP(80);
// Для файловой системы
File fsUploadFile;

// Определяем переменные wifi
String _ssidAP = "WiFi";   // SSID AP точки доступа
String _passwordAP = "88888888"; // пароль точки доступа

void setup() {
  Serial.begin(115200);
  Serial.println("");
  //Запускаем файловую систему 
  Serial.println("Start FS");
  FS_init();
  Serial.println("Start WIFI");
  //Запускаем WIFI
  StartAPMode();
  //Настраиваем и запускаем HTTP интерфейс
  Serial.println("Start WebServer");
  HTTP_init();

}

void loop() {
  HTTP.handleClient();
  delay(1);
}

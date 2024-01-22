// подключаем библиотеку для WiFi-связи:
#include <WiFi.h>

// здесь пишем учетные данные своей сети:
const char* ssid     = "Hall";
const char* password = "11111111";
IPAddress apIP(192, 168, 3, 2);

// задаем номер порта для веб-сервера («80»):
WiFiServer server(80);

// переменная для хранения HTTP-запроса:
String header;

// вспомогательные переменные
// для хранения текущего состояния выходных контактов:
String output26State = "off";
String output27State = "off";

// задаем номера для выходных GPIO-контактов:
const int output26 = 26;
const int output27 = 27;

void setup() {
  Serial.begin(115200);
  // делаем эти GPIO-контакты выходными:
  pinMode(output26, OUTPUT);
  pinMode(output27, OUTPUT);
  // присваиваем им значение «LOW»:
  digitalWrite(output26, LOW);
  digitalWrite(output27, LOW);

  // подключаемся к WiFi-сети при помощи заданных выше SSID и пароля:
  
  WiFi.mode(WIFI_AP);
  WiFi.softAPConfig(apIP, apIP, IPAddress(255, 255, 255, 0));
  WiFi.softAP("Hall", "11111111");
  Serial.println("");
  Serial.println("WiFi up AP");
  
  // печатаем в мониторе порта локальный IP-адрес
  // и запускаем веб-сервер:
  Serial.println("");
  Serial.println("WiFi connected.");  //  "WiFi подключен."
  Serial.println("IP address: ");  //  "IP-адрес: "
  Serial.println(WiFi.localIP());
  server.begin();
}

void loop(){
  // начинаем прослушивать входящих клиентов:
  WiFiClient client = server.available();

  if (client) {                     // если подключился новый клиент,     
    Serial.println("New Client.");  // печатаем сообщение
                                    // «Новый клиент.»
                                    // в мониторе порта;
    String currentLine = "";        // создаем строку для хранения
                                    // входящих данных от клиента;
    while (client.connected()) {    // цикл while() будет работать
                                    // все то время, пока клиент
                                    // будет подключен к серверу;
      if (client.available()) {     // если у клиента есть данные,
                                    // которые можно прочесть, 
        char c = client.read();     // считываем байт, а затем    
        Serial.write(c);            // печатаем его в мониторе порта 
        header += c;
        if (c == '\n') {            // если этим байтом является
                                    // символ новой строки
          // если мы получим два символа новой строки подряд 
          // то это значит, что текущая строчка пуста;
          // это конец HTTP-запроса клиента,
          // а значит – пора отправлять ответ:
          if (currentLine.length() == 0) {
            // HTTP-заголовки всегда начинаются
            // с кода ответа (например, «HTTP/1.1 200 OK»)
            // и информации о типе контента
            // (чтобы клиент понимал, что получает);
            // в конце пишем пустую строчку:
            client.println("HTTP/1.1 200 OK");
            client.println("Content-type:text/html");
            client.println("Connection: close");
                       //  "Соединение: отключено"
            client.println();
            
            // с помощью этого кода
            // включаем и выключаем GPIO-контакты:
            if (header.indexOf("GET /26/on") >= 0) {
              Serial.println("GPIO 26 on");  //  "GPIO26 включен"
              output26State = "on";
              digitalWrite(output26, HIGH);
            } else if (header.indexOf("GET /26/off") >= 0) {
              Serial.println("GPIO 26 off");  //  "GPIO26 выключен"
              output26State = "off";
              digitalWrite(output26, LOW);
            } else if (header.indexOf("GET /27/on") >= 0) {
              Serial.println("GPIO 27 on");  //  "GPIO27 включен"
              output27State = "on";
              digitalWrite(output27, HIGH);
            } else if (header.indexOf("GET /27/off") >= 0) {
              Serial.println("GPIO 27 off");  //  "GPIO27 выключен"
              output27State = "off";
              digitalWrite(output27, LOW);
            }
            
            // показываем веб-страницу с помощью этого HTML-кода:
            client.println("<!DOCTYPE html><html>");
            client.println("<head><meta name=\"viewport\" content=\"width=device-width, initial-scale=1\">");
            client.println("<link rel=\"icon\" href=\"data:,\">");
            // с помощью CSS задаем стиль кнопок «ON» и «OFF»;
            // если есть желание, можно поэкспериментировать
            // с фоновым цветом и атрибутами размера шрифта:
            client.println("<style>html { font-family: Helvetica; display: inline-block; margin: 0px auto; text-align: center;}");
            client.println(".button { background-color: #4CAF50; border: none; color: white; padding: 16px 40px;");
            client.println("text-decoration: none; font-size: 30px; margin: 2px; cursor: pointer;}");
            client.println(".button2 {background-color: #555555;}</style></head>");
            
            // заголовок веб-страницы:
            client.println("<body><h1>ESP32 Web Server</h1>");
            
            // рисуем кнопку для контакта GPIO26
            // и показываем его текущее состояние (ON/OFF): 
            client.println("<p>GPIO 26 - State " + output26State + "</p>");
            // если на контакте «output26State» значение «off»,
            // показываем кнопку «ON»:
            if (output26State=="off") {
              client.println("<p><a href=\"/26/on\"><button class=\"button\">ON</button></a></p>");
            } else {
              client.println("<p><a href=\"/26/off\"><button class=\"button button2\">OFF</button></a></p>");
            } 
               
            // рисуем кнопку для контакта GPIO27
            // и показываем его текущее состояние (ON/OFF): 
            client.println("<p>GPIO 27 - State " + output27State + "</p>");
            // если на контакте «output27State» значение «off»,
            // показываем кнопку «ON»:
            if (output27State=="off") {
              client.println("<p><a href=\"/27/on\"><button class=\"button\">ON</button></a></p>");
            } else {
              client.println("<p><a href=\"/27/off\"><button class=\"button button2\">OFF</button></a></p>");
            }
            client.println("</body></html>");
            
            // конец HTTP-ответа задается 
            // с помощью дополнительной пустой строки:
            client.println();
            // выходим из цикла while:
            break;
          } else {  // если получили символ новой строки,
                    // очищаем текущую строку «currentLine»:
            currentLine = "";
          }
        } else if (c != '\r') {  // если получили любые данные,
                                 // кроме символа возврата каретки,
          currentLine += c;      // добавляем эти данные 
                                 // в конец строки «currentLine»
        }
      }
    }
    // очищаем переменную «header»:
    header = "";
    // отключаем соединение:
    client.stop();
    Serial.println("Client disconnected.");
               //  "Клиент отключен."
    Serial.println("");
  }
}

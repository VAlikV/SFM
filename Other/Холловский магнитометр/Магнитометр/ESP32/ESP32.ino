#include <WiFi.h>

float zero[16] = {2.4964, 2.4435, 2.4435, 2.4435, 2.4435, 2.4435, 2.4435, 2.4435, 2.4435, 2.4435, 2.4435, 2.4435, 2.4435, 2.4435, 2.4435, 2.4964};
float sens[16] = {3.118, 3.418, 3.418, 3.418, 3.418, 3.418, 3.418, 3.418, 3.418, 3.418, 3.418, 3.418, 3.418, 3.418, 3.418, 3.118}; // мВ/Гаусс 

String _ssidAP = "Hall";   
String _passAP = "11111111"; 
IPAddress apIP(192, 168, 3, 2);

WiFiServer server(80);

// переменная для хранения HTTP-запроса:
String header;

void setup(void)
{  
  Serial.begin(115200);
  StartAPMode();
  server.begin();
}

void loop(void)
{
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
            
            client.println("<p>GPIO 26 - State </p>");
 
            client.println("<p>GPIO 27 - State </p>");

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

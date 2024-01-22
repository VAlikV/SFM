void HTTP_init(void) {
  HTTP.on("/restart", handle_Restart); // Перезагрузка модуля по запросу вида http://192.168.0.101/restart?device=ok
  HTTP.on("/Led", handle_sDiode); // Перезагрузка модуля по запросу вида http://192.168.0.101/restart?device=ok
  HTTP.on("/configs.json", handle_ConfigJSON); // формирование configs.json страницы для передачи данных в web интерфейс
  HTTP.on("/Zero", handle_setZero);
  HTTP.on("/Refresh", handle_Refresh);
  // Запускаем HTTP сервер
  HTTP.begin();

}

// Перезагрузка модуля по запросу вида http://192.168.0.101/restart?device=ok
void handle_Restart() {
  String restart = HTTP.arg("device");
  if (restart == "ok") ESP.restart();
  HTTP.send(200, "text/plain", "OK");
}

void handle_sDiode(){
  String st = HTTP.arg("status");
  if (st == "1") {
    digitalWrite(LED_BUILTIN, HIGH);
    HTTP.send(200, "text/plain", "______________On");
  }
  if (st == "0") {
    digitalWrite(LED_BUILTIN, LOW);
    HTTP.send(200, "text/plain", "________________Off");
  }
}

void handle_setZero(){
  String st = HTTP.arg("number");
  Serial.println(st);
  int num = st.substring(0, st.indexOf("=")).toInt();
  zero[num] = st.substring(st.indexOf("=")+1, st.length()).toFloat()/100000;
  Serial.println(zero[num],5);
  HTTP.send(200, "text/plain", "Ok");
}

void handle_ConfigJSON() {
  String json = "{";  // Формировать строку для отправки в браузер json формат
  //{"Zero0":"3.1111","Zero1":"3.2222",..."Sens0":"1.1111","Sens1":"1.2222"}
  // Имя SSDP
  json += "\"Number\":\"121\"";
  for(int i=0; i<121; i++){
    json += ",\"Zero";
    json += i;
    json += "\":\"";
    json += zero[i]*100000;
    json += "\"";
  }
 for(int i=0; i<121; i++){
    json += ",\"Sens";
    json += i;
    json += "\":\"";
    json += sens[i]*100000;
    json += "\"";
  }
  json += "}";
  HTTP.send(200, "text/json", json);
}

void handle_Refresh() {
  String Resp = String(volts, DEC);
  HTTP.send(200, "text/plain", Resp);
}

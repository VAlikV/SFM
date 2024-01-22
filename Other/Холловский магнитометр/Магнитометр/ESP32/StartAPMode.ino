bool StartAPMode()
{
  WiFi.disconnect();
  WiFi.mode(WIFI_AP);
  WiFi.softAPConfig(apIP, apIP, IPAddress(255, 255, 255, 0));
  WiFi.softAP("Hall", "11111111");
  Serial.println("");
  Serial.println("WiFi up AP");
  return true;
}

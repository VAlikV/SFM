function val(id){
 var v = document.getElementById(id).value;
 return v;
}

function createXmlHttpObject(){
  if(window.XMLHttpRequest){
    xmlHttp = new XMLHttpRequest();
  }else{
    xmlHttp = new ActiveXObject('Microsoft.XMLHTTP');
  }
  return xmlHttp;
}

function load(){

  var xmlHttp=createXmlHttpObject();

  if(xmlHttp.readyState==0 || xmlHttp.readyState==4){
    xmlHttp.open('GET','/configs.json',true);
    //xmlHttp.send(null);
    xmlHttp.responseType = "json";
    xmlHttp.onload = () => {
      console.log(xmlHttp.response);
      //data = JSON.parse(xmlHttp.response);
      loadBlock(xmlHttp.response);
    }
    xmlHttp.send()
  }
}

function loadBlock(data) {
  for (var key in data) {
    if (document.getElementById(key)){
      document.getElementById(key).value = Number(data[key])/100000;
    }
  }
}


/*
setInterval(WebRefresh,1000);

function WebRefresh(){

  var xml=createXmlHttpObject();
  xml.responseType = "text";

  xml.open('GET', '/Refresh', true);

  xml.onreadystatechange = function() {
    if (xml.readyState == 4 && xml.status == 200) {
      document.getElementById("test").value = xmlHttp.responseText;
      //console.log(xml.responseText);
    }
  }
  xml.send();
}
*/

// ------------------------------------------------------------------------------------------------------


function LedOn(){
  server = "/Led?status=1";
  send_request(server);
}

function LedOff(){
  server = "/Led?status=0";
  send_request(server);
}

function saveZero(number, zero){
  server = "/Zero?number=" + number + "=" + zero*100000;
  send_request(server);
}

function saveSens(number, sens){
  server = "/Sens?number=" + number + "=" + sens*100000;
  send_request(server);
}

function send_request(server){
  request = new XMLHttpRequest();
  request.open("GET", server, true);
  request.send();
}

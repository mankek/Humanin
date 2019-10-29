function openTab(div_name, elmnt){
    tabcontent = document.getElementsByClassName("content")
    for (var i=0; i < tabcontent.length; i++){
        tabcontent[i].style.display = "none";
    }

    tabs = document.getElementsByClassName("menu")
    for (var s=0; s < tabs.length; s++){
        tabs[s].style.backgroundColor = "";
    }


    document.getElementById(div_name).style.display = "block";
    elmnt.style.backgroundColor = "#ecf8f2";
}

document.getElementById("sets").click();
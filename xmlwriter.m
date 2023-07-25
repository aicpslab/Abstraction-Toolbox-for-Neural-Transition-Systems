import matlab.io.xml.dom.*
docNode = Document("nta");
docRootNode = getDocumentElement(docNode);
weekdays = ["Mon" "Tue" "Wed" "Thu" "Fri"];
weekdaysElement = createElement(docNode,"weekdays");
for i=1:5
    dayElement = createElement(docNode,"day");    
    appendChild(dayElement,createTextNode(docNode,weekdays(i)));
    appendChild(weekdaysElement,dayElement);
end
appendChild(docRootNode,weekdaysElement);
xmlFileName = "example.xml";
writer = matlab.io.xml.dom.DOMWriter;
writer.Configuration.FormatPrettyPrint = true;
writeToFile(writer,docNode,xmlFileName);
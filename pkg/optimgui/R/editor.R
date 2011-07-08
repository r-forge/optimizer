# Operator for convenience
"%+%" = function(x, y) sprintf("%s%s", x, y);

# The "EditorTab" class to describe a tab in the editor
setClass("EditorTab", representation(name = "character",
                                     title = "gLabel",
				                     label = "gTextBox",
				                     rcode = "gTextBox",
                                     output = "gTextBox",
                                     box = "gGroup"));
# The "constructors" of "EditorTab"
EditorTabNew = function(name, title.str, label.str, rcode.str)
{
    # Box to pack widgets
    tabBox = ggroup(horizontal = FALSE, spacing = 15);
    # Widget to show title
    title.str = if(is.null(title.str)) "Edit title here" else title.str;
    titleLabel = glabel(title.str, container = tabBox);
    font(titleLabel) = c(size = "xx-large");
    # Widget to show notes
    flag = is.null(label.str) | !length(label.str);
    label.str = if(flag) "Edit notes here." else label.str;
    docLabel = gtextbox(label.str, container = tabBox,
        			    font.attr = c(family = "sans", size = 11));
    visible(docLabel) = !flag;
    attr(docLabel, "visible") = !flag;
    # Widget to display R code
    flag = is.null(rcode.str);
    rcode.str = if(flag) "# Edit R code here." else rcode.str;
    codeText = gtextbox(rcode.str, container = tabBox,
        			    font.attr = c(family = "monospace", size = 11),
						frame = TRUE);
    visible(codeText) = !flag;
    attr(codeText, "visible") = !flag;
    # Widget to display output
    outputText = gtextbox("", container = tabBox,
        			      font.attr = c(family = "monospace", size = 11),
						  frame = TRUE, frameLabel = "Output");
    visible(outputText) = FALSE;
    
    val = new("EditorTab", name = name, title = titleLabel, label = docLabel,
              rcode = codeText, output = outputText, box = tabBox);
    return(val);
}
EditorTabNewFromXMLNode = function(tabNode)
{
    name = xmlAttrs(tabNode)["tabname"];
    names(name) = NULL;
    title.content = tabNode$children["title"][[1]]$children$text$value;
    label.content = tabNode$children["label"][[1]]$children$text$value;
    rcode.content = tabNode$children["rcode"][[1]]$children$text$value;
    label.content = gsub(" *\n *", "\n", label.content);
    label.content = gsub("\n\n", "#LiNeBrEaK#", label.content);
    label.content = gsub("\n", " ", label.content);
    label.content = gsub("#LiNeBrEaK#", "\n\n", label.content);
    val = EditorTabNew(name, title.content, label.content, rcode.content);
    return(val);
}

# Member functions
show.EditorTab = function(object)
{
    cat("An object of class \"EditorTab\":\n\n");
    cat("@Slot name:\n");
    print(object@name);
    cat("\n");
    cat("@Slot title:\n");
    cat("An object of class \"gLabel\"\n\n");
    cat("@Slot label:\n");
    cat("An object of class \"gTextBox\"\n\n");
    cat("@Slot rcode:\n");
    cat("An object of class \"gTextBox\"\n\n");
    cat("@Slot output:\n");
    cat("An object of class \"gTextBox\"\n");
}
setMethod("show", "EditorTab", show.EditorTab);
# End of definition of "EditorTab"





# The "Editor" class to describe the editor
setClass("Editor", representation(noteBook = "gNotebook",
                                  currentFile = "character",
                                  tabsList = "list",
                                  catalog = "data.frame"));
# The "constructor" of "Editor"
parseRop = function(filePath)
{
    RopFile = readLines(filePath);
    RopFile = RopFile[RopFile != ""];
    index.XML = grep("^#@", RopFile);
    index.rcode = (1:length(RopFile))[-index.XML];
    RopFile[index.XML] = gsub("^#@ ?", "", RopFile[index.XML]);
    RopFile[index.rcode] = gsub("<", "&lt;", RopFile[index.rcode]);
    RopFile[index.rcode] = gsub(">", "&gt;", RopFile[index.rcode]);
    RopTree = xmlTreeParse(RopFile, asText = TRUE);
    return(RopTree);
}
EditorNew = function(...)
{
    noteBook = gnotebook(closebuttons = FALSE, expand = TRUE, ...);
    noteBook@widget@widget$modifyBg(GtkStateType["normal"],
                                    gdkColorParse("white")$color);
    RopPath = system.file("resources", "Rop", package = "optimgui");
    RopFiles = list.files(RopPath, full.names = TRUE);
    getCatalog = function(RopFile)
    {
        d = parseRop(RopFile);
        catalog = xmlElementsByTagName(xmlRoot(d), "catalog")[[1]];
        description = xmlSApply(catalog, xmlAttrs);
        catalog = description[2, ];
        names(catalog) = description[1, ];
        return(c(FileName = basename(RopFile), catalog));
    }
    catalog = t(sapply(RopFiles, getCatalog, simplify = TRUE));
    rownames(catalog) = NULL;
    catalog = as.data.frame(catalog);
    val = new("Editor", noteBook = noteBook, currentFile = "", tabsList = list(),
              catalog = catalog);
    return(val);
}

# Member functions
clearAll.Editor = function(obj)
{
    visible(obj@noteBook) = FALSE;
    while(dispose(obj@noteBook)){}
    visible(obj@noteBook) = TRUE;
    obj@currentFile = character(0);
    obj@tabsList =  list();
    invisible(obj);
}
setGeneric("clearAll", function(obj, ...) standardGeneric("clearAll"));
setMethod("clearAll", "Editor", clearAll.Editor);

showWelcomePage.Editor = function(obj, param, ...)
{
    welcomePage = ggroup(horizontal = FALSE, expand = TRUE);
	welcomeLabel = glabel("Create a new project or open an existing Rop file.",
                          container = welcomePage);
	font(welcomeLabel) = c(size = "x-large");
	tmpBox1 = ggroup(container = welcomePage);
	tmpBox2 = ggroup(container = welcomePage);
	tmpBox3 = ggroup(container = welcomePage);
	gimage(system.file("resources", "images",  "wizard.png", package = "optimgui"), container = tmpBox1);
	wizardLabel = glinklabel("Create a new project by wizard");
    add(tmpBox1, wizardLabel);
	gimage(system.file("resources", "images",  "template.png", package = "optimgui"), container = tmpBox2);
	newLabel = glinklabel("Create a new project with template");
	add(tmpBox2, newLabel);
    gimage(system.file("resources", "images",  "open.png", package = "optimgui"), container = tmpBox3);
	openLabel = glinklabel("Open an existing Rop file");
	add(tmpBox3, openLabel);
    val = list(welcomePage = welcomePage, wizardLabel = wizardLabel,
               newLabel = newLabel, openLabel = openLabel);
	welcomePageAddEvent(val, param);
    add(obj@noteBook, welcomePage, label = "Welcome");
	invisible(obj);
}
setGeneric("showWelcomePage", function(obj, ...) standardGeneric("showWelcomePage"));
setMethod("showWelcomePage", "Editor", showWelcomePage.Editor);

showCatalogPage.Editor = function(obj, param, ...)
{
    catalogPage = ggroup(horizontal = FALSE, expand = TRUE);
	RopTable = gtable(obj@catalog, container = catalogPage, expand = TRUE);
	gseparator(container = catalogPage);
	tmpBox4 = ggroup(container = catalogPage);
	openRopButton = gbutton("OK", container = tmpBox4);
	val = list(catalogPage = catalogPage, RopTable = RopTable,
               openRopButton = openRopButton);
	catalogPageAddEvent(val, param);
    add(obj@noteBook, catalogPage, label = "Catalog");
    invisible(obj);
}
setGeneric("showCatalogPage", function(obj, ...) standardGeneric("showCatalogPage"));
setMethod("showCatalogPage", "Editor", showCatalogPage.Editor);

loadRopFile.Editor = function(obj, filePath)
{
    path = filePath;
	Encoding(path) = "UTF-8";
	RopTree = parseRop(path);
    obj@currentFile = path;
    tabNodes = xmlElementsByTagName(xmlRoot(RopTree), "tab");
    obj@tabsList = lapply(tabNodes, EditorTabNewFromXMLNode);
	names(obj@tabsList) = sapply(tabNodes, xmlAttrs);
    invisible(obj);
}
setGeneric("loadRopFile", function(obj, filePath, ...) standardGeneric("loadRopFile"));
setMethod("loadRopFile", "Editor", loadRopFile.Editor);

saveRopFile.Editor = function(obj, filePath)
{
    Rop = xmlTree("Roptimgui");
	for(tab in obj@tabsList)
	{
		Rop$addNode("tab", attrs = c(tabname = tab@name), close = FALSE);
        Rop$addNode("title", svalue(tab@title));
        if(tab@label@visible)
        {
            labelText = svalue(tab@label);
            labelText = strwrap(labelText, 70);
            labelText = "\n" %+% paste("    ", labelText, sep = "",
                                       collapse = "\n") %+% "\n    ";
            Rop$addNode("label", labelText);
        }
        if(tab@rcode@visible)
        {
            codeText = svalue(tab@rcode);
            codeText = gsub("\n", "\n##Rcode##", codeText);
            codeText = "\n##Rcode##" %+% codeText %+% "\n    ";
            Rop$addNode("rcode", codeText);
        }
		Rop$closeTag();
	}
    z = textConnection("buffer", open = "w");
    sink(z);
	cat(saveXML(Rop$value()));
    sink();
    close(z);
    buffer = paste("#@ ", buffer, sep = "");
    buffer = gsub("#@ ##Rcode##", "", buffer);
	buffer = gsub("&lt;", "<", buffer);
	buffer = gsub("&gt;", ">", buffer);
    writeLines(buffer, filePath);
    invisible(NULL);
}
setGeneric("saveRopFile", function(obj, filePath, ...) standardGeneric("saveRopFile"));
setMethod("saveRopFile", "Editor", saveRopFile.Editor);

buildWidgets.Editor = function(obj)
{
    visible(obj@noteBook) = FALSE;
    while(dispose(obj@noteBook)){}
    if(length(obj@tabsList))
    {
        for(tab in obj@tabsList) add(obj@noteBook, tab@box, label = tab@name);
    }
    visible(obj@noteBook) = TRUE;
    invisible(obj);
}
setGeneric("buildWidgets", function(obj, ...) standardGeneric("buildWidgets"));
setMethod("buildWidgets", "Editor", buildWidgets.Editor);

show.Editor = function(object)
{
    cat("An object of class \"Editor\":\n\n");
    cat("@Slot noteBook:\n");
    cat("An object of class \"gNotebook\"\n\n");
    cat("@Slot currentFile:\n");
    print(object@currentFile);
    cat("\n");
    cat("@Slot tabsList:\n");
    cat(sprintf("A list of length %s\n", length(object@tabsList)));
}
setMethod("show", "Editor", show.Editor);
# End of definition of "Editor"

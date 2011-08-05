# Operator for convenience
"%+%" = function(x, y) sprintf("%s%s", x, y);

#################################
#                               #
#   The "CatalogEntry" class    #
#                               #
#################################

CatalogEntry = setRefClass("CatalogEntry",
                           fields = list(tabname = "character",
                                         model = "RGtkDataFrame",
                                         box = "gGroup"));
# The "constructor" of "CatalogEntry"
CatalogEntryNew = function(tabname, entries)
{
    # Box to pack widgets
    tabBox = ggroup(horizontal = FALSE, spacing = 0);
    # Title
    addSpace(tabBox, 10, horizontal = FALSE);
    titleLabel = glabel("Catalog entries", container = tabBox);
    font(titleLabel) = c(size = "xx-large");
    addSpace(tabBox, 20, horizontal = FALSE);
    # Align box
    alignBox = gtkAlignmentNew(0.5, 0, 1, 1);
    alignBox$setPadding(0, 10, 20, 20);
    # Frame
    frame = gtkFrameNew();
    sw = gtkScrolledWindowNew();
    sw$setPolicy(GtkPolicyType["never"], GtkPolicyType["automatic"]);
    # Table
    textview = gtkTreeViewNew();
    dat = cbind(delete = FALSE, entries);
    model = rGtkDataFrame(dat);
    textview$setModel(model);
    # First column -- checkbox
    renderer = gtkCellRendererToggleNew();
    gSignalConnect(renderer, "toggled", onToggleCell, list(model = model));
    column = gtkTreeViewColumnNewWithAttributes("", renderer, active = 0);
    column$setSizing("fixed");
    column$setFixedWidth(50);
    textview$appendColumn(column);
    # Second column -- Entry Name
    renderer = gtkCellRendererTextNew();
    renderer$set(editable = TRUE);
    gSignalConnect(renderer, "edited", onEditCell, list(column = 2, model = model));
    column = gtkTreeViewColumnNewWithAttributes("Entry Name", renderer, text = 1);
    column$setSizing("fixed");
    column$setFixedWidth(250);
    textview$appendColumn(column);
    # Third column -- Entry Value
    renderer = gtkCellRendererTextNew();
    renderer$set(editable = TRUE);
    gSignalConnect(renderer, "edited", onEditCell, list(column = 3, model = model));
    column = gtkTreeViewColumnNewWithAttributes("Entry Value", renderer, text = 2);
    column$setSizing("fixed");
    column$setFixedWidth(250);
    textview$appendColumn(column);
    # Add widgets
    tabBox@widget@widget$packStart(alignBox);
    alignBox$add(frame);
    frame$add(sw);
    sw$add(textview);
    # ActionsBox
    actionsBox = ggroup(container = tabBox);
    addSpace(actionsBox, 20);
    selectAllChBox = gcheckbox("Select All", container = actionsBox,
                               anchor = c(-1, 0));
    DeleteEntryButton = gbutton("Delete entries", container = actionsBox);
    AddEntryButton = gbutton("Add entry", container = actionsBox);
    addSpace(tabBox, 30, horizontal = FALSE);
    # Add events
    addHandlerChanged(selectAllChBox, onSelectAll, list(model = model));
    addHandlerClicked(AddEntryButton, onAddEntry, list(model = model));
    addHandlerClicked(DeleteEntryButton, onDeleteEntry, list(model = model));

    val = new("CatalogEntry", tabname = tabname, model = model,
              box = tabBox);
    return(val);
}
# Member functions
getItems.CatalogEntry = function(...)
{
    items = .self$model[, 2];
    return(items);
}
CatalogEntry$methods(getItems = getItems.CatalogEntry);

getValues.CatalogEntry = function(...)
{
    values = .self$model[, 3];
    return(values);
}
CatalogEntry$methods(getValues = getValues.CatalogEntry);

haveEntry.CatalogEntry = function(...)
{
    return(nrow(.self$model) > 0);
}
CatalogEntry$methods(haveEntry = haveEntry.CatalogEntry);


#################################
#                               #
#     The "EditorTab" class     #
#                               #
#################################

# The "EditorTab" class to describe a tab in the editor
EditorTab = setRefClass("EditorTab", fields = list(name = "character",
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
    # Align box
    alignBox = gtkAlignmentNew(0.5, 0, 1, 1);
	alignBox$setPadding(10, 20, 0, 2);
    tabBox@widget@widget$packStart(alignBox);
    # Scroll window
    scroll = gtkScrolledWindowNew();
    scroll$setPolicy(GtkPolicyType["never"], GtkPolicyType["automatic"]);
    alignBox$add(scroll);
    # Scroll box
    scrollBox = ggroup(horizontal = FALSE, spacing = 15);
    scroll$addWithViewport(scrollBox@widget@widget);
    viewport = scroll$getChildren()[[1]];
    viewport$setShadowType(GtkShadowType["none"]);
    viewport$modifyBg(GtkStateType["normal"],
                      gdkColorParse("white")$color);
    # Widget to show title
    title.str = if(is.null(title.str)) "Edit title here" else title.str;
    titleLabel = glabel(title.str, editable = TRUE, container = scrollBox);
    font(titleLabel) = c(size = "xx-large");
    # Widget to show notes
    flag = is.null(label.str) | !length(label.str);
    label.str = if(flag) "Edit notes here." else label.str;
    docLabel = gtextbox(label.str, container = scrollBox,
        			    font.attr = c(family = "sans", size = 11));
    visible(docLabel) = !flag;
    # Widget to display R code
    flag = is.null(rcode.str);
    rcode.str = if(flag) "# Edit R code here." else rcode.str;    
    codeText = gtextbox(rcode.str, container = scrollBox,
        			    font.attr = c(family = "monospace", size = 11),
						frame = TRUE);
    visible(codeText) = !flag;
    # Widget to display output
    outputText = gtextbox("", container = scrollBox,
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
    cat("name:\n");
    cat(object$name, "\n\n");
    cat("title:\n");
    cat(svalue(object$title), "\n\n");
    cat("label:\n");
    cat("An object of class \"gTextBox\"\n\n");
    cat("rcode:\n");
    cat("An object of class \"gTextBox\"\n\n");
    cat("output:\n");
    cat("An object of class \"gTextBox\"\n");
}
setMethod("show", "EditorTab", show.EditorTab);
# End of definition of "EditorTab"


#################################
#                               #
#      The "Editor" class       #
#                               #
#################################

# The "Editor" class to describe the editor
Editor = setRefClass("Editor", fields = list(noteBook = "gNotebook",
                                             currentFile = "character",
                                             CatalogEntry = "CatalogEntry",
                                             tabsList = "list",
                                             catalogSet = "data.frame"));
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
    RopFiles = list.files(RopPath, "*\\.[Rr][Oo][Pp]", full.names = TRUE);
    getCatalog = function(RopFile)
    {
        d = parseRop(RopFile);
        catalog = xmlElementsByTagName(xmlRoot(d), "catalog")[[1]];
        description = xmlSApply(catalog, xmlAttrs);
        catalog = description[2, ];
        names(catalog) = description[1, ];
        return(c(FileName = basename(RopFile), catalog));
    }
    catalogSet = lapply(RopFiles, getCatalog);
    itemnames = unique(unlist(lapply(catalogSet, names)));
    tmp = matrix("", length(catalogSet), length(itemnames));
    colnames(tmp) = itemnames;
    for(i in 1:length(catalogSet))
    {
        tmp[i, names(catalogSet[[i]])] = catalogSet[[i]];
    }                               
    catalogSet = as.data.frame(tmp);
    val = new("Editor", noteBook = noteBook, currentFile = "",
              CatalogEntry = new("CatalogEntry"), tabsList = list(),
              catalogSet = catalogSet);
    return(val);
}

# Member functions
clearAll.Editor = function(...)
{
    visible(.self$noteBook) = FALSE;
    while(dispose(.self$noteBook)){}
    visible(.self$noteBook) = TRUE;
    .self$currentFile = character(0);
    .self$CatalogEntry = new("CatalogEntry");
    .self$tabsList = list();
    invisible(.self);
}
Editor$methods(clearAll = clearAll.Editor);


showWelcomePage.Editor = function(param, ...)
{
    welcomePage = ggroup(horizontal = FALSE, expand = TRUE);
    addSpace(welcomePage, 10, horizontal = FALSE);
	welcomeLabel = glabel("Create a new project or open an existing Rop file",
                          container = welcomePage);
	font(welcomeLabel) = c(size = "x-large");
    addSpace(welcomePage, 10, horizontal = FALSE);
	tmpBox1 = ggroup(container = welcomePage);
	tmpBox2 = ggroup(container = welcomePage);
    addSpace(tmpBox1, 10);
	gimage(system.file("resources", "images",  "template.png", package = "optimgui"), container = tmpBox1);
	newLabel = glinklabel("Catalog of existing problems and templates");
	add(tmpBox1, newLabel);
    addSpace(tmpBox2, 10);
    gimage(system.file("resources", "images",  "open.png", package = "optimgui"), container = tmpBox2);
	openLabel = glinklabel("Open an existing Rop file");
	add(tmpBox2, openLabel);
    val = list(welcomePage = welcomePage,
               newLabel = newLabel, openLabel = openLabel);
	welcomePageAddEvent(val, param);
    add(.self$noteBook, welcomePage, label = "Welcome");
	invisible(.self);
}
Editor$methods(showWelcomePage = showWelcomePage.Editor);


showCatalogPage.Editor = function(param, ...)
{
    catalogPage = ggroup(horizontal = FALSE, expand = TRUE);
	RopTable = gtable(.self$catalogSet, container = catalogPage, expand = TRUE);
	gseparator(container = catalogPage);
	tmpBox4 = ggroup(container = catalogPage);
	openRopButton = gbutton("  OK  ", container = tmpBox4);
	val = list(catalogPage = catalogPage, RopTable = RopTable,
               openRopButton = openRopButton);
	catalogPageAddEvent(val, param);
    add(.self$noteBook, catalogPage, label = "Catalog");
    invisible(.self);
}
Editor$methods(showCatalogPage = showCatalogPage.Editor);


loadRopFile.Editor = function(filePath, ...)
{
    path = filePath;
	Encoding(path) = "UTF-8";
	RopTree = parseRop(path);
    .self$currentFile = path;
    catalogNode = xmlElementsByTagName(xmlRoot(RopTree), "catalog");
    if(length(catalogNode))
    {
        catalogNode = catalogNode[[1]];
        descrip = xmlElementsByTagName(catalogNode, "description");
        entries = sapply(descrip, xmlAttrs);
        colnames(entries) = NULL;
        entries = as.data.frame(t(entries), stringsAsFactors = FALSE);
        .self$CatalogEntry = CatalogEntryNew("Catalog", entries);
    }else{
        catalog = .self$catalogSet;
        entries = colnames(catalog)[-1];
        entries = data.frame(itemname = entries, value = "",
                             stringsAsFactors = FALSE);
        .self$CatalogEntry = CatalogEntryNew("Catalog", entries);
    }
    tabNodes = xmlElementsByTagName(xmlRoot(RopTree), "tab");
    if(length(tabNodes))
    {
        .self$tabsList = lapply(tabNodes, EditorTabNewFromXMLNode);
        names(.self$tabsList) = sapply(tabNodes, xmlAttrs);  
    }
    invisible(path);
}
Editor$methods(loadRopFile = loadRopFile.Editor);


saveRopFile.Editor = function(filePath, ...)
{
    Rop = xmlTree("Roptimgui");
    if(.self$CatalogEntry$haveEntry())
    {
        Rop$addNode("catalog", close = FALSE);
        items = .self$CatalogEntry$getItems();
        values = .self$CatalogEntry$getValues();
        for(i in 1:length(items))
        {
            Rop$addNode("description", attrs = c(itemname = items[i],
                                                 value = values[i]));
        }
        Rop$closeTag();
    }
	for(tab in .self$tabsList)
	{
		Rop$addNode("tab", attrs = c(tabname = tab$name), close = FALSE);
        Rop$addNode("title", svalue(tab$title));
        if(visible(tab$label))
        {
            labelText = svalue(tab$label);
            labelText = strwrap(labelText, 70);
            labelText = "\n" %+% paste("    ", labelText, sep = "",
                                       collapse = "\n") %+% "\n    ";
            Rop$addNode("label", labelText);
        }
        if(visible(tab$rcode))
        {
            codeText = svalue(tab$rcode);
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
    invisible(filePath);
}
Editor$methods(saveRopFile = saveRopFile.Editor);


showCurrentCatalog.Editor = function(...)
{
    add(.self$noteBook, .self$CatalogEntry$box, label = .self$CatalogEntry$tabname);
    invisible(.self);
}
Editor$methods(showCurrentCatalog = showCurrentCatalog.Editor);


buildWidgets.Editor = function(...)
{
    visible(.self$noteBook) = FALSE;
    while(dispose(.self$noteBook)){}
    .self$showCurrentCatalog();
    if(length(.self$tabsList))
    {
        for(tab in .self$tabsList) add(.self$noteBook, tab$box, label = tab$name);
    }
    visible(.self$noteBook) = TRUE;
    invisible(.self);
}
Editor$methods(buildWidgets = buildWidgets.Editor);


# Insert an EditorTab at a given position
# CatalogEntry is on position 0
insertTab.Editor = function(tab, after, ...)
{
    notebook = .self$noteBook@widget@widget;
    notebook$insertPage(tab$box@widget@widget, gtkLabelNew(tab$name), after + 1);
    newTab = list(tab);
    names(newTab) = tab$name;
    .self$tabsList = append(.self$tabsList, newTab, after = after);
    invisible(.self);
}
Editor$methods(insertTab = insertTab.Editor);


getTabByName.Editor = function(tabname, ...)
{
    return(.self$tabsList[tabname][[1]]);
}
Editor$methods(getTabByName = getTabByName.Editor);


reportTest.Editor = function(...)
{
    cat("========== Test Report ==========\n");
    objFunTab = .self$getTabByName("Objective");
    objFun = analyzeAssignment(svalue(objFunTab$rcode));
    objFunName = objFun$parName;
    objFunBody = objFun$parVal;
    cat("Objective function:", objFunName, "\n");
    if(!("Gradient" %in% names(.self$tabsList)))
    {
        grFunName = "Not available";
    }else{
        grFunTab = .self$getTabByName("Gradient");
        if(!visible(grFunTab$rcode))
        {
            grFunName = "Not available";
        }else grFunName = analyzeAssignment(svalue(grFunTab$rcode))$parName;
    }
    if(grFunName == "") grFunName = "Not available";
    cat("Gradient function:", grFunName, "\n");
    parTab = .self$getTabByName("Parameters");
    parReport = reportPar(svalue(parTab$rcode));
    constrType = parReport$constrType;
    bound = parReport$withinBound;
    cat("Initial parameters value:", parReport$assignment, "\n");
    cat("Initial function evaluation:", objFunBody(parReport$parVal), "\n");
    cat("Constraints:\n");
    for(i in 1:length(parReport$constrType))
    {
        cat("    ", i, ". ", constrType[i], " constraint", sep = "");
        if(!bound[i]) cat(", initial value doesn't satisfy this constraint!\n") else cat(".\n");
    }
    cat("=================================\n\n");
    return(NULL);
}
Editor$methods(reportTest = reportTest.Editor);


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

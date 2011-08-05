#################################
#                               #
#   Events for "CatalogEntry"   #
#                               #
#################################

# Add entry -- click the button to add a new entry
# h$action = list(model) for onAddEntry(), onDeleteEntry() and onSelectAll()
onAddEntry = function(h, ...)
{
    newEntry = data.frame(delete = FALSE,
                          itemname = "Click here to edit",
                          value = "",
                          stringsAsFactors = FALSE);
    h$action$model$appendRows(newEntry);
    return(NULL);
}
# Delete entries -- click the button to delete selected entries
onDeleteEntry = function(h, ...)
{
    newDat = as.data.frame(h$action$model);
    if(!nrow(newDat)) return(NULL);
    newDat = newDat[!newDat[, 1], ];
    h$action$model$setFrame(newDat);
    return(NULL);
}
# Select all -- select/deselect all entries
onSelectAll = function(h, ...)
{
    model = h$action$model;
    if(!nrow(model)) return(NULL);
    model[, 1] = svalue(h$obj);
    return(NULL);
}

# Toggle checkbox -- Edit cell in column 1
# param = list(model)
onToggleCell = function(cell, path, param)
{
    index = as.integer(path) + 1;
    param$model[index, 1] = !param$model[index, 1];
    return(NULL);
}
# Edit cell -- Edit cell in column 2 and column 3
# param = list(column, model)
onEditCell = function(cell, path, new.text, param)
{
    index = as.integer(path) + 1;
    param$model[index, param$column] = new.text;
    return(NULL);
}


#################################
#                               #
#   Events for menu items       #
#                               #
#################################

# h$action = param = list(mainEnvir, mainWin, buttonGroup, mSave) for all the
# handler functions in this section, except 

# Exit -- exit optimgui
onExit = function(h, ...) dispose(h$action$mainWin);
# Open catalog page -- click the link label to open catalog page
onOpenCatalog = function(h, ...)
{
    editor = h$action$editor;
    editor$showCatalogPage(h$action);
	return(NULL);
}
# Default event -- temporary handler for events that are not implemented
onDefaultEvent = function(h, ...) gmessage("Not implemented yet. :(", "Message");
# Open Rop file -- not a handler, but a function called by some handlers
openRopFile = function(filePath, param)
{
    editor = param$editor;
    visible(param$buttonGroup) = FALSE;
    path = filePath;
	Encoding(path) = "UTF-8";
    editor$loadRopFile(path);
    editor$buildWidgets();
	visible(param$buttonGroup) = TRUE;
    enabled(param$mSave) = TRUE;
    svalue(param$mainWin) = "optimgui [" %+% path %+% "]";
	return(NULL);
}
# Open built-in Rop file -- open the specified template Rop file
onOpenBuiltinRopFile = function(h, ...)
{
	fileName = as.character(svalue(h$action$catalogPage$RopTable));
	filePath = system.file("resources", "Rop", fileName, package = "optimgui");
	if(!length(filePath)) return(NULL);
    openRopFile(filePath, h$action$param);
    h$action$param$editor$currentFile = character(0);
    svalue(h$action$param$mainWin) = "optimgui [New File]";
	return(NULL);
}
# Open Rop file -- open an Rop file using dialog
onOpenRopFile = function(h, ...)
{
	filePath = gfile(text = "Choose an Rop file...", type = "open",
			         filter = list("Rop files (*.Rop)" = list(patterns = "*.Rop")),
			         handler = function(h,...) return(NULL));
	if(is.na(filePath)) return(NULL);
	openRopFile(filePath, h$action);
	return(NULL);
}
# Save Rop file -- save the Rop file using dialog
onSaveRopFile = function(h, ...)
{
    editor = h$action$editor;
    filePath = if(length(editor$currentFile))
    {
        editor$currentFile;
    }else{
        gfile(text = "Save an Rop file...", type = "save",
              filter = list("Rop files (*.Rop)" = list(patterns = "*.Rop")),
              handler = function(h,...) return(NULL));
    }
	if(is.na(filePath)) return(NULL);
    Encoding(filePath) = "UTF-8";
    filePath = gsub("\\.[Rr][Oo][Pp]$", "", filePath);
	editor$saveRopFile(filePath %+% ".Rop");
    svalue(h$action$mainWin) = "optimgui [" %+% filePath %+% ".Rop]";
	return(NULL);
}
# Close Rop file -- close the Rop file and return to the welcome page
onCloseRopFile = function(h, ...)
{
    editor = h$action$editor;
	visible(h$action$buttonGroup) = FALSE;
    editor$clearAll();
    editor$showWelcomePage(h$action);
    enabled(h$action$mSave) = FALSE;
    svalue(h$action$mainWin) = "optimgui";
	return(NULL);
}


#################################
#                               #
#   Events for notebook widget  #
#                               #
#################################

# Focus tab event, to set the value of "noteCheckBox" and "codeCheckBox"
onFocusTab = function(h, ...)
{
    editor = h$action$editor;
    if(!length(editor$tabsList)) return(NULL);
    currentPos = svalue(editor$noteBook) - 1;
    if(currentPos < 1) return(NULL);
    currentTab = editor$tabsList[[currentPos]];
    svalue(h$action$noteCheckBox) = visible(currentTab$label);
    svalue(h$action$codeCheckBox) = visible(currentTab$rcode);
    return(NULL);
}

# Add tab event
onAddTab = function(h, ...)
{
     editor = h$action$editor;
     currentPos = svalue(editor$noteBook) - 1;
     tabName = ginput("Please input the name of the tab:");
     if(is.na(tabName)) return(NULL);
     tab = EditorTabNew(tabName, NULL, NULL, NULL);
     visible(tab$label) = TRUE;
     visible(tab$rcode) = TRUE;
     editor$insertTab(tab, currentPos);
     svalue(editor$noteBook) = currentPos + 2;
     return(NULL);
}
# Delete tab event
onDeleteTab = function(h, ...)
{
    editor = h$action$editor;
    currentPos = svalue(editor$noteBook) - 1;
    if(currentPos < 1)
    {
        gmessage("This tab could not be closed!");
        return(NULL);
    }
    tabname = editor$tabsList[[currentPos]]$name;
    if(tabname %in% c("Objective", "Run"))
    {
        gmessage("This tab could not be closed!");
        return(NULL);
    }
    editor$tabsList = editor$tabsList[-currentPos];
    dispose(editor$noteBook);
    return(NULL);
}
# Test function event
onTestFunction = function(h, ...)
{
    editor = h$action$editor;
    z = textConnection("output", open = "w");
    sink(z);
	editor$reportTest();
	sink();
	close(z);
	output = paste(output, collapse = "\n");
	outputWidget = editor$tabsList[["Run"]]$output;
    visible(outputWidget) = TRUE;
	svalue(outputWidget) = paste(svalue(outputWidget), output, "\n", sep = "");
    svalue(editor$noteBook) = which(names(editor$tabsList) == "Run") + 1;
    return(NULL);
}
# Run code event
onRunCode = function(h, ...)
{
    editor = h$action$editor;
    tmpfile = tempfile();
    editor$saveRopFile(tmpfile);
	z = textConnection("output", open = "w");
	sink(z);
    cat("=========== Run report ==========\n");
	source(tmpfile, local = TRUE, echo = FALSE, print.eval = TRUE);
    cat("=================================\n\n");
	sink();
	close(z);
	output = paste(output, collapse = "\n");
	outputWidget = editor$tabsList[["Run"]]$output;
    visible(outputWidget) = TRUE;
	svalue(outputWidget) = paste(svalue(outputWidget), output, "\n", sep = "");
    svalue(editor$noteBook) = which(names(editor$tabsList) == "Run") + 1;
    return(NULL);
}
# Show or hide notes
toggleShowNote = function(h, ...)
{
    editor = h$action$editor;
    currentPos = svalue(editor$noteBook) - 1;
    if(currentPos < 1) return(NULL);
    tab = editor$tabsList[[currentPos]];
    visible(tab$label) = svalue(h$obj);
    return(NULL);
}
# Show or hide code
toggleShowCode = function(h, ...)
{
    editor = h$action$editor;
    currentPos = svalue(editor$noteBook) - 1;
    if(currentPos < 1) return(NULL);
    tab = editor$tabsList[[currentPos]];
    if(tab$name %in% c("Objective", "Run"))
    {
        svalue(h$obj) = TRUE;
        return(NULL);
    }
    visible(tab$rcode) = svalue(h$obj);
    return(NULL);
}

welcomePageAddEvent = function(welcomePage, param)
{
	addHandlerClicked(welcomePage$newLabel, onOpenCatalog, param);
	addHandlerClicked(welcomePage$openLabel, onOpenRopFile, param);
}
catalogPageAddEvent = function(catalogPage, param)
{
    addHandlerClicked(catalogPage$openRopButton, onOpenBuiltinRopFile,
                      action = list(catalogPage = catalogPage, param = param));
    addHandlerDoubleclick(catalogPage$RopTable, onOpenBuiltinRopFile,
                          action = list(catalogPage = catalogPage, param = param));
}


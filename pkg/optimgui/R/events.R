# Passing parameters
# h$action = param = list(mainEnvir = environment(), mainWin = mainWin,
#                         wizardPage = wizardPage, buttonGroup = buttonGroup,
#                         mSave = mSave);

# Events
# Exit
onExit = function(h, ...) dispose(h$action$mainWin);
# Open wizard
onOpenWizard = function(h, ...)
{
    h$action$wizardPage$show();
	return(NULL);
}
# Wizard confirmed
onWizardConfirmed = function(widget, param)
{
	widget$hide();
	fileName = "shobbs.Rop";
	filePath = system.file("resources", "Rop", fileName, package = "optimgui");
	if(!length(filePath)) return(NULL);
	openRopFile(filePath, param);
	return(FALSE);
}
# Open Catalog Page event
onOpenCatalog = function(h, ...)
{
    editor = get("editor", envir = h$action$mainEnvir);
    showCatalogPage(editor, h$action);
	return(NULL);
}
# Default event
onDefaultEvent = function(h, ...) gmessage("Not implemented yet. :(", "Message");
# Open Rop file
openRopFile = function(filePath, param)
{
    editor = get("editor", envir = param$mainEnvir);
    visible(param$buttonGroup) = FALSE;
    path = filePath;
	Encoding(path) = "UTF-8";
    editor = loadRopFile(editor, path);
    buildWidgets(editor);
    assign("editor", editor, envir = param$mainEnvir);
	visible(param$buttonGroup) = TRUE;
    enabled(param$mSave) = TRUE;
	return(NULL);
}
# Open built-in Rop file
onOpenBuiltinRopFile = function(h, ...)
{
	fileName = as.character(svalue(h$action$catalogPage$RopTable));
	filePath = system.file("resources", "Rop", fileName, package = "optimgui");
	if(!length(filePath)) return(NULL);
	openRopFile(filePath, h$action$param);
	return(NULL);
}
# Open Rop file event
onOpenRopFile = function(h, ...)
{
	filePath = gfile(text = "Choose an Rop file...", type = "open",
			         filter = list("Rop files (*.Rop)" = list(patterns = "*.Rop")),
			         handler = function(h,...) return(NULL));
	if(is.na(filePath)) return(NULL);
	openRopFile(filePath, h$action);
	return(NULL);
}
# Save Rop file event
onSaveRopFile = function(h, ...)
{
    editor = get("editor", envir = h$action$mainEnvir);
    filePath = gfile(text = "Save an Rop file...", type = "save",
			         filter = list("Rop files (*.Rop)" = list(patterns = "*.Rop")),
			         handler = function(h,...) return(NULL));
	if(is.na(filePath)) return(NULL);
    Encoding(filePath) = "UTF-8";
	saveRopFile(editor, filePath %+% ".Rop");
	return(NULL);
}
# Close Rop file event
onCloseRopFile = function(h, ...)
{
    editor = get("editor", envir = h$action$mainEnvir);
	visible(h$action$buttonGroup) = FALSE;
    editor = clearAll(editor);
    showWelcomePage(editor, h$action);
    assign("editor", editor, envir = h$action$mainEnvir);
    enabled(h$action$mSave) = FALSE;
	return(NULL);
}
# Add tab event
onAddTab = function(h, ...)
{
     editor = get("editor", envir = h$action$mainEnvir);
     currentPos = svalue(editor@noteBook);
     tabName = ginput("Please input the name of the tab:");
     if(is.na(tabName)) return(NULL);
     tab = EditorTabNew(tabName, NULL, NULL, NULL);
     visible(tab@label) = TRUE;
     visible(tab@rcode) = TRUE;
     editor = insertTab(editor, tab, currentPos);
     svalue(editor@noteBook) = currentPos + 1;
     assign("editor", editor, envir = h$action$mainEnvir);
     return(NULL);
}
# Delete tab event
onDeleteTab = function(h, ...)
{
    editor = get("editor", envir = h$action$mainEnvir);
    currentPos = svalue(editor@noteBook);
    editor@tabsList = editor@tabsList[-currentPos];
    dispose(editor@noteBook);
    assign("editor", editor, envir = h$action$mainEnvir);
    return(NULL);
}
# Run code event
onRunCode = function(h, ...)
{
    editor = get("editor", envir = h$action$mainEnvir);
	if(length(editor@currentFile))
	{
		z = textConnection("output", open = "w");
		sink(z);
		source(editor@currentFile, echo = FALSE, print.eval = TRUE);
		sink();
		close(z);
		output = paste(output, collapse = "\n");
		outputWidget = editor@tabsList[["Run"]]@output;
        visible(outputWidget) = TRUE;
		svalue(outputWidget) = output;
        svalue(editor@noteBook) = which(names(editor@tabsList) == "Run");
	}
    return(NULL);
}

welcomePageAddEvent = function(welcomePage, param)
{
    addHandlerClicked(welcomePage$wizardLabel, onOpenWizard, param);
	addHandlerClicked(welcomePage$newLabel, onOpenCatalog, param);
	addHandlerClicked(welcomePage$openLabel, onOpenRopFile, param);
}
catalogPageAddEvent = function(catalogPage, param)
{
    addHandlerClicked(catalogPage$openRopButton, onOpenBuiltinRopFile,
                      action = list(catalogPage = catalogPage, param = param));
}


# Events
# Exit
onExit = function(h, ...) dispose(mainWin);
# Open wizard
onOpenWizard = function(h, ...)
{
    wizardPage$show();
	return(NULL);
}
# Wizard confirmed
onWizardConfirmed = function(widget, param)
{
	widget$hide();
	fileName = "shobbs.Rop";
	filePath = system.file("resources", "Rop", fileName, package = "optimgui");
	if(!length(filePath)) return(NULL);
	openRopFile(filePath);
	return(FALSE);
}
# Open Catalog Page event
onOpenCatalog = function(h, ...)
{
    showCatalogPage.Editor(editor);
	return(NULL);
}
# Default event
onDefaultEvent = function(h, ...) gmessage("Not implemented yet. :(", "Message");
# Open built-in Rop file
# Open Rop file
openRopFile = function(filePath)
{
    visible(buttonGroup) = FALSE;
    path = filePath;
	Encoding(path) = "UTF-8";
    editor <<- openRopFile.Editor(editor, path);
    editor <<- buildWidgets.Editor(editor);
	visible(buttonGroup) = TRUE;
    enabled(mSave) = TRUE;
	return(NULL);
}
onOpenBuiltinRopFile = function(h, ...)
{
	fileName = as.character(svalue(h$action$RopTable));
	filePath = system.file("resources", "Rop", fileName, package = "optimgui");
	if(!length(filePath)) return(NULL);
	openRopFile(filePath);
	return(NULL);
}
# Open Rop file event
onOpenRopFile = function(h, ...)
{
	filePath = gfile(text = "Choose an Rop file...", type = "open",
			         filter = list("Rop files (*.Rop)" = list(patterns = "*.Rop")),
			         handler = function(h,...) return(NULL));
	if(is.na(filePath)) return(NULL);
	openRopFile(filePath);
	return(NULL);
}
# Save Rop file event
onSaveRopFile = function(h, ...)
{
    filePath = gfile(text = "Save an Rop file...", type = "save",
			         filter = list("Rop files (*.Rop)" = list(patterns = "*.Rop")),
			         handler = function(h,...) return(NULL));
	if(is.na(filePath)) return(NULL);
    Encoding(filePath) = "UTF-8";
	saveRopFile.Editor(editor, filePath %+% ".Rop");
	return(NULL);
}
# Close Rop file event
onCloseRopFile = function(h, ...)
{
	visible(buttonGroup) = FALSE;
    editor <<- clearAll.Editor(editor);
    editor <<- showWelcomePage.Editor(editor);
    enabled(mSave) = FALSE;
	return(NULL);
}
# Run code event
onRunCode = function(h, ...)
{
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

welcomePageAddEvent = function(welcomePage)
{
    addHandlerClicked(welcomePage$wizardLabel, onOpenWizard);
	addHandlerClicked(welcomePage$newLabel, onOpenCatalog);
	addHandlerClicked(welcomePage$openLabel, onOpenRopFile);
}
catalogPageAddEvent = function(catalogPage)
{
    addHandlerClicked(catalogPage$openRopButton, handler = onOpenBuiltinRopFile,
                   action = catalogPage);
}


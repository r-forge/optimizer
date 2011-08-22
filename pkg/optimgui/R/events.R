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
    enabled(param$mSave) = enabled(param$mSaveAs) =
        enabled(param$mSaveAsT) = TRUE;
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
    h$action$editor$currentFile = filePath;
	return(NULL);
}
# Save Rop file -- save the Rop file using dialog
onSaveRopFile = function(h, ...)
{
    editor = h$action$editor;
    oldwd = setwd(editor$repoPath);
    filePath = if(length(editor$currentFile))
    {
        editor$currentFile;
    }else{
        gfile(text = "Save an Rop file...", type = "save",
              filter = list("Rop files (*.Rop)" = list(patterns = "*.Rop")),
              handler = function(h,...) return(NULL));
    }
    setwd(oldwd);
	if(is.na(filePath)) return(NULL);
    Encoding(filePath) = "UTF-8";
    filePath = gsub("\\.[Rr][Oo][Pp]$", "", filePath) %+% ".Rop";
	editor$saveRopFile(filePath);
    svalue(h$action$mainWin) = "optimgui [" %+% filePath %+% "]";
    editor$currentFile = filePath;
	return(NULL);
}
onSaveAsRopFile = function(h, ...)
{
    editor = h$action$editor;
    oldwd = setwd(editor$repoPath);
    filePath = gfile(text = "Save an Rop file...", type = "save",
                     filter = list("Rop files (*.Rop)" = list(patterns = "*.Rop")),
                     handler = function(h,...) return(NULL));
    setwd(oldwd);
    if(is.na(filePath)) return(NULL);
    Encoding(filePath) = "UTF-8";
    filePath = gsub("\\.[Rr][Oo][Pp]$", "", filePath) %+% ".Rop";
	editor$saveRopFile(filePath);
    svalue(h$action$mainWin) = "optimgui [" %+% filePath %+% "]";
    editor$currentFile = filePath;
	return(NULL);
}
onSaveAsTRopFile = function(h, ...)
{
    editor = h$action$editor;
    fileName = ginput("Please enter the file name:",
                      parent = h$action$mainWin);
    if(is.na(fileName)) return(NULL);
    Encoding(fileName) = "UTF-8";
    fileName = gsub("\\.[Rr][Oo][Pp]$", "", fileName) %+% ".Rop";
    filePath = system.file("resources", "Rop", package = "optimgui");
    filePath = file.path(filePath, fileName);
    editor$saveRopFile(filePath);
    svalue(h$action$mainWin) = "optimgui [" %+% filePath %+% "]";
    editor$currentFile = filePath;
	return(NULL);
}
# Set the user repository to store Rop files
onSetUserRepo = function(h, ...)
{
    editor = h$action$editor;
    oldwd = setwd(editor$repoPath);
    repoPath = gfile(text = "Select a directory to store Rop files", type = "selectdir");
    setwd(oldwd);
    if(is.na(repoPath)) return(NULL);
    Encoding(repoPath) = "UTF-8";
    editor$repoPath = repoPath;
    return(NULL);
}
# Close Rop file -- close the Rop file and return to the welcome page
onCloseRopFile = function(h, ...)
{
    editor = h$action$editor;
	visible(h$action$buttonGroup) = FALSE;
    editor$clearAll();
    editor$showWelcomePage(h$action);
    enabled(h$action$mSave) = enabled(h$action$mSaveAs) =
        enabled(h$action$mSaveAsT) = FALSE;
    svalue(h$action$mainWin) = "optimgui";
	return(NULL);
}
# Open the help manual
onOpenManual = function(h, ...)
{
    helpFile = system.file("doc", "optimgui.pdf", package = "optimgui");
    shell.exec(helpFile);
    return(NULL);
}
# The About dialog
onAboutDialog = function(h, ...)
{
    win = h$action$mainWin@widget@widget;
    comments = packageDescription("optimgui", fields = "Title");
    comments = paste(strwrap(comments, 50), collapse = "\n");
    logo = gdkPixbufNewFromFile(system.file("resources", "images",
                                            "Logo-png.png",
                                            package = "optimgui"));
    gtkShowAboutDialog(win,
                       authors = c("John C. Nash <nashjc@uottawa.ca>",
                                   "Ben Bolker <bbolker@gmail.com>",
                                   "Yixuan Qiu <yixuan.qiu@cos.name>"),
                       comments = comments,
                       logo = logo,
                       "program-name" = "optimgui",
                       version = as.character(packageVersion("optimgui")),
                       title = "About optimgui",
                       website = "https://r-forge.r-project.org/scm/viewvc.php/pkg/optimgui/?root=optimizer",
                       "website-label" = "Website on R-Forge");
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
     tabName = ginput("Please input the name of the tab:",
                      parent = h$action$mainWin);
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
        gmessage("This tab could not be closed!", "Message",
                 parent = h$action$mainWin);
        return(NULL);
    }
    tabname = editor$tabsList[[currentPos]]$name;
    if(tabname %in% c("Objective", "Parameters", "Run"))
    {
        gmessage("This tab could not be closed!", "Message",
                 parent = h$action$mainWin);
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
    error = function(e)
    {
        gmessage("Error occurs when testing function!", "Error", icon = "error",
                 parent = h$action$mainWin);
    }
    sink(z);
	tryCatch(editor$reportTest(), error = error);
    sink();
    close(z);
	output = paste(output, collapse = "\n");
	outputWidget = editor$tabsList[["Run"]]$output;
    visible(outputWidget) = TRUE;
	svalue(outputWidget) = paste(svalue(outputWidget), output, "\n", sep = "");
    svalue(editor$noteBook) = which(names(editor$tabsList) == "Run") + 1;
    return(NULL);
}
# Choose method event
chooseMethod = function(optInfo)
{
    cmd1 = "library(optimx)
ans <- optimx(par = %s, fn = %s, gr = %s,
              lower = %s, upper = %s,
              control = list(trace = 1))";
    cmd2 = "ans <- constrOptim(theta = %s, f = %s, grad = %s,
                   ui = %s, ci = %s,
                   control = list(trace = 1))";
    cmd3 = "library(Rsolnp)
ans <- solnp(pars = %s, fun = %s, ineqfun = %s,
             ineqLB = %s, ineqUB = %s,
             LB = %s, UB = %s)";
    runCode = "";
    grad = gsub("Not available", "NULL", optInfo$grFunName);
    if(!is.null(optInfo$nonlinearConstr2$assign))
    {
        if(!is.null(optInfo$boxConstr$assign))
        {
            runCode = paste(optInfo$boxConstr$assign, collapse = "\n") %+% "\n";
        } else {
            runCode = "LB <- NULL\nUB <- NULL\n";
        }
        ineqLB = rep(-Inf, length(optInfo$nonlinearConstr2$ineqUB));
        runCode = runCode %+% paste(optInfo$nonlinearConstr2$assign, collapse = "\n") %+%
            "\n" %+% sprintf("ineqLB <- %s", deparse(ineqLB)) %+% "\n";
        runCode = runCode %+% sprintf(cmd3, optInfo$parName, optInfo$objFunName,
                                      "ineqFun", "ineqLB", "ineqUB", "LB", "UB");
    } else if(!is.null(optInfo$linearConstr$assign)){
        runCode = paste(optInfo$boxlinearConstr$assign, collapse = "\n") %+% "\n";
        runCode = runCode %+% sprintf(cmd2, optInfo$parName, optInfo$objFunName,
                                      grad, "-linMat", "-linUB");
    } else if(!is.null(optInfo$boxConstr$assign)){
        runCode = paste(optInfo$boxConstr$assign, collapse = "\n") %+% "\n";
        runCode = runCode %+% sprintf(cmd1, optInfo$parName, optInfo$objFunName,
                                      grad, "LB", "UB");
    } else {
        runCode = runCode %+% sprintf(cmd1, optInfo$parName, optInfo$objFunName,
                                      grad, "-Inf", "Inf");
    }
    return(runCode);
}
onChooseMethod = function(h, ...)
{
    editor = h$action$editor;
    if(!length(editor$optInfo))
    {
        gmessage("You should first test the function by clicking the \"Test function\" button!",
                 "Message", parent = h$action$mainWin);
        return(NULL);
    }
    runTab = editor$getTabByName("Run");
    svalue(runTab$rcode) = chooseMethod(editor$optInfo);
    return(NULL);
}
# Run code event
onRunCode = function(h, ...)
{
    editor = h$action$editor;
    tmpfile = tempfile();
    editor$saveRopFile(tmpfile);
    error = function(e)
    {
        gmessage("Error occurs when running code!", "Error", icon = "error",
                 parent = h$action$mainWin);
    }
	z = textConnection("output", open = "w");
	sink(z);
    cat("===================== Run Report ====================\n");
	tryCatch(source(tmpfile, local = TRUE, echo = FALSE, print.eval = TRUE),
             error = error);
    cat("=====================================================\n\n");
	sink();
	close(z);
	output = paste(output, collapse = "\n");
	outputWidget = editor$tabsList[["Run"]]$output;
    visible(outputWidget) = TRUE;
	svalue(outputWidget) = paste(svalue(outputWidget), output, "\n", sep = "");
    svalue(editor$noteBook) = which(names(editor$tabsList) == "Run") + 1;
    return(NULL);
}
# Save output event
onSaveOutput = function(h, ...)
{
    editor = h$action$editor;
    runTab = editor$getTabByName("Run");
    filePath = gfile(text = "Save ouput...", type = "save",
                     filter = list("Text files (*.txt)" = list(patterns = "*.txt")),
                     handler = function(h,...) return(NULL));
    if(is.na(filePath)) return(NULL);
    Encoding(filePath) = "UTF-8";
    filePath = gsub("\\.[Tt][Xx][Tt]$", "", filePath) %+% ".txt";
    write(sprintf("optimgui report\n\nDate: %s\n\n", date()), filePath);
    write(svalue(runTab$output), filePath, append = TRUE);
	return(NULL);
}
# Clear output event
onClearOutput = function(h, ...)
{
    editor = h$action$editor;
    runTab = editor$getTabByName("Run");
    svalue(runTab$output) = "";
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
    addHandlerClicked(welcomePage$manualLabel, onOpenManual, param);
}
catalogPageAddEvent = function(catalogPage, param)
{
    addHandlerClicked(catalogPage$openRopButton, onOpenBuiltinRopFile,
                      action = list(catalogPage = catalogPage, param = param));
    addHandlerDoubleclick(catalogPage$RopTable, onOpenBuiltinRopFile,
                          action = list(catalogPage = catalogPage, param = param));
}


# Uncomment this before building package
optimgui = function(){

# Comment this before building package
#library(XML);
#library(RGtk2);
#library(gWidgets);
#library(gWidgetsRGtk2);
#source("aaaWidgets.R");
#source("editor.R");
#source("events.R");
#source("parameters.R");
#source("zzz.R");



# Main Window
mainWin = gwindow("optimgui", visible = FALSE);
size(mainWin) = c(900, 600);
# Horizontal layout box
hGroup = ggroup(container = mainWin, expand = TRUE, spacing = 15);
# The editor
editor = EditorNew(container = hGroup);
# Button box
buttonGroup = ggroup(container = hGroup, horizontal = FALSE);
visible(buttonGroup) = FALSE;
# Buttons
testButton = gbutton("Test function", container = buttonGroup);
chooseButton = gbutton("Choose method", container = buttonGroup);
runButton = gbutton("Run", container = buttonGroup);
addSpace(buttonGroup, 25, horizontal = FALSE);
saveOutputButton = gbutton("Save output", container = buttonGroup);
clearOutputButton = gbutton("Clear output", container = buttonGroup);
addSpace(buttonGroup, 25, horizontal = FALSE);
addButton = gbutton("Add tab", container = buttonGroup);
deleteButton = gbutton("Delete tab", container = buttonGroup);
noteCheckBox = gcheckbox("Show note", checked = TRUE, container = buttonGroup);
codeCheckBox = gcheckbox("Show code", checked = TRUE, container = buttonGroup);

# Add events
# Parameters to be passed to event handlers
# "mSave" must be created before "param" is assigned.
mSave = gaction(label = "Save", handler = onSaveRopFile,
                action = list(mainWin = mainWin, editor = editor), icon = "save");
mSaveAs = gaction(label = "Save As...", handler = onSaveAsRopFile,
                  action = list(mainWin = mainWin, editor = editor), icon = "save-as");
mSaveAsT = gaction(label = "Save As Template...", handler = onSaveAsTRopFile,
                   action = list(mainWin = mainWin, editor = editor));
param = list(mainWin = mainWin, editor = editor,
             buttonGroup = buttonGroup, mSave = mSave, mSaveAs = mSaveAs,
             mSaveAsT = mSaveAsT, noteCheckBox = noteCheckBox,
             codeCheckBox = codeCheckBox);
editor$showWelcomePage(param);
addHandlerFocus(editor$noteBook, onFocusTab, param);
addHandlerClicked(testButton, onTestFunction, param);
addHandlerClicked(chooseButton, onChooseMethod, param);
addHandlerClicked(runButton, onRunCode, param);
addHandlerClicked(saveOutputButton, onSaveOutput, param);
addHandlerClicked(clearOutputButton, onClearOutput, param);
addHandlerClicked(addButton, onAddTab, param);
addHandlerClicked(deleteButton, onDeleteTab, param);
addHandlerChanged(noteCheckBox, toggleShowNote, param);
addHandlerChanged(codeCheckBox, toggleShowCode, param);

# Menu list
mOpen = gaction(label = "Open File...",
                handler = onOpenRopFile, action = param, icon = "open");
mSetDir = gaction(label = "Set User Repository...",
                  handler = onSetUserRepo, action = param)
mClose = gaction(label = "Close File",
                 handler = onCloseRopFile, action = param, icon = "close");
mExit = gaction(label = "Exit optimgui",
                handler = onExit, action = param, icon = "quit");
mHelp = gaction(label = "optimgui Help", handler = onOpenManual, icon = "help");
mAbout = gaction(label = "About", handler = onAboutDialog, action = param,
                 icon = "about");

menuList = list(File = list(open = mOpen, save = mSave, saveAs = mSaveAs,
                            saveAsT = mSaveAsT, mSetDir = mSetDir,
                            close = mClose, exit = mExit),
		        Help = list(help = mHelp, about = mAbout));
# Main menu
mainMenu = gmenu(menuList, container = mainWin);
menuFolder1 = mainMenu@widget@widget$getChildren();
menuFolder2 = lapply(menuFolder1, function(menuItem) menuItem$getSubmenu()$getChildren());
accelGroup = gtkAccelGroupNew();
mainWin@widget@widget$addAccelGroup(accelGroup);
mapply(addKeyAccel, menuFolder2[[1]], list("o", "s", NULL, NULL, NULL, "w", "q"),
       MoreArgs = list(accelGroup = accelGroup));
enabled(mSave) = enabled(mSaveAs) = enabled(mSaveAsT) = FALSE;
# Show main window
visible(mainWin) = TRUE;

# Uncomment this before building package
}

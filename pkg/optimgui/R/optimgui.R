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
addButton = gbutton("Add tab", container = buttonGroup);
deleteButton = gbutton("Delete tab", container = buttonGroup);
runButton = gbutton("Run", container = buttonGroup);
noteCheckBox = gcheckbox("Show note", checked = TRUE, container = buttonGroup);
codeCheckBox = gcheckbox("Show code", checked = TRUE, container = buttonGroup);
# synButton = gbutton("Syntex check", container = buttonGroup);
# anoButton = gbutton("Another", container = buttonGroup);

# Add events
# Parameters to be passed to event handlers
# "mSave" must be created before "param" is assigned.
mSave = gaction(label = "Save", handler = onSaveRopFile,
                action = list(mainWin = mainWin, editor = editor), icon = "save");
param = list(mainWin = mainWin, editor = editor,
             buttonGroup = buttonGroup, mSave = mSave,
             noteCheckBox = noteCheckBox, codeCheckBox = codeCheckBox);
editor$showWelcomePage(param);
addHandlerFocus(editor$noteBook, onFocusTab, param);
addHandlerClicked(addButton, onAddTab, param);
addHandlerClicked(deleteButton, onDeleteTab, param);
addHandlerClicked(runButton, onRunCode, param);
addHandlerChanged(noteCheckBox, toggleShowNote, param);
addHandlerChanged(codeCheckBox, toggleShowCode, param);

# Menu list
mOpen = gaction(label = "Open", handler = onOpenRopFile, action = param, icon = "open");
mClose = gaction(label = "Close", handler = onCloseRopFile, action = param, icon = "close");
mExit = gaction(label = "Exit", handler = onExit, action = param, icon = "quit");
mCopy = gaction(label = "Copy", handler = onDefaultEvent);
mCut = gaction(label = "Cut", handler = onDefaultEvent);
mTool = gaction(label = "Tool", handler = onDefaultEvent);
mHelp = gaction(label = "optimgui Help", handler = onDefaultEvent);
mAbout = gaction(label = "About", handler = onDefaultEvent);

menuList = list(File = list(open = mOpen, save = mSave, close = mClose,
                            exit = mExit),
		        Edit = list(copy = mCopy, cut = mCut),
		        Tools = list(tool = mTool),
		        Help = list(help = mHelp, about = mAbout));
# Main menu
gmenu(menuList, container = mainWin);
enabled(mSave) = FALSE;
# Show main window
visible(mainWin) = TRUE;

# Uncomment this before building package
}
